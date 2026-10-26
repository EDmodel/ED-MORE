# SPDX-License-Identifier: Apache-2.0
"""The workers that run trials and chains: each is one MEDS process, and a fit runs hundreds at a time.

A run uses as many threads as its length earns (trials.threads_for), and takes that many of its
machine's cores: a long run (a seasonal run, a chain) finishes as soon as the ten-day windows beside
it, and the cores the windows leave idle are not wasted. The long runs of a batch start first.

- LocalWorkers: N cores on this machine (a workstation, or one node).
- QueueWorkers: a directory queue shared by workers on many nodes. The driver writes one task file
  per run; each node's `calibrate_fast.py worker` (one per node, started inside the same Slurm
  allocation) claims tasks by an atomic rename, one at a time and only when its cores have room,
  runs them and writes a done file. One allocation holds the workers for the whole fit, so no run
  waits in the Slurm queue.

Both run a batch of tasks and return {task id: (status, seconds)}, status "ok", "timeout" or
"exit <code>". A task is (id, command, directory, log file, timeout in seconds, threads).
"""
from __future__ import annotations

import json
import os
import socket
import subprocess
import threading
import time
import uuid
from concurrent.futures import ThreadPoolExecutor
from dataclasses import dataclass
from pathlib import Path


@dataclass
class Task:
    id: str
    argv: list
    cwd: str
    log: str
    timeout: float
    threads: int = 1


def run_task(t: Task) -> tuple[str, float]:
    t0 = time.time()
    env = dict(os.environ, OMP_NUM_THREADS=str(t.threads))
    try:
        with open(t.log, "w") as fh:
            res = subprocess.run(t.argv, cwd=t.cwd, stdout=fh, stderr=subprocess.STDOUT,
                                 timeout=t.timeout, env=env)
        status = "ok" if res.returncode == 0 else f"exit {res.returncode}"
    except subprocess.TimeoutExpired:
        status = "timeout"
    return status, time.time() - t0


def longest_first(tasks: list[Task]) -> list[Task]:
    """The order a batch starts in: the most threads (the longest runs) first."""
    return sorted(tasks, key=lambda t: -t.threads)


class Cores:
    """A machine's cores: a run takes as many as it has threads (at most all), and runs take them
    in the order they asked, so a long run is not passed over by the short ones behind it."""

    def __init__(self, n: int):
        self.n = self.free = max(1, int(n))
        self.cond = threading.Condition()
        self.next_ticket = self.serving = 0

    def take(self, threads: int) -> int:
        k = min(max(1, int(threads)), self.n)
        with self.cond:
            ticket, self.next_ticket = self.next_ticket, self.next_ticket + 1
            self.cond.wait_for(lambda: self.serving == ticket and self.free >= k)
            self.free -= k
            self.serving += 1
            self.cond.notify_all()
        return k

    def give(self, k: int):
        with self.cond:
            self.free += k
            self.cond.notify_all()


class LocalWorkers:
    """N cores on this machine, shared by every caller (the chains start together)."""

    def __init__(self, n: int):
        self.capacity = max(1, int(n))           # the most threads one run can use
        self.cores = Cores(self.capacity)
        self.ex = ThreadPoolExecutor(self.capacity)

    def run_one(self, t: Task) -> tuple[str, float]:
        k = self.cores.take(t.threads)
        try:
            return run_task(t)
        finally:
            self.cores.give(k)

    def run(self, tasks: list[Task]) -> dict:
        futures = {t.id: self.ex.submit(self.run_one, t) for t in longest_first(tasks)}
        return {t.id: futures[t.id].result() for t in tasks}

    def close(self):
        self.ex.shutdown(wait=True)


class QueueWorkers:
    """The driver's side of the directory queue at `root` (queue/, claimed/, done/)."""

    capacity = 1 << 30                           # the workers' nodes set the threads they can give

    def __init__(self, root, poll: float = 0.5):
        self.root = Path(root)
        for sub in ("queue", "claimed", "done"):
            (self.root / sub).mkdir(parents=True, exist_ok=True)
        self.poll = poll

    def run(self, tasks: list[Task]) -> dict:
        batch = uuid.uuid4().hex[:8]
        ids = {}
        for t in longest_first(tasks):              # the queue is claimed in name order
            name = f"{batch}-{len(ids):06d}"
            ids[name] = t.id
            tmp = self.root / "queue" / f".{name}.tmp"
            tmp.write_text(json.dumps(t.__dict__))
            tmp.rename(self.root / "queue" / f"{name}.json")
        out = {}
        while len(out) < len(ids):
            for name in list(ids):
                if ids[name] in out:
                    continue
                done = self.root / "done" / f"{name}.json"
                if done.exists():
                    try:
                        rec = json.loads(done.read_text())
                    except json.JSONDecodeError:
                        continue                     # still being written
                    out[ids[name]] = (rec["status"], rec["seconds"])
                    done.unlink()
            if len(out) < len(ids):
                time.sleep(self.poll)
        return out

    def close(self):
        (self.root / "STOP").write_text("stop\n")


def queue_worker(root, slots: int, poll: float = 0.5) -> None:
    """A node's worker: claim tasks from root/queue until root/STOP exists. `slots` is the node's
    cores; a run takes as many as it has threads (at most all). The worker claims one task at a time
    and holds at most one that does not fit yet, so the other nodes get the rest."""
    root = Path(root)
    me = f"{socket.gethostname()}-{os.getpid()}"
    for sub in ("queue", "claimed", "done"):
        (root / sub).mkdir(parents=True, exist_ok=True)
    running, waiting = {}, None                  # name -> (future, cores); a claimed task with no room yet
    with ThreadPoolExecutor(slots) as ex:
        while True:
            for name, (fut, _) in list(running.items()):
                if fut.done():
                    status, secs = fut.result()
                    tmp = root / "done" / f".{name}.tmp"
                    tmp.write_text(json.dumps({"status": status, "seconds": secs, "worker": me}))
                    tmp.rename(root / "done" / f"{name}.json")
                    (root / "claimed" / f"{me}.{name}.json").unlink(missing_ok=True)
                    del running[name]
            free = slots - sum(k for _, k in running.values())
            while free > 0:
                if waiting is None:
                    for q in sorted((root / "queue").glob("*.json")):
                        claim = root / "claimed" / f"{me}.{q.name}"
                        try:
                            q.rename(claim)          # atomic: one worker wins
                        except (FileNotFoundError, OSError):
                            continue
                        waiting = (q.stem, Task(**json.loads(claim.read_text())))
                        break
                    if waiting is None:
                        break                        # the queue is empty
                name, t = waiting
                k = min(max(1, t.threads), slots)
                if k > free:
                    break                            # wait for cores
                running[name] = (ex.submit(run_task, t), k)
                free -= k
                waiting = None
            if ((root / "STOP").exists() and not running and waiting is None
                    and not any((root / "queue").glob("*.json"))):
                return
            time.sleep(poll)


def make_workers(kind: str, n: int, queue_dir=None):
    """"local": n processes on this machine; "queue": the directory queue at queue_dir."""
    if kind == "local":
        return LocalWorkers(n)
    if kind == "queue":
        return QueueWorkers(queue_dir)
    raise ValueError(f"unknown workers '{kind}'")
