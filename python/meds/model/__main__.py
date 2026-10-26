# SPDX-License-Identifier: Apache-2.0
"""``python -m meds.model CONFIG`` -- run a config to its end through the Python API, in place of
``meds_main CONFIG``: the same output files and the same lines on stdout. The exit status is 0 when
the run completes and 1 when it stops on a bad state."""
import sys

from ._run import run


def main(argv=None):
    argv = sys.argv[1:] if argv is None else argv
    if len(argv) != 1 or argv[0].startswith("-"):
        print("usage: python -m meds.model CONFIG", file=sys.stderr)
        return 2
    try:
        run(argv[0])
    except RuntimeError as e:
        print(f" ERROR: {e}", flush=True)
        return 1
    return 0


if __name__ == "__main__":
    sys.exit(main())
