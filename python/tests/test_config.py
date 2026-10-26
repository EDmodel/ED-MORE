# SPDX-License-Identifier: Apache-2.0
"""meds.config reads a run's main file and the PFT (plant trait) file it names, changes keys, and
writes a pair of files that read back as set. No compiled library is needed."""
from pathlib import Path

import pytest

try:
    import tomllib
except ModuleNotFoundError:  # py3.10 and older
    import tomli as tomllib

from meds.config import RunConfig, dumps, load_toml, read_record

ROOT = Path(__file__).resolve().parents[2]


def test_load_reads_both_files_and_resolves_inputs():
    cfg = RunConfig.load(ROOT / "meds_config_main.toml")
    assert cfg.pft_path == ROOT / "meds_config_pft.toml"
    assert Path(cfg.get("init.census_file")).is_absolute()
    assert Path(cfg.get("init.census_file")).exists()
    assert cfg.get("output.io_config") == ""          # "no file" is left alone
    assert cfg.n_pft == len(cfg.get("pft.vcmax25", file="pft")) > 1


def test_set_one_pft_keeps_the_others(tmp_path):
    cfg = RunConfig.load(ROOT / "meds_config_main.toml")
    before = list(cfg.get("pft.vcmax25", file="pft"))
    trial = cfg.copy()
    trial.set("pft.vcmax25", 61.5, file="pft", pft=2)
    trial.set("aerodynamics.z0m_ratio", 0.11)
    trial.set("run.end_time", "2001-06-03 00:00:00")
    assert cfg.get("pft.vcmax25", file="pft") == before          # the copy is independent
    main = trial.write(tmp_path / "run")
    m, p = load_toml(main), load_toml(tmp_path / "run" / "pft.toml")
    assert m["init"]["pft_config"] == str(tmp_path / "run" / "pft.toml")
    assert m["aerodynamics"]["z0m_ratio"] == 0.11
    assert p["pft"]["vcmax25"] == before[:1] + [61.5] + before[2:]
    assert RunConfig.load(main).get("pft.vcmax25", file="pft", pft=2) == 61.5


def test_a_missing_trait_needs_a_single_pft():
    cfg = RunConfig(main={}, pft={"pft": {"wood_density": [0.6, 0.7]}})
    with pytest.raises(KeyError):
        cfg.set("pft.leaf_width", 0.05, file="pft", pft=1)
    one = RunConfig(main={}, pft={"pft": {"wood_density": [0.6]}})
    one.set("pft.leaf_width", 0.05, file="pft", pft=1)
    assert one.get("pft.leaf_width", file="pft") == [0.05]



def test_a_key_for_every_pft_reads_and_sets_one_pft():
    cfg = RunConfig(main={}, pft={"pft": {"wood_density": [0.6, 0.7]}, "hydraulics": {"root_beta": 0.1}})
    assert cfg.get("hydraulics.root_beta", file="pft", pft=2) == 0.1
    cfg.set("hydraulics.root_beta", 0.2, file="pft", pft=2)
    assert cfg.get("hydraulics.root_beta", file="pft") == [0.1, 0.2]
    one = RunConfig(main={}, pft={"pft": {"wood_density": [0.6]}, "hydraulics": {"root_beta": 0.1}})
    one.set("hydraulics.root_beta", 0.2, file="pft", pft=1)
    assert one.get("hydraulics.root_beta", file="pft") == 0.2


def test_dumps_writes_every_digit_and_nested_tables():
    d = {"run": {"dt_slow": "1d", "slow_on": False}, "output": {"enabled": True, "fast": {"enabled": True}},
         "pft": {"vcmax25": [45.0, 0.1 + 0.2]}}
    assert tomllib.loads(dumps(d)) == d


def test_dumps_writes_the_sites_as_tables_of_an_array():
    d = {"run": {"dt_slow": "1d"}, "soil": {"depth": 2.0},
         "sites": [{"area_fraction": 0.25, "depth": 0.6}, {"area_fraction": 0.75}]}
    text = dumps(d)
    assert text.count("[[sites]]") == 2
    assert tomllib.loads(text) == d


def test_read_record_labels_the_files(tmp_path):
    rec = tmp_path / "r.csv"
    rec.write_text('source,key,index,present,value\n"/a/main.toml",run.n_threads,0,false,"1"\n'
                   '"/a/pft.toml",pft.vcmax25,2,true,"6.1500000000000000E+001"\n'
                   '"/a/main.toml",run.mode,0,false,"site"\n')
    r = read_record(rec, {"/a/main.toml": "main", "/a/pft.toml": "pft"})
    assert r[("main", "run.n_threads", 0)] == (False, 1.0)
    assert r[("pft", "pft.vcmax25", 2)] == (True, 61.5)
    assert r[("main", "run.mode", 0)] == (False, "site")
