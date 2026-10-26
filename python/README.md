# `meds` — Python interface to MEDS

`meds` wraps the compiled Fortran model in plain Python: call the leaf and phenology kernels, drive
the demography engine, run the whole coupled model with the time loop in Python, or edit and write
the config files of a run. One shared library, `libmeds.so`, backs every sub-module, so there is one
binary and one ABI.

| Sub-module | What it gives you | Needs `libmeds.so` |
|---|---|---|
| [`meds.plant.leaf`](meds/plant/leaf.py) | leaf gas exchange: FvCB C3 and Collatz C4 photosynthesis, Leuning, Medlyn and Katul stomata, the coupled A–gs–Ci solver | yes |
| [`meds.plant.pheno`](meds/plant/pheno.py) | leaf phenology: the flush and shed tendencies, the coupled model's leaf rule, four leaf-habit presets (see `examples/example02_canopy_phenology/`) | yes |
| [`meds.canopy`](meds/canopy/__init__.py) | the canopy's leaf solve and two-stream radiation over a frozen stand, parameterized by a run's own config | yes |
| [`meds.demography`](meds/demography/_site.py) | `Config` and `Site`: load a TOML config, build a site, and step the slow dynamics or feed them externally computed rates (`apply_rates`); the engine `meds_main` runs, behind an opaque handle | yes |
| [`meds.model`](meds/model/_run.py) | `Run` and `run`: the whole coupled model as `meds_main` runs it, with the time loop in Python | yes |
| [`meds.config`](meds/config.py) | `RunConfig`: a run's main TOML and the PFT (plant functional type) TOML it names | no |

`import meds` loads no library; each sub-module loads `libmeds.so` on first import.

## Examples

Leaf gas exchange and phenology:

```python
import meds.plant.leaf as leaf
params = leaf.c3_params(vcmax25=60.0, jmax25=108.0)          # or leaf.c4_params(...)
flux = leaf.gas_exchange(par=1500.0, leaf_temp=298.15,       # leaf_temp in KELVIN
                         vpd=1000.0, ca=400.0, params=params,
                         stomata=leaf.Stomata.MEDLYN)
print(flux.A_net, flux.gs, flux.ci, flux.limitation, flux.converged)

import meds.plant.pheno as pheno
ph = pheno.Phenology(pheno.temperate_deciduous())            # params + cue memory + leaf cover
day = ph.step(temp_day=290.0, par_hours=14.5, doy=150)       # par_hours: hours of PAR above par_min
print(day.leaf_flush_tendency, day.leaf_cover)               # tendency [0-1], leaf / full leaf
```

The whole model: `run(config)` runs a config to its end, and `python -m meds.model CONFIG` does the
same from the shell, in place of `meds_main CONFIG`. A `Run` steps one slow step at a time:

```python
from meds.model import Run

with Run("meds_config_main.toml") as run:
    for step in run:
        if step.is_new_year:
            print(step.date, run.total_agb, run.soil_carbon)
```

On a dev build the output is `meds_main`'s, bit for bit; on a production build it agrees to round-off
([`docs/building.md`](../docs/building.md), "Build kinds").

Editing a run's config: `RunConfig` reads the main TOML and the PFT TOML it names, changes keys (one
PFT's element of a trait array included), writes the pair for a run, and `meds.config.read_record` reads
back a run's parameter record to check what the model read. It needs no compiled library, and `scripts/calibrate_fast`
builds every trial with it.

```python
from meds.config import RunConfig
from meds.model import run
cfg = RunConfig.load("examples/example04_column_biophysics/meds_config_eval.toml")
cfg.set("pft.vcmax25", 50.0, file="pft", pft=1)              # the first PFT's top-of-canopy Vcmax25
cfg.set("run.end_time", "2012-09-01 00:00:00")
run(cfg.write("trial"))                                       # trial/main.toml + trial/pft.toml
```

## Install

`pip install python/` compiles the Fortran and bundles the library: scikit-build-core drives the
top-level CMake, and the wheel ships `meds/libmeds.so` beside the Python modules. Installing needs a
Fortran compiler, CMake and the netCDF C library; running needs `numpy` and `netCDF4` (and `tomli`
below Python 3.11).

```bash
source /opt/intel/oneapi/setvars.sh                    # so CMake finds ifx
CMAKE_PREFIX_PATH=$CONDA_PREFIX pip install python/    # netCDF prefix; builds and bundles libmeds.so
python -m meds.plant                                    # round-trip self-test
```

The wheel is built serial (`MEDS_OPENMP=OFF`), so `libmeds.so` never brings a second OpenMP runtime
into a Python process.

**The installed library needs no `LD_LIBRARY_PATH`.** Its RPATH is written at build time from the
directories it linked (netCDF, HDF5 and the Fortran runtime), so it loads in a bare venv. The
corollary is that **the wheel is machine-local**: those paths are this machine's conda prefix and
oneAPI install. It is not a redistributable manylinux wheel, which would need a gfortran build,
`auditwheel` and the netCDF dependency closure (about 67 MB) vendored.

### Dev loop: edit Python without recompiling

```bash
# 1. Build the library once from the repository root, with the dev preset:
FC=ifx cmake --preset dev
cmake --build build-dev --target meds_py                # -> build-dev/libmeds.so

# 2. Run against the source tree (the Fortran runtime must be on LD_LIBRARY_PATH here):
source /opt/intel/oneapi/setvars.sh
PYTHONPATH=python python -m meds.plant
PYTHONPATH=python pytest python/tests
```

`meds._libmeds` looks for the library in this order:

1. the path in the `MEDS_LIB` environment variable;
2. `meds/libmeds.so` beside the package (how a wheel ships it);
3. `build-dev/libmeds.so` in the source tree, the dev preset's build.

Both paths above need no configuration. To load another build, name it:
`MEDS_LIB=build-production/libmeds.so`.

## License

Apache License 2.0, like the rest of MEDS. See [`LICENSE`](LICENSE) and [`NOTICE`](NOTICE).
