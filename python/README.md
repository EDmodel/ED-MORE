# `edmore` — Python interface to ED-MORE

`edmore` wraps the compiled Fortran model in plain Python: call the leaf and phenology kernels, drive
the demography engine, run the whole coupled model with the time loop in Python, or edit and write
the config files of a run. One shared library, `libedmore.so`, backs every sub-module, so there is one
binary and one ABI.

| Sub-module | What it gives you | Needs `libedmore.so` |
|---|---|---|
| [`edmore.plant.leaf`](edmore/plant/leaf.py) | leaf gas exchange: FvCB C3 and Collatz C4 photosynthesis, Leuning, Medlyn and Katul stomata, the coupled A–gs–Ci solver | yes |
| [`edmore.plant.pheno`](edmore/plant/pheno.py) | leaf phenology: the flush and shed tendencies, the coupled model's leaf rule, four leaf-habit presets (see `examples/example02_canopy_phenology/`) | yes |
| [`edmore.canopy`](edmore/canopy/__init__.py) | the canopy's leaf solve and two-stream radiation over a frozen stand, parameterized by a run's own config | yes |
| [`edmore.demography`](edmore/demography/_site.py) | `Config` and `Site`: load a TOML config, build a site, and step the slow dynamics or feed them externally computed rates (`apply_rates`); the engine `edmore` runs, behind an opaque handle | yes |
| [`edmore.model`](edmore/model/_run.py) | `Run` and `run`: the whole coupled model as `edmore` runs it, with the time loop in Python | yes |
| [`edmore.config`](edmore/config.py) | `RunConfig`: a run's main TOML and the PFT (plant functional type) TOML it names | no |
| [`edmore.forcing`](edmore/forcing.py) | `write_forcing_file`: the one writer of a forcing file (`[forcing].format = "default"`), which every tool that makes one uses | no |

`import edmore` loads no library; each sub-module loads `libedmore.so` on first import.

## Examples

Leaf gas exchange and phenology:

```python
import edmore.plant.leaf as leaf
params = leaf.c3_params(vcmax25=60.0, jmax25=108.0)          # or leaf.c4_params(...)
flux = leaf.gas_exchange(par=1500.0, leaf_temp=298.15,       # leaf_temp in KELVIN
                         vpd=1000.0, ca=400.0, params=params,
                         stomata=leaf.Stomata.MEDLYN)
print(flux.A_net, flux.gs, flux.ci, flux.limitation, flux.converged)

import edmore.plant.pheno as pheno
ph = pheno.Phenology(pheno.temperate_deciduous())            # params + cue memory + leaf cover
day = ph.step(temp_day=290.0, par_hours=14.5, doy=150)       # par_hours: hours of PAR above par_min
print(day.leaf_flush_tendency, day.leaf_cover)               # tendency [0-1], leaf / full leaf
```

The whole model: `run(config)` runs a config to its end, and `python -m edmore.model CONFIG` does the
same from the shell, in place of `edmore CONFIG`. A `Run` steps one slow step at a time:

```python
from edmore.model import Run

with Run("ed_config_main.toml") as run:
    for step in run:
        if step.is_new_year:
            print(step.date, run.total_agb, run.soil_carbon)
```

On a dev build the output is `edmore`'s, bit for bit; on a production build it agrees to round-off
([`docs/building.md`](../docs/building.md), "Build kinds").

Editing a run's config: `RunConfig` reads the main TOML and the PFT TOML it names, changes keys (one
PFT's element of a trait array included), writes the pair for a run, and `edmore.config.read_record` reads
back a run's parameter record to check what the model read. It needs no compiled library, and `scripts/calibrate_fast`
builds every trial with it.

```python
from edmore.config import RunConfig
from edmore.model import run
cfg = RunConfig.load("examples/example04_column_biophysics/ed_config_eval.toml")
cfg.set("pft.vcmax25", 50.0, file="pft", pft=1)              # the first PFT's top-of-canopy Vcmax25
cfg.set("run.end_time", "2012-09-01 00:00:00")
run(cfg.write("trial"))                                       # trial/main.toml + trial/pft.toml
```

## Install

`pip install python/` compiles the Fortran and bundles the library: scikit-build-core drives the
top-level CMake, and the wheel ships `edmore/libedmore.so` beside the Python modules. Installing needs a
Fortran compiler, CMake and the netCDF C library; running needs `numpy` and `netCDF4` (and `tomli`
below Python 3.11).

```bash
source /opt/intel/oneapi/setvars.sh                    # so CMake finds ifx
CMAKE_PREFIX_PATH=$CONDA_PREFIX pip install python/    # netCDF prefix; builds and bundles libedmore.so
python -m edmore.plant                                 # round-trip self-test
```

The wheel is built serial (`ED_OPENMP=OFF`), so `libedmore.so` never brings a second OpenMP runtime
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
cmake --build build-dev --target ed_py                # -> build-dev/libedmore.so

# 2. Run against the source tree (the Fortran runtime must be on LD_LIBRARY_PATH here):
source /opt/intel/oneapi/setvars.sh
PYTHONPATH=python python -m edmore.plant
PYTHONPATH=python pytest python/tests
```

`edmore._libedmore` looks for the library in this order:

1. the path in the `EDMORE_LIB` environment variable;
2. `edmore/libedmore.so` beside the package (how a wheel ships it);
3. `build-dev/libedmore.so` in the source tree, the dev preset's build.

Both paths above need no configuration. To load another build, name it:
`EDMORE_LIB=build-production/libedmore.so`.

## License

Apache License 2.0, like the rest of ED-MORE. See [`LICENSE`](LICENSE) and [`NOTICE`](NOTICE).
