# ED-MORE: a Modular, Open, Reproducible, and Extensible Ecosystem Demography model

ED-MORE is a size- and age-structured terrestrial biosphere model with fully coupled carbon, energy
and water dynamics, written in Fortran 2018. A forest is cohorts of plants on patches of land, with
fusion and fission keeping both bounded. At every sub-daily step, ED-MORE solves the carbon, energy
and water dynamics of the canopy air, the leaves, the wood, the snow and the soil together, so leaf
temperature, canopy humidity and soil moisture are computed rather than read from a driver file. The
demography then grows, kills and recruits plants on the carbon this produces. A run can be a single
site, a polygon divided into sites with their own soil, elevation and slope, or a region in which
every ERA5-Land grid cell within a box is its own polygon.

**Status: v0.4.1 in preparation, and not yet benchmarked.** The examples fit some parameters to
observations (a flux tower, a forest census, leaf area and litter) and set a regional run against
satellite leaf area, but no regression suite against reference runs exists yet, and no systematic
evaluation against observations has been made. What is verified is internal: the test suite on two
compilers, the per-step conservation budgets, and thread invariance. Treat the numbers a run produces
as a working model's numbers. The [changelog](CHANGELOG-v0.md) records which numbers each release
moved.

## Quick start

```bash
./scripts/install_netcdf.sh                  # if you do not have netCDF; prints the prefix

FC=ifx cmake --preset dev                    # the dev build, in build-dev/; netCDF from $CONDA_PREFIX
cmake --build build-dev -j
ctest --test-dir build-dev --output-on-failure

# A forest regrowing from near-bare ground for 300 years (example 03; about a minute)
PYTHONPATH=python python examples/example03_demography/run_demography.py --start bare
python post_proc/plot_site_timeseries.py examples/example03_demography/output/bare.nc -o timeseries.png
```

netCDF is required; gfortran, ifx and nvfortran are supported. [`docs/building.md`](docs/building.md) covers the four
build kinds (debug, dev, profile, and production for long runs), the parallel builds and the compiler
traps. The dev build also makes the shared library the Python front end loads ([`python/README.md`](python/README.md));
the plots need numpy, matplotlib and netCDF4. Each example's README says which data it needs and how
to run it.

## Where to go next

| If you want to | Read |
|---|---|
| learn what is in the model: each process and its equations | [`docs/science/`](docs/science/), one page per process, indexed in [`docs/README.md`](docs/README.md) |
| run ED-MORE, coming from ED2 | [`docs/ed2_comparison.md`](docs/ed2_comparison.md): the main differences first, then process by process, keyed to `ED2IN` |
| set up a run | [`docs/configuration.md`](docs/configuration.md): the two TOML files and the three ways to start |
| see worked examples | `examples/`: leaf gas exchange, phenology, demography, a flux tower, forest regrowth, a region |
| prepare forcing and census inputs | `scripts/prepare_era5/`, `scripts/prepare_flux_tower/`, `scripts/prepare_census/` |
| drive it from Python | [`python/README.md`](python/README.md): the `edmore` package |
| plot the output | `post_proc/README.md` |
| change the code | [Developing ED-MORE](#developing-ed-more) below, then [`CONTRIBUTING.md`](CONTRIBUTING.md) |
| know what changed | [`CHANGELOG-v0.md`](CHANGELOG-v0.md) |

## Developing ED-MORE

The name states four principles. The rules that follow from them, each with the check that enforces
it where one exists, are in [`AGENTS.md`](AGENTS.md): written so that AI coding agents can follow
them, and binding on everyone who changes the code. How to propose, send and merge a change is in
[`CONTRIBUTING.md`](CONTRIBUTING.md).

- **Modular: one job per piece, and one place for each piece.** The source splits by timescale, a
  fast sub-daily loop and a slow daily one, over a state layer both share. A process kernel sees only
  the stores it works on, dependencies point one way, and there is no hidden global state: config and
  state are passed as arguments, so every routine can be tested alone. [`src/README.md`](src/README.md)
  says where a new file goes.
- **Open: anyone can read the model, check it and reuse it.** The code is under Apache 2.0. Every
  equation is published in [`docs/science/`](docs/science/), versioned with the code that implements it, and the code
  cites it. Comments say what the code does and why, in plain words. Every parameter value carries
  its source, and decisions are made in the open, on issues and pull requests.
- **Reproducible: the same inputs give the same answer, and every change to the answer is
  declared.** Every parameter comes from the config files, and a run refuses what it cannot do right
  rather than fall back in silence. Output is identical at any thread count, a state file restarts a
  run bit for bit, and every budget is checked every step. Every pull request states its effect on
  results.
- **Extensible: new science goes in without rewiring what is there.** A process is a module behind
  the interface its driver already calls, a choice of scheme is a selector, and adding a config key, a
  per-cohort field or an output variable follows a fixed recipe. Outside code drives the model through
  its C API, with the same operators the model uses.

The code is modern Fortran 2018: modules, derived types, explicit interfaces, `allocatable`
ownership, parameterized kinds, `pure` and `elemental` kernels, `implicit none` everywhere, and at
most 132 columns. [`docs/building.md`](docs/building.md) lists the compiler traps: code that one compiler builds
correctly and another miscompiles.

## Origins

ED-MORE carries forward, from version 0.4.1, the development of its prototype, MEDS. The MEDS
repository, with its history, issues and pull requests up to v0.4.0, is archived at
<https://github.com/xiangtaoxu/MEDS>.

## Citing ED-MORE

Cite the version you used, as [`CITATION.cff`](CITATION.cff) gives it; GitHub's "Cite this
repository" button formats it. The science ED-MORE builds on:

- Moorcroft, P. R., Hurtt, G. C. and Pacala, S. W. (2001). A method for scaling vegetation dynamics:
  the Ecosystem Demography model (ED). *Ecological Monographs* 71(4), 557–586.
  <https://doi.org/10.1890/0012-9615(2001)071[0557:AMFSVD]2.0.CO;2>
- Medvigy, D., Wofsy, S. C., Munger, J. W., Hollinger, D. Y. and Moorcroft, P. R. (2009).
  Mechanistic scaling of ecosystem function and dynamics in space and time: Ecosystem Demography
  model version 2. *Journal of Geophysical Research: Biogeosciences* 114.
  <https://doi.org/10.1029/2008JG000812>
- Xu, X., Medvigy, D., Powers, J. S., Becknell, J. M. and Guan, K. (2016). Diversity in plant
  hydraulic traits explains seasonal and inter-annual variations of vegetation dynamics in seasonally
  dry tropical forests. *New Phytologist* 212(1), 80–95. <https://doi.org/10.1111/nph.14009>
- Longo, M., et al. (2019). The biophysics, ecology, and biogeochemistry of functionally diverse,
  vertically and horizontally heterogeneous ecosystems: the Ecosystem Demography model, version 2.2 –
  Part 1: Model description. *Geoscientific Model Development* 12, 4309–4346.
  <https://doi.org/10.5194/gmd-12-4309-2019>

These papers describe ED and its successor ED2, which is developed at
<https://github.com/EDmodel/ED2>.

## Licence

ED-MORE is released under the [Apache License 2.0](LICENSE). [`NOTICE`](NOTICE) carries the
copyright notice and the attributions to MEDS and ED2, and [`AUTHORS`](AUTHORS) lists the copyright
holders. Contributions are welcome under the same licence; see [`CONTRIBUTING.md`](CONTRIBUTING.md).
