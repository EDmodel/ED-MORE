# Building ED-MORE

How to build ED-MORE, choose a build kind, run the tests and run in parallel. The first two sections
are all a first build needs; the rest is reference: the compiler traps, the build kinds and what
each costs, and the test suite.

## What you need

- A Fortran 2018 compiler: Intel `ifx` (the everyday compiler), GNU `gfortran`, or NVIDIA
  `nvfortran` ([Compilers](#compilers)).
- CMake 3.20 or later (3.21 or later to use the presets below).
- The netCDF **C** library, with its CMake config. conda-forge ships it, and
  `environment.yml` is a conda environment that holds it and the Python tools
  ([Installing dependencies](#installing-dependencies)).

## Quick start

`ed_main` is the single entry point: it reads a config, runs the simulation, writes netCDF output
and exits. Three commands build it and check it:

```bash
mamba env create -f environment.yml    # once: netCDF and the Python tools (or conda env create)
conda activate ed
source /opt/intel/oneapi/setvars.sh    # so CMake finds ifx

FC=ifx cmake --preset dev              # 1. configure the dev kind, in build-dev/
cmake --build build-dev -j             # 2. build ed_main, libed.so and the tests
ctest --test-dir build-dev --output-on-failure   # 3. run the suite

./build-dev/ed_main examples/example03_demography/example_config_main.toml
```

Run these from the repository root. The preset takes netCDF from the active conda environment
(`$CONDA_PREFIX`); with netCDF elsewhere, add `-DCMAKE_PREFIX_PATH=<prefix>`. Without presets, the
same kind is
`cmake -S . -B build-dev -DCMAKE_Fortran_COMPILER=ifx -DED_BUILD_PYLIB=ON -DCMAKE_PREFIX_PATH=<prefix>`.

`build-dev` is where the example scripts look for `ed_main` and the Python package for
`libed.so`. The other kinds and the presets that name them are under [Build kinds](#build-kinds).

If the run cannot find `libnetcdf` or the compiler runtime, put both on the loader path:
`LD_LIBRARY_PATH=$CONDA_PREFIX/lib:$LD_LIBRARY_PATH`. Keep the `:$LD_LIBRARY_PATH`: replacing the
variable instead of extending it drops the Fortran runtime, and every test fails to load.

## netCDF is a hard dependency

There is no netCDF-free build. Every configure needs `CMAKE_PREFIX_PATH` (the presets set it to
`$CONDA_PREFIX`) to point at the netCDF-C install, and because the netCDF CMake config pulls in HDF5,
the build enables the C language too.

**ED-MORE does not use netCDF-Fortran.** It writes netCDF through the C API via
`iso_c_binding`, because netCDF-Fortran's compiled module format is gfortran-only and would not
load under ifx or nvfortran. The C ABI is compiler-agnostic, so the same `libnetcdf` links
identically under all three compilers. For the same reason ED-MORE reads TOML with its own minimal
reader (`src/config/ed_toml.f90`) rather than toml-f: a compiled Fortran module is
compiler-specific.

## Compilers

| Compiler | Role | Activation |
|---|---|---|
| Intel `ifx` | The everyday compiler. Strict standards checking; runs the full suite. | `source /opt/intel/oneapi/setvars.sh` |
| NVIDIA `nvfortran` | The second back end, and the host-multicore path. (There is no GPU build — see [Parallel builds](#parallel-builds).) | put the HPC SDK `compilers/bin` on `PATH` |
| GNU `gfortran` | Supported; ED2's reference toolchain. | usually already on `PATH` |

**A green ifx run is not sufficient.** Build new modules with nvfortran (the multicore back end) too,
and with gfortran, which on most machines is the second compiler at hand. These eight traps are
invisible to an ifx dev build:

- **Never pass an array-valued function result straight into a call.** nvfortran's whole-program
  optimizer miscompiles the temporary descriptor — silently wrong values at `-O2`, a segfault at
  `-O0` — while `ifx -stand f18 -check all` emits only a remark. Bind to a named array first. The
  same optimizer returned 0 for a flux the fast driver formed from fields it read back after the
  call that filled them, so `fast_step` returns the sensible heat flux `h_flux` itself.
- **nvfortran rejects a `BLOCK` construct** anywhere inside a parallel region.
- **ifx builds `private` and `firstprivate` copies of a derived type through a compiler-generated
  static mold** that every thread writes. Nor does ifx default-initialise a derived type declared
  in a `BLOCK` inside an outlined region: its allocatable components hold garbage, an `allocated`
  test reads them as allocated, and the first deallocate segfaults, even on one thread. This is why
  the fast loop's per-patch scratch is an explicit per-thread pool indexed by thread number rather
  than an OpenMP data-sharing clause or a block-local.
- **Never pass a component section `a(:)%c` as an actual argument when `c`'s type has
  allocatable components.** gfortran copies the section into a temporary and back, and the
  copy-out leaves the originals' allocatable components dangling: a segfault a step later, with
  nothing reported by `-fcheck=all`. Keep such objects in a contiguous array of their own, as the
  region keeps its polygons' output buffers.
- **Never write an array expression over a two-level component section `a(:)%b%c`**, directly or
  through an `associate`. ifx's AVX-512 code (`ED_CPU_TARGET=cascadelake`, which the production
  preset sets; IPO is not needed) reads such a section from the wrong offsets, while the dev build
  reads it right: a check of every PFT's leaf apoplastic fraction saw other PFTs' elastic modulus
  and osmotic potential instead, and refused every PFT file with more than one PFT. The stored
  values were right. Read such values element by element.
- **Never nest an implied-do inside a `pack` inside an array constructor**, as in
  `[pack([(i, i = 0, n - 1)], m1), pack([(i, i = 0, n - 1)], m2)]`. nvfortran 25.11 returns wrong
  elements at `-O2` and segfaults at `-O0`. Bind the implied-do to a named array and pack that.
- **Keep `findloc` off LOGICAL arrays, and never search a character array for a shorter value.**
  nvfortran 25.11's runtime aborts with "FINDLOC: unimplemented for data type" on a LOGICAL
  array; search an integer mask instead, `findloc(merge(1, 0, mask), 1)`. On a character array it
  returns 0 for a value shorter than the elements, such as `findloc(names, trim(s))`, where the
  standard compares blank-padded; pass the untrimmed value.
- **Never reset with an empty structure constructor, `x = t()`, when `t` has a fixed-size array
  component whose own type has allocatable components.** nvfortran 25.11 compiles it to an
  ALLOCATE of a garbage size, and a small program using the same pattern is an internal compiler
  error. Assign a default-initialised local that is never written, or reset through an
  `intent(out)` dummy.

## Build kinds

| Kind | Use it for | Build type and flags | Builds | Configure |
|---|---|---|---|---|
| `debug` | writing or fixing engine code: every bounds, uninitialised-value and floating-point check is on | `Debug`: `-O0 -g` and the compiler's runtime checks | `ed_main`, tests | `cmake --preset debug` |
| `dev` | everyday work: the tests, short runs, the Python package, comparing outputs between commits | `Release`: `-O2`, any x86-64 CPU | `ed_main`, `libed.so`, tests | `cmake --preset dev` |
| `profile` | finding where the time goes (VTune, perf) | dev plus `-g` | as dev | `cmake --preset profile` |
| `production` | long runs: the examples, calibration, regional runs | dev plus `ED_IPO=ON` and `ED_CPU_TARGET=cascadelake` | `ed_main`, `libed.so` | `cmake --preset production` |

Each preset builds in `build-<kind>` and takes netCDF from the active conda environment
(`$CONDA_PREFIX`). Run it from the repository root, and choose the compiler with `FC=ifx` (or
`gfortran`) before `cmake --preset`; presets need CMake 3.21 or later. Then build and test as in the
[quick start](#quick-start), with `build-<kind>` for `build-dev`.

A `-D` after `--preset` overrides the preset, and `-B` its directory:
`cmake --preset dev -DCMAKE_PREFIX_PATH=<netCDF prefix>`. Without presets, the same kinds are
`-DCMAKE_BUILD_TYPE=Debug|Release`, plus `-DCMAKE_Fortran_FLAGS=-g` for profile and
`-DED_IPO=ON -DED_CPU_TARGET=cascadelake` for production, and `-DED_BUILD_PYLIB=ON` (which adds
`libed.so`) for every kind but debug. **A configure that names no build type builds Release (dev)**,
without `libed.so` unless it adds `-DED_BUILD_PYLIB=ON`. `-DED_BUILD_TESTS=OFF` skips the tests
(the production preset does). The per-compiler flags live in the `ed_fortran_flags()` function in
`CMakeLists.txt`.

**`build-dev` is the default build.** dev builds `ed_main`, `libed.so` and the tests, and the
example scripts run `build-dev/ed_main` and the Python package loads `build-dev/libed.so`, unless
`--ed-main` (or `ED_MAIN`) or `ED_LIB` names another build. The library costs dev no
measurable build time, because the tests compile its sources anyway and it adds one link: with or
without it, a clean configure and build at `-j 16` takes 44–50 s with ifx and 52–56 s with gfortran,
and the rebuild after an edit 3–4 s. It compiles every library position-independent (`-fPIC`),
which leaves every output value of `ed_main` as it is, on ifx and gfortran, and its run time
within noise. The suite includes the Python package's tests: 86 tests.

**dev is the reference.** It is the project's `-O2`, for any x86-64 CPU. Outputs are compared between
commits on dev builds; gfortran dev builds of two commits that change no arithmetic give
byte-identical output. A profile build gives dev's output byte for byte. A production build gives
the same results to round-off ([Production's output against dev's](#productions-output-against-devs)),
so do not compare a production run's output byte for byte with a dev run's.

**production is 26–29 % faster than dev** (ifx, example 04's 120-patch stand, 1 and 20 threads). Two
things give it: inlining across modules (`ED_IPO`), and code for the CPU family
(`ED_CPU_TARGET=cascadelake`, which gains more with more threads; it covers cbsuxu 01–10). Under
IPO the files compile at `-O2`, and ifx generates the code at the link at `-O3`, which the build
names; gfortran keeps each file's `-O2` through the link. A production binary may not start on an
older CPU: build with `-DED_CPU_TARGET=` (empty) for a portable one, or `host` for the machine you
build on. Use ifx: gfortran's production build is only 4–8 % faster than ifx dev. IPO and the CPU
target are offered for ifx and gfortran; nvfortran builds debug and dev.

**production builds `ed_main` and `libed.so`, not the tests.** Under IPO each of the 47 test
programs is a whole-program link of its own, and several tests assert results to the last bit, which
a production build moves at round-off. With `-DED_BUILD_TESTS=ON`, 4 tests fail with ifx and 2 with
gfortran ([The suite runs on dev](#the-suite-runs-on-dev)). Run the suite on dev. For the same
reason `ed.model` and `ed_main` agree to round-off on a production build, not bit for bit: each
is optimised as a whole program of its own.
To run the examples or a calibration on a production build, point them at it: a script that runs
`ed_main` takes `--ed-main` (or `ED_MAIN`) and otherwise runs `build-dev/ed_main`, and the
Python API takes `ED_LIB=build-production/libed.so` and otherwise loads `build-dev/libed.so`.

**Why development uses dev.** With ifx, a production build relinks the whole program at every link:
an edit to one file takes ~35 s to rebuild instead of ~4 s, and rebuilding the tests ~85 s instead of
~11 s. It needs under 3 GB of memory, tests included. For a long run without that cost, add
`-DED_CPU_TARGET=cascadelake` to a dev build in a directory of its own (`-B build-dev-cascadelake`):
most of the gain at 20 threads, none of the relink.

**profile** is dev with `-g`, which VTune and perf read; gprof also needs `-pg`
(`cmake --preset profile -DCMAKE_Fortran_FLAGS="-g -pg"`). **debug** is many times slower and prints
one harmless `arg_temp_created` remark per netCDF call; prefer dev when the output layer is what you
are exercising.

## Parallel builds

**Host threading** is compiled in by default; a run asks for it. A site run threads its patches, a
region run its polygons (`docs/configuration.md`, "Regional runs"):

```bash
FC=ifx cmake --preset dev
# then set [run].n_threads in the TOML (default 1), and give the job that many cores.
```

Output is **byte-identical at any thread count**, and the test suite asserts it.

- `-DED_OPENMP=OFF` builds serial. A compiler with no Fortran OpenMP falls back to serial with a
  CMake warning.
- A build directory whose CMake cache holds `ED_OPENMP=OFF` keeps it. Pass `-DED_OPENMP=ON` or
  configure a fresh directory.
- The Python wheel builds serial (`python/pyproject.toml`), so `libed.so` never brings a second
  OpenMP runtime into a Python process.
- gfortran with a conda netCDF prefix warns that the prefix's `libgomp.so.1` may hide the
  compiler's. That is harmless when the conda copy is the newer one, as it is in the `ed`
  environment (GCC 16 against the system's GCC 11), and the test suite passes that way.

**An OpenMP build does not reproduce a serial build bit for bit.** Compiling the patch loop as a
parallel region changes its rounding. On the BCI census example the two builds agree for 62 hours,
then differ in the 13th significant digit, and the drift reaches about 2% of the largest hourly flux
after five years. The tower statistics and the stand agree to every printed digit. The stack flag
described next is not the cause: a serial build with only `-auto` reproduces the serial build exactly.

`ED_OPENMP` does two things, and the second is load-bearing: it puts the OpenMP flag on the
fast-loop target, **and it adds the per-compiler "all locals on the stack" flag** (`-auto`,
`-frecursive`, `-Mrecursive`) to *every* target. Intel Fortran defaults to `-auto-scalar`, which
places local arrays and derived types in static storage shared by every thread; without that flag
the kernels race and return plausible, silently thread-count-dependent numbers.

**What threads buy.** The BCI flux-tower example (130 census patches, falling to about 57 over the
five years), ifx Release, each run alone on an idle `R128C40` node; the 60-day times are the second
of two runs:

| `n_threads` | 60 days | five years |
|---|---|---|
| 1 | 82 s | 26 min 0 s |
| 4 | 30 s | 7 min 53 s |
| 8 | 20 s | 4 min 34 s |
| 16 | 16 s | 3 min 12 s |
| 32 | 14 s | 2 min 37 s |

- **The patches are the parallel axis,** so a site run gains while it has more patches than
  threads: at BCI 16 threads run 8× faster than one, and 32 only 20 % faster than 16.
- **No routine passes a procedure as an argument.** ifx allocates a lock-guarded record on every
  call of a routine that hands one of its contained functions to another routine, so the threads
  queue on it: two such routines in the fast loop are enough to make 8 and 16 threads slower than 4
  (41.5 s and 56.7 s on the 60-day case). The `no_procedure_arguments` test keeps every ED-MORE source free of
  them.

**Large runs need a big stack.** The fast loop keeps automatic arrays on the stack, and a threaded
run takes each thread's stack from `OMP_STACKSIZE`, whose small default overflows on a large region:
`export OMP_STACKSIZE=512M` before it (and `ulimit -s unlimited` for a serial run with very many
cohorts). `-DED_HEAP_ARRAYS=ON` moves large automatic arrays to the heap instead, but trips an ifx
fault in a debug build, so prefer the runtime setting.

**There is no GPU build.** OpenMP `target` offload under nvfortran measured 1.4× slower than the
CPU (49.8 s against 36.2 s), with one kernel at 0.4 % occupancy. The measurement and what to do
instead are in the archived GPU evaluation:
<https://github.com/xiangtaoxu/MEDS/blob/v0.4.0/docs/dev_plans/archive/MEDS_GPU_EVALUATION.md>.
nvfortran threads the host through `ED_OPENMP`, as the other compilers do.

## Installing dependencies

Helper scripts check whether a dependency is already present and prompt before changing anything
(`-y` to skip the prompt, `-h` for help). They target Debian, Ubuntu and WSL.

```bash
./scripts/install_netcdf.sh            # netCDF C library, via conda-forge or apt (--conda or --apt forces one)
./scripts/install_ifx.sh               # Intel ifx via the oneAPI APT repository
./scripts/install_gfortran.sh          # GNU gfortran
./scripts/install_gfortran.sh --hdf5   # ... plus HDF5 dev files
```

`install_netcdf.sh` detects an existing netCDF via `nc-config` and prints the
`-DCMAKE_PREFIX_PATH` to pass. conda-forge is recommended because it ships the CMake config.

## Python and post-processing

The post-processing scripts need Python with numpy, matplotlib, netCDF4 and pillow. The tools'
tests also need pandas and pytest. The conda environment in `environment.yml`
provides all of them:

```bash
mamba env create -f environment.yml    # or: conda env create -f environment.yml
conda activate ed
```

To drive the model from Python, `pip install python/` compiles the Fortran and bundles
`libed.so` beside the package; from the source tree, `PYTHONPATH=python` and a dev build are
enough. See `python/README.md`.

## Running the tests

```bash
ctest --test-dir build-dev --output-on-failure                    # all of them
ctest --test-dir build-dev -L unit --output-on-failure            # one label (below)
ctest --test-dir build-dev -L integration --output-on-failure
ctest --test-dir build-dev -LE long --output-on-failure           # all but the long runs
ctest --test-dir build-debug -R slow_update --output-on-failure   # one, by regex
```

The suite covers allometry round-trips, carbon and area conservation, container integrity, the
individual physical kernels, the fast step on both integrators, the output round-trip, the
conservation ledgers, whole runs and restarts. Every C-API shim is compiled by a default test
target, so an ABI change fails the build instead of breaking the optional Python library unseen.

**One test per `src` module, laid out as `src` is.** ctest `slow_update` is
`test/slow/test_slow_update.f90`, the test of `src/slow_dynamics/driver/ed_slow_update.f90`; the
folders are `test/share/`, `config/`, `state/`, `fast/`, `slow/`, `io/`, `driver/` and `c_api/`, with
the shared fixtures in `test/support/`, the CMake checks in `test/checks/` and the restart test in
`test/python/`. `test/CMakeLists.txt` registers every test through `ed_add_test` and builds the
test programs in `<build>/test/`. Each test carries one label:

| Label | What it holds | Tests (dev) | Time (dev, one core) |
|---|---|---|---|
| `unit` | each module's test, and `ed_main`'s refusals (`*_refused`, `main_failure_finalizes`) | 67 | 15 s |
| `integration` | whole runs: `sites`, `region`, `fast_driver`, `c_api_run`, `restart_exact` | 5 | 44 s |
| `lint` | the source-tree rules (`no_issue_refs`, `constants_named`, `kernels_patch_only`, ...) and the tracked copies (`io_config_example`, `license_copy`, `notice_copy`) | 10 | 2 s |
| `python` | the Python tools' and the package's pytest suites (`prepare_flux_tower`, `calibrate_fast`, `python_api`, `prepare_census`) | 4 | 26 s |
| `long` | opt-in: example 05's Palo Verde inventory and its 150-year regrowth | 2 | 34 min on 12 threads |

A debug build has 84 tests: no `python_api` (it builds no `libed.so`), and no
`soil_nan_layer_refused` (it traps the NaN that case feeds the loader). A `-L` takes a regular
expression, so `-L 'unit|lint'` runs both.

**The long tests are opt-in**, and the runs to make before merging a change to the physics or to the
tissue energy. Configure with them, and run them alone:

```bash
FC=ifx cmake --preset dev -DED_TEST_LONG=ON -DED_ERA5LAND_ARCHIVE=/path/to/ED_ERA5land
ctest --test-dir build-dev -L long --output-on-failure
```

Both runs are in `examples/example05_forest_regeneration/palo_verde/`:

- `palo_verde_inventory` runs the 2009 inventory (`ed_config_inventory.toml`; 85 s on the
  example's 12 threads, ifx dev), reading the archive.
- `palo_verde_regrowth` runs 150 years from bare ground (`ed_config_regrowth.toml`; 31 minutes) on
  Palo Verde's cell, 2003–2022. It cuts that cell from the archive into its directory once (under
  2 minutes) unless `-DED_TEST_PV_FORCING=` names a cut file.

Each runs a copy of the example's config in `<build>/long/palo_verde_<case>/`, which holds its
output and its log, and passes when the run completes and its slow ledger closes
(`test/checks/check_long_run.cmake`). `ED_ERA5LAND_ARCHIVE` defaults to `ERA5LAND_ARCHIVE` from the
environment.

### The suite runs on dev

**Run the suite on a dev build** (or debug); dev passes all 86 tests with ifx and with gfortran.
Several tests compare two computations of the same number bit for bit, or within a few units in the
last place, and that holds only when both are compiled the same way. dev compiles each file on its
own at `-O2`, for any x86-64 CPU. production optimises across files (IPO) and generates AVX-512 code,
so the two sides of such a check can be compiled differently and part in the last bits. Built with
`-DED_BUILD_TESTS=ON`, production fails 4 tests with ifx and 2 with gfortran:

| Test | Compiler | The check that fails | Production's value | Tolerance | Cause |
|---|---|---|---|---|---|
| `fast_driver` | ifx | the fast loop stages each sub-step's forcing sample exactly | the test's own recomputation of a sample does not equal the staged one; recomputed in a routine of its own, all eight fields match exactly, and printing the difference in place makes the check pass | exact | IPO |
| `plant_hydraulics` | ifx | a leaf above its closing water content transpires its demand | 1.9999999999999998e-5 against 2e-5: one unit in the last place (3e-21, 2e-16 relative); dev gives 2e-5 exactly | exact | IPO |
| `plant_hydraulics` | ifx | transpiration's derivative matches finite differences | where the derivative is 0, the finite difference is 2.6e-14, from round-off in the two transpirations it differences 2e-7 apart; dev's is 0 | 1e-15 | IPO |
| `plant_hydraulics` | gfortran | transpiration's derivative matches finite differences | the same, 5.8e-15 | 1e-15 | FMA under LTO |
| `python_api` | ifx | `ed.model` writes `ed_main`'s output bit for bit | 101 of 301 variables differ, by up to 1.3e-12 of a variable's largest value (patch sensible heat), the near-zero conservation residuals aside | exact | IPO |
| `python_api` | gfortran | `ed.model` writes `ed_main`'s output bit for bit | 104 of 301 variables differ, by up to 9.5e-14 (patch sensible heat), the residuals aside | exact | FMA under LTO |
| `sites` | ifx | each of two identical sites holds the one-site stand's plant carbon | 9.696054 against 9.696198 kgC m⁻² after 60 days, 1.5e-5 relative; the two sites agree with each other to 2e-16 | 1e-12 relative | AVX-512 vectorisation |

The causes come from building each part of production on its own: IPO alone, the CPU target alone,
and each with fused multiply-add (FMA) turned off (`-no-fma`, `-ffp-contract=off`).

- **IPO**: the test fails with IPO alone and passes with the CPU target alone. Each program is
  optimised as a whole, so one routine compiles differently where it is inlined: in the test and in
  the driver, or in `ed_main` and in `libed.so`. With ifx, turning FMA off leaves the same four
  failures.
- **FMA under LTO**: gfortran fails only with LTO and the CPU target together, and passes with FMA
  off. Inlining across files lets FMA contract expressions that each file alone leaves apart.
- **AVX-512 vectorisation**: the test fails with the CPU target alone, with or without FMA, and
  passes with IPO alone. The rounding of a vectorised loop depends on how many elements it covers,
  so a site in a two-site polygon parts from the same site run alone in the last bit, and the run's
  fusion and disturbance carry that to 1.5e-5 in 60 days. It is the size of the difference between
  dev builds with the two compilers (1.6e-5 in the same quantity), and production's one-site value
  lies 3.6e-5 from dev's.

### Production's output against dev's

Ten runs, each with a dev and a production build, one thread: the largest difference in any hourly
value (daily where the run writes no hourly output), over ifx and gfortran, and the largest change
of a run mean.

| Run | LE (W m⁻²) | H (W m⁻²) | Canopy air (K) | GPP (µmol m⁻² s⁻¹) | Largest run-mean change |
|---|---|---|---|---|---|
| A census on the reference climate, 2 years (daily) | 3.9e-5 | 4.4e-5 | 1.8e-6 | 3.3e-5 | 3.5e-6 (NEE) |
| The same, 1 year, hourly | 6.9e-4 | 5.3e-4 | 1.2e-6 | 4.6e-4 | 8.4e-7 (NEE) |
| A restart from the 2-year run, 6 months (daily) | 5.9e-7 | 4.4e-5 | 1.8e-6 | 5.8e-7 | 1.3e-5 (GPP) |
| Example 04, the BCI census on tower forcing, 1 month | 7.7e-5 | 7.4e-5 | 1.2e-6 | 1.2e-7 | 2.2e-9 (H) |
| The same on RK45 | 0.10 | 0.093 | 3.5e-3 | 1.8e-4 | 6.6e-6 (H) |
| Example 05, BCI from bare ground on ERA5-Land, 1 year | 0.085 | 0.031 | 1.5e-3 | 6.4e-4 | 1.8e-5 (NEE) |
| The same polygon as two sites | 0.059 | 0.025 | 1.1e-3 | 5.0e-4 | 1.4e-5 (NEE) |
| A region of ERA5-Land cells at Palo Verde, 2 months (daily) | 0.23 | 0.089 | 3.0e-3 | monthly only | 4.0e-4 (NEE) |
| A restart of the region, 1 month (daily) | 0.12 | 0.063 | 3.3e-3 | monthly only | 1.5e-3 (NEE) |
| Ithaca from bare ground under snow, 2 winter months | 1.6e-6 | 8.2e-7 | 3.4e-8 | 0 in both | 2.6e-10 (LE) |

Every run mean of LE, H and GPP stays within 2.1e-5 of dev's, relative; NEE, the small difference
of two large fluxes, moves most. Most values differ far less than the largest: the median
difference in hourly LE is below 3e-5 W m⁻² in every run that writes hourly output. A dev build with
ifx and one with gfortran differ by the same amounts. **Production's results are dev's to round-off,
carried through a run by its thresholds and adaptive steps, so compare a production run with dev's
within tolerances, never byte for byte.**
