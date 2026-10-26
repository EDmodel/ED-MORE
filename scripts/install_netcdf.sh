#!/usr/bin/env bash
# SPDX-License-Identifier: Apache-2.0
#==========================================================================================#
# install_netcdf.sh — check for the netCDF C library and offer to install it.              #
#                                                                                          #
# ED-MORE writes its output through the netCDF *C* library (via iso_c_binding). EVERY build  #
# needs it (netCDF is a hard dependency — always compiled); CMake locates it with            #
# `find_package(netCDF CONFIG)`, which you point at the install via -DCMAKE_PREFIX_PATH=<prefix>. #
#                                                                                          #
# Two routes, both Linux (e.g. WSL):                                                        #
#   * conda / conda-forge (recommended): `libnetcdf` ships the CMake config that            #
#     find_package(netCDF CONFIG) consumes directly. This matches environment.yml.          #
#   * Debian/Ubuntu apt: `libnetcdf-dev`. Note Debian's package may NOT ship the upstream    #
#     CMake config, so `find_package(netCDF CONFIG)` can fail to locate it — the conda route #
#     is the smoother path for this project.                                                #
#                                                                                          #
# Usage:                                                                                    #
#   ./scripts/install_netcdf.sh          # detect; prompt before installing (conda if a     #
#                                        #   conda env is active, otherwise apt)            #
#   ./scripts/install_netcdf.sh --conda  # force the conda route (into the active conda env) #
#   ./scripts/install_netcdf.sh --apt    # force the Debian/Ubuntu apt route                #
#   ./scripts/install_netcdf.sh -y       # install without prompting (for CI)              #
#   ./scripts/install_netcdf.sh -h       # help                                            #
#                                                                                          #
# Notes:                                                                                    #
#   * The apt route needs root; the script uses sudo and will prompt for your password.     #
#   * For the full output + post-processing stack at once, prefer: mamba env create -f      #
#     environment.yml  (installs libnetcdf plus the Python tools).                          #
#==========================================================================================#
set -euo pipefail

ASSUME_YES="${EDMORE_ASSUME_YES:-0}"
ROUTE="auto"               # auto | conda | apt

usage() {
   cat <<'EOF'
install_netcdf.sh — check for the netCDF C library and offer to install it.

Usage:
  ./scripts/install_netcdf.sh           detect; prompt before installing (conda if active, else apt)
  ./scripts/install_netcdf.sh --conda   force the conda route (into the active conda env)
  ./scripts/install_netcdf.sh --apt     force the Debian/Ubuntu apt route (libnetcdf-dev)
  ./scripts/install_netcdf.sh -y        install without prompting (CI)
  ./scripts/install_netcdf.sh -h        show this help

Environment:
  EDMORE_ASSUME_YES=1   same as -y

ED-MORE's CMake uses find_package(netCDF CONFIG); conda-forge's libnetcdf ships that config,
so the conda route is recommended. Build ED-MORE with -DCMAKE_PREFIX_PATH=<prefix> (the printed
prefix). netCDF is required for every build.
EOF
}

#----- Parse arguments. -------------------------------------------------------------------#
for arg in "$@"; do
   case "${arg}" in
      -y|--yes)  ASSUME_YES=1 ;;
      --conda)   ROUTE="conda" ;;
      --apt)     ROUTE="apt" ;;
      -h|--help) usage; exit 0 ;;
      *) echo "Unknown argument: ${arg} (try -h)"; exit 2 ;;
   esac
done

have_nc_config() { command -v nc-config >/dev/null 2>&1; }

confirm() {                                   # $1 = prompt; honor ASSUME_YES
   [[ "${ASSUME_YES}" == "1" ]] && return 0
   local reply
   if ! read -r -p "$1 [y/N] " reply; then reply=""; fi
   case "${reply}" in
      y|Y|yes|YES) return 0 ;;
      *) echo "Aborted — no changes made. Re-run when ready."; exit 0 ;;
   esac
}

report_netcdf() {                             # report version/prefix + how to point CMake at it
   local prefix version
   version="$(nc-config --version 2>/dev/null || true)"
   prefix="$(nc-config --prefix  2>/dev/null || true)"
   echo "✓ netCDF-C available: ${version:-version unknown}"
   [[ -n "${prefix}" ]] && echo "  prefix: ${prefix}"
   if [[ -n "${prefix}" && -f "${prefix}/lib/cmake/netCDF/netCDFConfig.cmake" ]]; then
      echo "  CMake config present -> find_package(netCDF CONFIG) will find it."
      echo "  Build ED-MORE with:  -DCMAKE_PREFIX_PATH=${prefix}"
   else
      echo "  NOTE: no CMake config under <prefix>/lib/cmake/netCDF — find_package(netCDF CONFIG)"
      echo "        may not locate it. conda-forge's libnetcdf ships one; install that (netCDF is"
      echo "        a required ED-MORE build dependency)."
   fi
}

#------------------------------------------------------------------------------------------#
# 1. Already available?                                                                     #
#------------------------------------------------------------------------------------------#
if have_nc_config; then
   report_netcdf
   exit 0
fi
echo "netCDF C library not found (no 'nc-config' on PATH)."

#------------------------------------------------------------------------------------------#
# 2. Pick a route when auto: conda if a conda toolchain is present, else apt.               #
#------------------------------------------------------------------------------------------#
if [[ "${ROUTE}" == "auto" ]]; then
   if command -v mamba >/dev/null 2>&1 || command -v conda >/dev/null 2>&1; then
      ROUTE="conda"
   else
      ROUTE="apt"
   fi
fi

#------------------------------------------------------------------------------------------#
# 3a. conda route — install libnetcdf from conda-forge into the active environment.         #
#------------------------------------------------------------------------------------------#
if [[ "${ROUTE}" == "conda" ]]; then
   CONDA_BIN=""
   command -v mamba >/dev/null 2>&1 && CONDA_BIN="mamba"
   [[ -z "${CONDA_BIN}" ]] && command -v conda >/dev/null 2>&1 && CONDA_BIN="conda"
   if [[ -z "${CONDA_BIN}" ]]; then
      echo "Requested the conda route, but neither 'mamba' nor 'conda' is on PATH." >&2
      echo "Install miniforge/conda first, or re-run with --apt." >&2
      exit 1
   fi
   if [[ -z "${CONDA_PREFIX:-}" ]]; then
      cat >&2 <<EOF
No conda environment is active (\$CONDA_PREFIX is empty). Activate one first, e.g.:
    conda activate base        # or your project env
then re-run, or create the full ED-MORE env in one step:
    mamba env create -f environment.yml && conda activate edmore
EOF
      exit 1
   fi
   echo "Will install 'libnetcdf' (conda-forge) into the active env: ${CONDA_PREFIX}"
   confirm "Run: ${CONDA_BIN} install -c conda-forge libnetcdf ?"
   echo "==> Installing libnetcdf via ${CONDA_BIN}"
   "${CONDA_BIN}" install -y -c conda-forge libnetcdf

#------------------------------------------------------------------------------------------#
# 3b. apt route — Debian/Ubuntu system package.                                             #
#------------------------------------------------------------------------------------------#
elif [[ "${ROUTE}" == "apt" ]]; then
   if ! command -v apt-get >/dev/null 2>&1; then
      cat >&2 <<EOF
The apt route supports Debian/Ubuntu (apt) only — e.g. WSL Ubuntu.
On macOS:  brew install netcdf
On Fedora/RHEL:  sudo dnf install netcdf-devel
Or use a conda environment (recommended): re-run with --conda, or
    mamba env create -f environment.yml
EOF
      exit 1
   fi
   echo "Will install 'libnetcdf-dev' (Debian/Ubuntu apt)."
   echo "(If find_package(netCDF CONFIG) later can't locate it, the conda route ships the CMake config.)"
   confirm "Install libnetcdf-dev via apt?"
   SUDO=""
   if [[ "${EUID}" -ne 0 ]]; then
      if command -v sudo >/dev/null 2>&1; then
         SUDO="sudo"
      else
         echo "Root privileges are required to install system packages, and sudo is absent." >&2
         exit 1
      fi
   fi
   echo "==> Updating package lists"
   ${SUDO} apt-get update
   echo "==> Installing libnetcdf-dev"
   ${SUDO} apt-get install -y libnetcdf-dev
fi

#------------------------------------------------------------------------------------------#
# 4. Verify.                                                                                #
#------------------------------------------------------------------------------------------#
echo
echo "==> Verifying installation"
if have_nc_config; then
   report_netcdf
else
   echo "netCDF installed, but 'nc-config' is not on PATH for this shell." >&2
   echo "(conda: re-activate the env; apt: open a new shell.)" >&2
fi

cat <<'EOF'

Done. Build ED-MORE (netCDF is required), pointing CMake at the prefix above:
    cmake -B build-dev -DCMAKE_Fortran_COMPILER=ifx -DCMAKE_PREFIX_PATH=<prefix>
    cmake --build build-dev -j
    LD_LIBRARY_PATH=<prefix>/lib ./build-dev/edmore ed_config_main.toml
EOF
