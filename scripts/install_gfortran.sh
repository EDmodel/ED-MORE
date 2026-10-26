#!/usr/bin/env bash
# SPDX-License-Identifier: Apache-2.0
#==========================================================================================#
# install_gfortran.sh — check for the GNU gfortran Fortran compiler and offer to install it.#
#                                                                                          #
# gfortran is the GNU Fortran compiler (part of GCC) and ED2's reference toolchain.        #
# This script targets Debian/Ubuntu (apt) — e.g. WSL Ubuntu — and installs gfortran from   #
# the standard distribution repositories (no third-party repo needed).                     #
#                                                                                          #
# Usage:                                                                                    #
#   ./scripts/install_gfortran.sh         # check; prompt before installing gfortran       #
#   ./scripts/install_gfortran.sh --hdf5  # also install HDF5 dev files (libhdf5-dev)       #
#   ./scripts/install_gfortran.sh -y      # install without prompting (for CI)             #
#   ./scripts/install_gfortran.sh -h      # help                                           #
#                                                                                          #
# Notes:                                                                                    #
#   * System install needs root; the script uses sudo and will prompt for your password.   #
#   * Ubuntu's HDF5 Fortran bindings land under /usr/include/hdf5/serial. ED-MORE needs    #
#     none (its I/O is netCDF); a build that does points there or uses h5fc/pkg-config.    #
#==========================================================================================#
set -euo pipefail

ASSUME_YES="${EDMORE_ASSUME_YES:-0}"
WITH_HDF5=0
PACKAGES=(gfortran)

usage() {
   cat <<'EOF'
install_gfortran.sh — check for the GNU gfortran compiler and offer to install it.

Usage:
  ./scripts/install_gfortran.sh           check; prompt before installing
  ./scripts/install_gfortran.sh --hdf5    also install HDF5 dev files (libhdf5-dev)
  ./scripts/install_gfortran.sh -y        check; install without prompting (CI)
  ./scripts/install_gfortran.sh -h        show this help

Environment:
  EDMORE_ASSUME_YES=1   same as -y

Targets Debian/Ubuntu (apt). Installs gfortran (and optionally libhdf5-dev) from
the standard repositories; the install step needs sudo (prompts for password).
EOF
}

#----- Parse arguments. -------------------------------------------------------------------#
for arg in "$@"; do
   case "${arg}" in
      -y|--yes)  ASSUME_YES=1 ;;
      --hdf5)    WITH_HDF5=1; PACKAGES+=(libhdf5-dev) ;;
      -h|--help) usage; exit 0 ;;
      *) echo "Unknown argument: ${arg} (try -h)"; exit 2 ;;
   esac
done

have_gfortran() { command -v gfortran >/dev/null 2>&1; }

report_gfortran() {
   echo "✓ gfortran available: $(command -v gfortran)"
   gfortran --version | head -1
}

#------------------------------------------------------------------------------------------#
# 1. Already on PATH?  (Still install HDF5 if requested and absent.)                        #
#------------------------------------------------------------------------------------------#
if have_gfortran && [[ "${WITH_HDF5}" -eq 0 ]]; then
   report_gfortran
   exit 0
fi

if have_gfortran; then
   echo "gfortran is already installed; proceeding to install HDF5 dev files as requested."
   PACKAGES=(libhdf5-dev)
else
   echo "gfortran (GNU Fortran compiler) is not installed."
fi

if ! command -v apt-get >/dev/null 2>&1; then
   cat >&2 <<EOF
This installer supports Debian/Ubuntu (apt) only — e.g. WSL Ubuntu.
On other systems install the GNU Fortran compiler with your package manager
(macOS: 'brew install gcc'; Fedora/RHEL: 'sudo dnf install gcc-gfortran').
EOF
   exit 1
fi

if [[ "${ASSUME_YES}" != "1" ]]; then
   if ! read -r -p "Install: ${PACKAGES[*]} ? [y/N] " reply; then
      reply=""
   fi
   case "${reply}" in
      y|Y|yes|YES) ;;
      *) echo "Aborted — no changes made. Re-run when ready."; exit 0 ;;
   esac
fi

#----- Resolve sudo. ----------------------------------------------------------------------#
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
echo "==> Installing: ${PACKAGES[*]}"
${SUDO} apt-get install -y "${PACKAGES[@]}"

echo
echo "==> Verifying installation"
if have_gfortran; then
   report_gfortran
else
   echo "gfortran still not found on PATH — installation may have failed." >&2
   exit 1
fi

if [[ "${WITH_HDF5}" -eq 1 ]]; then
   echo "✓ HDF5 dev files installed (Ubuntu: headers/.mod under /usr/include/hdf5/serial,"
   echo "  libs under /usr/lib/x86_64-linux-gnu/hdf5/serial)."
fi

cat <<'EOF'

Done. gfortran is ready. To build ED-MORE with it via CMake (the dev kind; the
others are in docs/building.md, "Build kinds"):
    cmake -B build-dev -DCMAKE_Fortran_COMPILER=gfortran -DCMAKE_PREFIX_PATH=<netCDF prefix>
    cmake --build build-dev -j

Ubuntu 20.04 ships gfortran 9 (good Fortran 2018 coverage). For newer standards
support you can install a later series, e.g.: sudo apt-get install gfortran-11
EOF
