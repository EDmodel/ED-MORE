#!/usr/bin/env bash
# SPDX-License-Identifier: Apache-2.0
#==========================================================================================#
# install_ifx.sh — check for the Intel ifx Fortran compiler and offer to install it.       #
#                                                                                          #
# ifx is Intel's LLVM-based Fortran compiler, shipped in the Intel oneAPI toolkit.         #
# This script targets Debian/Ubuntu (apt) — e.g. WSL Ubuntu — and installs only the        #
# Fortran compiler package (intel-oneapi-compiler-fortran).                                 #
#                                                                                          #
# Usage:                                                                                    #
#   ./scripts/install_ifx.sh            # check; prompt before installing                  #
#   ./scripts/install_ifx.sh -y         # check; install without prompting (for CI)        #
#   MEDS_ASSUME_YES=1 ./scripts/install_ifx.sh   # same as -y                              #
#   ./scripts/install_ifx.sh -h         # help                                             #
#                                                                                          #
# Notes:                                                                                    #
#   * System install needs root; the script uses sudo and will prompt for your password.   #
#   * Intel MPI / MKL are separate oneAPI packages (intel-oneapi-mpi-devel, intel-oneapi-  #
#     mkl-devel); MEDS uses neither, so this script does not install them.                 #
#==========================================================================================#
set -euo pipefail

ONEAPI_ROOT="${ONEAPI_ROOT:-/opt/intel/oneapi}"
SETVARS="${ONEAPI_ROOT}/setvars.sh"
KEYRING="/usr/share/keyrings/oneapi-archive-keyring.gpg"
APT_LIST="/etc/apt/sources.list.d/oneAPI.list"
KEY_URL="https://apt.repos.intel.com/intel-gpg-keys/GPG-PUB-KEY-INTEL-SW-PRODUCTS.PUB"
FORTRAN_PKG="intel-oneapi-compiler-fortran"

ASSUME_YES="${MEDS_ASSUME_YES:-0}"

usage() {
   cat <<'EOF'
install_ifx.sh — check for the Intel ifx Fortran compiler and offer to install it.

Usage:
  ./scripts/install_ifx.sh          check; prompt before installing
  ./scripts/install_ifx.sh -y       check; install without prompting (CI)
  ./scripts/install_ifx.sh -h       show this help

Environment:
  MEDS_ASSUME_YES=1   same as -y
  ONEAPI_ROOT=<dir>   oneAPI install prefix (default /opt/intel/oneapi)

Targets Debian/Ubuntu (apt). Installs intel-oneapi-compiler-fortran via the
Intel oneAPI APT repository; the install step needs sudo (prompts for password).
EOF
}

#----- Parse arguments. -------------------------------------------------------------------#
for arg in "$@"; do
   case "${arg}" in
      -y|--yes) ASSUME_YES=1 ;;
      -h|--help) usage; exit 0 ;;
      *) echo "Unknown argument: ${arg} (try -h)"; exit 2 ;;
   esac
done

have_ifx() { command -v ifx >/dev/null 2>&1; }

# Source the oneAPI environment without tripping 'set -u' (setvars.sh references unset vars).
source_oneapi() {
   [[ -f "${SETVARS}" ]] || return 1
   set +u
   # shellcheck disable=SC1090
   source "${SETVARS}" >/dev/null 2>&1 || true
   set -u
}

report_ifx() {
   echo "✓ ifx available: $(command -v ifx)"
   ifx --version | head -1
}

#------------------------------------------------------------------------------------------#
# 1. Already on PATH?                                                                       #
#------------------------------------------------------------------------------------------#
if have_ifx; then
   report_ifx
   exit 0
fi

#------------------------------------------------------------------------------------------#
# 2. Installed but environment not sourced?                                                 #
#------------------------------------------------------------------------------------------#
if [[ -f "${SETVARS}" ]]; then
   echo "oneAPI is installed at ${ONEAPI_ROOT}, but ifx is not on your PATH."
   source_oneapi
   if have_ifx; then
      report_ifx
      echo
      echo "Activate it in every new shell by adding this line to ~/.bashrc:"
      echo "    source ${SETVARS}"
      exit 0
   fi
   echo "Could not activate the compiler from ${SETVARS}; reinstall may be needed." >&2
   exit 1
fi

#------------------------------------------------------------------------------------------#
# 3. Not installed — offer to install via the Intel oneAPI APT repository.                  #
#------------------------------------------------------------------------------------------#
echo "ifx (Intel Fortran compiler) is not installed."

if ! command -v apt-get >/dev/null 2>&1; then
   cat >&2 <<EOF
This installer supports Debian/Ubuntu (apt) only — e.g. WSL Ubuntu.
For other platforms, see Intel's instructions:
  https://www.intel.com/content/www/us/en/developer/tools/oneapi/fortran-compiler-download.html
EOF
   exit 1
fi

if [[ "${ASSUME_YES}" != "1" ]]; then
   if ! read -r -p "Install the Intel oneAPI Fortran compiler (ifx) now? [y/N] " reply; then
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

echo "==> Refreshing base packages and TLS prerequisites (wget, gpg, ca-certificates)"
${SUDO} apt-get update
${SUDO} apt-get install -y --no-install-recommends ca-certificates wget gnupg

echo "==> Adding the Intel oneAPI APT repository"
tmpkey="$(mktemp)"
trap 'rm -f "${tmpkey}"' EXIT
wget -qO "${tmpkey}" "${KEY_URL}"
gpg --dearmor < "${tmpkey}" | ${SUDO} tee "${KEYRING}" >/dev/null
echo "deb [signed-by=${KEYRING}] https://apt.repos.intel.com/oneapi all main" \
   | ${SUDO} tee "${APT_LIST}" >/dev/null

echo "==> Installing ${FORTRAN_PKG} (large download — this may take several minutes)"
${SUDO} apt-get update
${SUDO} apt-get install -y "${FORTRAN_PKG}"

echo
echo "==> Verifying installation"
source_oneapi
if have_ifx; then
   report_ifx
else
   echo "Installed, but ifx is not yet on PATH for this shell." >&2
fi

cat <<EOF

Done. To use ifx, activate the oneAPI environment in your shell:
    source ${SETVARS}
Add that line to ~/.bashrc to make ifx available in every new terminal.

To build MEDS with ifx via CMake (needs CMake >= 3.20; the dev kind, the others
are in docs/building.md, "Build kinds"):
    source ${SETVARS}
    cmake -B build-dev -DCMAKE_Fortran_COMPILER=ifx -DCMAKE_PREFIX_PATH=<netCDF prefix>
    cmake --build build-dev -j
EOF
