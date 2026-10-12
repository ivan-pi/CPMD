#!/usr/bin/env bash
#
# Evaluate a CPMD machine configuration file (configure/<NAME>) the same way
# configure.sh does, and print the resulting settings as CMake code
# (one "set(CPMD_CFG_<VAR> ...)" per variable).
#
# The CMake build calls this script automatically when CPMD_CONFIG is given,
# either on the command line or in the environment:
#
#   cmake -B build -DCPMD_CONFIG=LINUX-X86_64-GFORTRAN-MPI [-DCPMD_ENABLE_OPENMP=ON]
#   CPMD_CONFIG=LINUX-X86_64-GFORTRAN-MPI cmake -B build
#
# It can also be run by hand to see what CMake will get from a configuration.
#
# Usage: import_config.sh [-omp] [-debug] [-coverage] [-qmmm] CONFIG
#
#   CONFIG     name of a file in the configure/ directory, or a path to a file
#   -omp, -debug, -coverage, -qmmm
#              same meaning as for configure.sh
#

CPMD_ROOT=$(cd "$(dirname "${BASH_SOURCE[0]}")/../.." && pwd)

unset omp debug coverage qmmm
Configuration=
for option in "$@"; do
  case $option in
    -omp)           omp=1 ;;
    -debug)         debug=1 ;;
    -coverage)      coverage=1 ;;
    -qmmm)          qmmm=1 ;;
    -*)             echo "import_config.sh: unknown option $option" >&2; exit 2 ;;
    *)              Configuration=$option ;;
  esac
done

if [ -z "$Configuration" ]; then
  echo "usage: import_config.sh [-omp] [-debug] [-coverage] [-qmmm] CONFIG" >&2
  exit 2
fi
if [ -f "$Configuration" ]; then
  cfg_file=$Configuration
elif [ -f "${CPMD_ROOT}/configure/${Configuration}" ]; then
  cfg_file=${CPMD_ROOT}/configure/${Configuration}
else
  echo "import_config.sh: unknown configuration '${Configuration}'" >&2
  exit 2
fi

# Helpers that configuration files may use (see configure.sh)
if [ -x /usr/bin/nawk ]; then AWK='/usr/bin/nawk'; else AWK='awk'; fi
GREP='grep'

# Variables used by the CMake build (configuration files may set others,
# such as LD, CPP, LIBS, AR or RANLIB, which CMake does not need)
VARS="FC CC CPPFLAGS FFLAGS CFLAGS LFLAGS NVCC NVCCFLAGS NOOPT_FLAG NOOPT_OBJS
      FFLAGS_GROMOS FFLAGS_GROMOS_MODULES OMP3_DISABLED"
for v in $VARS; do unset "$v"; done

# Source the configuration. Anything it prints goes to stderr; an 'exit'
# inside the file (e.g. "Non-OMP settings not available!") aborts the import.
exec 3>&1 1>&2
. "$cfg_file"
exec 1>&3 3>&-

# Resolve make-style references such as $(MKLROOT) from the environment.
# $(FFLAGS) inside the *_GROMOS flags is dropped: in the CMake build those
# flags are appended to the regular Fortran flags.
expand() {
  local name=$1 val=$2 ref
  if [ "${name#FFLAGS_GROMOS}" != "$name" ]; then
    val=${val//\$(FFLAGS)/}
  fi
  while [[ $val =~ \$\(([A-Za-z_][A-Za-z0-9_]*)\) ]]; do
    ref=${BASH_REMATCH[1]}
    val=${val//\$($ref)/${!ref}}
  done
  # collapse line continuations and whitespace (without globbing)
  val=${val//\\$'\n'/ }
  set -f
  echo $val
  set +f
}

for v in $VARS; do
  [ -n "${!v+x}" ] || continue
  val=$(expand "$v" "${!v}")
  printf 'set(CPMD_CFG_%s [==[%s]==])\n' "$v" "$val"
done
printf 'set(CPMD_CFG_FILE [==[%s]==])\n' "$cfg_file"
