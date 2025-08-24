#!/usr/bin/env bash
set -euo pipefail

# Build and run the MSA program on Linux/macOS
# Usage:
#   ./run_msa.sh

cd "$(dirname "$0")"

if command -v g++ >/dev/null 2>&1; then
  CXX=g++
elif command -v clang++ >/dev/null 2>&1; then
  CXX=clang++
else
  echo "Error: No C++ compiler found (g++ or clang++)." >&2
  exit 1
fi

# Determine OpenMP flag support
OMP_FLAG="-fopenmp"
if [[ "$CXX" == "clang++" ]]; then
  # On some macOS/clang setups, -fopenmp requires libomp installed; attempt and fallback
  OMP_TEST=$(echo "#include <omp.h>\nint main(){return 0;}" | $CXX -x c++ - -o /dev/null -fopenmp 2>&1 || true)
  if [[ -n "$OMP_TEST" ]]; then
    echo "Warning: OpenMP not available for $CXX, proceeding without -fopenmp" >&2
    OMP_FLAG=""
  fi
fi

CXXFLAGS=(-O3 -DNDEBUG -std=c++17)
if [[ -n "$OMP_FLAG" ]]; then CXXFLAGS+=("$OMP_FLAG"); fi

echo "Compiling with $CXX ${CXXFLAGS[*]}..."
"$CXX" "${CXXFLAGS[@]}" main_final_syntax_ok.cpp blosum62.cpp -o msa_program
echo "Compilation successful -> ./msa_program"

if [[ ! -f input.fasta ]]; then
  echo "Warning: input.fasta not found in project root. The program may error if the file is missing." >&2
fi

# Set OpenMP defaults if not already set
: "${OMP_NUM_THREADS:=$(getconf _NPROCESSORS_ONLN 2>/dev/null || sysctl -n hw.ncpu 2>/dev/null || echo 1)}"
: "${OMP_PROC_BIND:=close}"
: "${OMP_SCHEDULE:=static}"
export OMP_NUM_THREADS OMP_PROC_BIND OMP_SCHEDULE

echo "Running ./msa_program..."
./msa_program

