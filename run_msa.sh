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

echo "Compiling with $CXX (-std=c++17 -g)..."
"$CXX" -std=c++17 -g -fopenmp main_final_syntax_ok.cpp blosum62.cpp -o msa_program
echo "Compilation successful -> ./msa_program"

if [[ ! -f input.fasta ]]; then
  echo "Warning: input.fasta not found in project root. The program may error if the file is missing." >&2
fi

echo "Running ./msa_program..."
./msa_program

