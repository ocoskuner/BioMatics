## 2025-08-24

### Added
- CHANGELOG.md to document performance-focused, behavior-preserving refinements.

### Optimizations (behavior preserved)
- BLOSUM matrix caching:
  - `blosum62.cpp`: `getRealBlosumLogOddsMatrix()` now computes the matrix once and returns a cached copy on subsequent calls. Thread-safe static initialization; numerical values unchanged.

- One-hot EMD lookup for pairwise distances:
  - `main_final_syntax_ok.cpp`: For raw sequence distance matrix, per-position distributions are one-hot. The EMD of two one-hot distributions equals `tanh(clamp(cost[i][j]))` with the same clamping as the original path. A precomputed `emd_lookup` table is used to replace repeated EMD calls, producing identical results with less compute.
  - Unknown characters (not in `AA_INDEX`) still contribute `EMD_ERROR` to the average, preserving prior behavior.

- Parallelization of pairwise distances:
  - `main_final_syntax_ok.cpp`: Built a flat list of `(i,j)` pairs with `i < j` and processed them using OpenMP `parallel for`. Each pair writes only to `D[i][j]` and `D[j][i]`, eliminating data races. The inner per-position accumulation remains sequential to preserve bitwise-identical results.

- Global distance reuse in iterative refinement:
  - `main_final_syntax_ok.cpp`: The full distance matrix computed once at startup is cached in `GLOBAL_DIST_MATRIX`. During `iterativeRefinement`, submatrices for rest-sets are sliced from the cache instead of recomputing distances, with the same CSV side-effect kept via `writeCSV` for consistency.

### Build and scripts
- PowerShell (Windows):
  - `run_msa.ps1` now compiles with `-O3 -DNDEBUG -std=c++17 -fopenmp` and sets OpenMP runtime defaults if absent (`OMP_NUM_THREADS`=logical cores, `OMP_PROC_BIND=close`, `OMP_SCHEDULE=static`). Warns if `input.fasta` is missing.

- Bash (Linux/macOS):
  - `run_msa.sh` now compiles with `-O3 -DNDEBUG -std=c++17` and attempts `-fopenmp` when available (detects clang OpenMP support and falls back cleanly). Sets the same OpenMP defaults when running.

### Notes on correctness
- All changes were designed to preserve outputs exactly:
  - BLOSUM caching affects only performance, not values.
  - One-hot lookup is mathematically equivalent to the previous EMD path for one-hot inputs, using the same clamp and `tanh` transformation.
  - Parallelization processes independent pairs; the per-pair reduction order is unchanged by keeping inner loops sequential.
  - Refinement now reuses the originally computed distances; no recomputation with different code paths, so values remain the same.

### Future optional improvements (also behavior-preserving)
- Parallelize column-based scoring (`computeCSScore`, `computeSPSScore`, `findStrongColumns`).
- Parallelize the minimum-pair scan in guide tree construction with a parallel reduction.
- Reduce debug log flushing frequency inside tight loops or guard with a verbosity flag.


