# Dependency Audit (Phase 1)

## Scope
This audit reconstructs dependencies for the legacy `VirHostMatcher-Net` Python project without refactoring core logic.

## Method
- Reviewed declared dependencies in `README.md` and `setup.py`.
- Scanned imports across Python modules under `src/`, `predictor.py`, and `VirHostMatcher-Net.py`.
- Identified runtime external executables invoked from Python.

## Reconstructed Direct Python Dependencies
- `numpy`
  - Used throughout (`predictor.py`, `src/s2star.py`, `setup.py` include path for C++ extension build).
- `pandas`
  - Used for all table I/O and feature/prediction matrix operations.
- `biopython`
  - Uses `Bio.Blast.Applications.NcbiblastnCommandline` in `src/crispr.py` (and `src/blast.py`).

## Version Constraints
- `numpy<2`
  - Required because code uses `np.float` (`predictor.py`), removed in NumPy 2.0.
- `pandas`
  - No strict pin required for import/runtime in this phase; installed latest compatible with interpreter.
- `biopython`
  - Pinned to `<1.85` because this code imports `Bio.Blast.Applications`, which is not available in newer Biopython releases used during validation.

## Non-Python / System Dependencies
- `BLAST+` (`blastn` executable required at runtime)
  - Invoked by `NcbiblastnCommandline` for CRISPR/BLAST feature generation.
- C++ toolchain (`g++`/compatible compiler)
  - Required to build local extension module `tools` from `setup.py`.

## Local Compiled Extension
- Module name: `tools`
- Sources:
  - `src/tools/tools.cpp`
  - `src/tools/WIsH.cpp`
  - `src/tools/mm.cpp`
- Include paths include NumPy headers and `src/tools/`.

## Data Dependency
- Input dataset tarball used for run validation:
  - `data_VirHostMatcher-Net_complete_genome_mode_alone.tar.gz`
- Expected extracted structure under `data/` includes at least:
  - `tables/`
  - `crispr_db_prefix/`
  - `host_wish_model/`

## Risks / Legacy Notes
- `numpy` 2.x is incompatible without code changes.
- `biopython` 1.85+ is incompatible without code changes (`Bio.Blast.Applications` import path).
- Runtime depends on external `blastn` being installed and on PATH.
- No formal unit tests are present; validation relies on execution of the CLI workflow.
