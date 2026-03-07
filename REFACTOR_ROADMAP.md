# Refactor Roadmap

## Scope and Constraints
- Preserve CLI compatibility exactly as shown in `README.md`.
- Prioritize high-ROI structural cleanup before deeper optimization.
- Do not require the large downloaded dataset for unit tests.
- Keep current BLAST feature state unchanged in initial phases (currently not active in scoring pipeline).
- Focus on safe, incremental refactors with runnable checkpoints.

## Current Baseline (What We Keep Stable)
- Main command and flags remain unchanged:
  - `python VirHostMatcher-Net.py -q ... -o ... [-t] [--short-contig] [-n] [-i] [-l] [-d]`
- Existing model/data behavior and output file formats remain unchanged during Phase 1-2.
- Native extension and BLAST/WIsH external tool usage remain functionally equivalent.

## Phase 1: Structural Cleanup (Highest ROI, Lowest Risk)
Goal: untangle architecture without changing user-facing behavior.

1. Introduce explicit runtime config object.
- Replace import-time argument side effects in `src/Variables.py`.
- Build a `Config` object from CLI args and pass it through functions/classes.

2. Split orchestration from script entrypoint.
- Keep `VirHostMatcher-Net.py` as compatibility wrapper.
- Move core run flow to a dedicated pipeline module (`run_prediction(config)` style).

3. Isolate validation and file-system setup.
- Move path checks and output/intermediate directory creation into helper modules.
- Keep existing error messages/exit behavior as close as possible.

4. Keep BLAST feature state unchanged.
- Do not re-enable BLAST scoring yet.
- Mark BLAST module as legacy/inactive in docs/comments to reduce confusion.

Acceptance criteria:
- Existing CLI command still works unchanged.
- Output directory structure and filenames remain the same.
- No regression in basic run completion.

## Phase 2: Readability and Maintainability Pass
Goal: make code easier to reason about and safer to modify.

1. Naming and API cleanup (non-breaking).
- Standardize function/variable naming style.
- Correct obvious typos like `s2star_caclculator` (keep compatibility shim if needed).

2. Remove dead/commented legacy blocks.
- Delete large commented-out sections that are superseded.
- Keep history in git; avoid retaining code tombstones in active files.

3. Add lightweight typing and docstrings.
- Add type hints to public functions and pipeline interfaces.
- Keep docs concise and practical.

4. Centralize constants.
- Separate static model constants from runtime path config.

Acceptance criteria:
- Module responsibilities are clear.
- Reduced file complexity and fewer side effects.
- No CLI behavior change.

## Phase 3: Test Foundation (Dataset-Independent)
Goal: protect refactors with fast tests that do not require large data packs.

1. Unit tests on pure logic.
- Scoring composition and score transform.
- Matrix/index alignment behavior.
- Neighborhood feature math on synthetic inputs.

2. Validation tests with tiny fixtures.
- Required-path validation logic.
- Error cases for malformed configs.

3. Optional integration smoke test (not default).
- Gate with env var/marker.
- Uses real dataset only when explicitly enabled.

Acceptance criteria:
- `pytest` default run works without downloaded dataset.
- Core mathematical/shape logic covered by small tests.

## Phase 4: Packaging and Developer UX
Goal: make builds/runs consistent on Linux/macOS and easier for users.

1. Add `pyproject.toml`.
- Standardize package/build metadata and extension build entry.
- Keep generated artifacts and workflows consistent with current behavior.

2. Standard commands.
- Add `Makefile` (or equivalent) with:
  - `make setup`
  - `make build-ext`
  - `make test`
  - `make smoke` (optional, data-dependent)

3. Documentation refresh.
- Update README setup section with concise Linux/macOS quick start.
- Clarify external dependencies (`blastn`, compiler) and expected data layout.

Acceptance criteria:
- New users can build and run with fewer manual steps.
- Linux/macOS setup is reproducible and documented.

## Phase 5: Performance and Reliability (After Structure Is Stable)
Goal: improve runtime and robustness with lower regression risk.

1. Profile hotspot stages.
- Measure `s2star`, CRISPR BLAST calls, and large dataframe operations.

2. Optimize targeted bottlenecks.
- Replace expensive concat/groupby patterns where safe.
- Add bounded retries/timeouts for external tool calls.

3. Improve observability.
- Add stage-level timing logs and clear failure diagnostics.

Acceptance criteria:
- Measurable runtime improvements on representative workloads.
- Better failure handling for external command stages.

## Proposed Target Layout (Incremental Migration)
```text
src/
  virhostmatcher/
    cli.py
    config.py
    pipeline.py
    predictors/
      host_predictor.py
      scoring.py
    features/
      s2star.py
      crispr.py
      wish.py
      neighborhood.py
    io/
      validators.py
      writers.py
      loaders.py
```

Notes:
- Keep compatibility entry script `VirHostMatcher-Net.py`.
- Migrate internals gradually; do not force big-bang rewrite.

## Recommended Execution Order
1. Phase 1 (config + pipeline split)  
2. Phase 2 (readability cleanup)  
3. Phase 3 (dataset-independent tests)  
4. Phase 4 (pyproject + standard commands)  
5. Phase 5 (performance tuning)

## Risk Controls
- Small PR-sized changes, one phase at a time.
- Preserve CLI and output schema until explicitly changed.
- Add/update tests before behavior-affecting edits.
- Keep fallbacks for legacy runtime dependencies during transition.
