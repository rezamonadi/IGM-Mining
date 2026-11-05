# Branch comparison: `posterior-correction` vs `master`

## Comparison used

There is no local branch named `main`. The repository's local default branch is `master`, so this report compares the current `posterior-correction` branch with local `master` at `7de9942084482e7b28a137e8df543e4afea17842`.

The merge base is that same `master` commit. The branch contains two commits after it:

| Commit | Date | Subject |
|---|---|---|
| `f86ff0b567805b5f3e6b204006c5ab67021195d0` | 2026-11-05 | Added recent work after Ming-Feng's bug cath: repitted rows and many more bugs |
| `2d72d40a036a8a73e87d6d84208acfd325d52d3f` | 2026-11-05 | Added recent work after Ming-Feng's bug cath: repitted rows and many more bugs |

The local `master` is one commit ahead of `origin/master` (`7de9942` vs `835bb24`). This report uses local `master`, not a freshly fetched remote reference. The worktree also contains an unrelated modified `.gitignore` and untracked Markdown reports; these are not part of the branch commit comparison.

## Summary

Compared with local `master`, `posterior-correction` changes **114 tracked paths**, with **9,192,726 insertions and 329 deletions**. The large insertion count is dominated by captured `.out` files, not executable source. The branch adds a set of revised C IV processing and diagnostic tools, changes several existing DR16/DR7 scripts and run defaults, adds a DR7-specific path and transition-profile C/MEX sources, and checks in diagnostic plots and HPCC logs. It also removes one existing plot and one tracked entry named `CDDF-CIV-DR16`.

## Code and workflow changes

### DR16 C IV processing and posterior outputs

New processing variants and supporting code:

- `process_qsos_dr16_par_fixed.m` adds a separate fixed processing path. It stores likelihood information by absorber component, uses a fine Voigt grid and samples the convolved profile back onto native SDSS pixels, and saves compact, matched 256-draw posterior samples for absorber redshift, column density, width, and both equivalent widths.
- The fixed path adds a doublet coverage gate: L1 requires support around C IV 1548, while L2 requires support for both 1548 and 1550. The documented defaults require at least three surviving pixels within 300 km/s. The source explicitly describes this as a data-coverage criterion, not a cut on the measured REW ratio.
- `process_qsos_dr16_par_fixed_precoverage.m` provides a matched version with the coverage gate disabled for comparison. `process_qsos_dr16_par_nofine.m` provides a native-grid variant.
- `process_qsos_dr16_par.m` is also edited. Among the functional changes are processing an indexed QSO range, handling the final 250,001+ chunk, storing L1 likelihood samples with an absorber-component dimension, and adjusting MAP likelihood selection. It also comments out portions of its earlier REW and pixel-count calculations; the fixed processor is a separate file rather than a rename of this script.
- New orchestration and support files include `run_script_dr16_fixed.m`, `mergeProcessedDR16.m`, `preload_qsos_dr16_fixed.m`, `learn_qso_model_fixed.m`, `objective_fixed.m`, `mean_posterior_continuum.m`, and `reporter.m`.

The existing `set_parameters_dr16.m` changes the fine-grid setting from `nAVG = 20` to `nAVG = 8`, defines `fine_factor = nAVG + 1`, changes the sample count from 50,000 to 10,000, changes the testing-set name to `posterior_corrected`, and changes the special plate-2256 spectrum path. `run_script_dr16.m` changes its defaults to process 100 QSOs with plotting off and saving on, removes `clear`, changes a hard-coded MEX include/library invocation, and calls `mergeProcessedDR16` when merging is enabled.

`generate_c4_samples.m` changes the sigma offset to use the second quasi-random sequence dimension, changes the output location to `data/dr16/processed/`, and comments out the equivalent-width sampling and histogram/plot generation block. `build_catalog_dr16.m` now saves to the fixed path `data/dr16/processed/catalogDR16.mat`.

### DR7 and DR12 paths

- New DR7-specific files: `set_parameters_dr7_fixed.m`, `preload_qsos_dr7_fixed.m`, `read_spec_dr7_fixed.m`, and `plot_learned_covariance_dr7.m`.
- `run_script_dr7.m` is substantially reworked to select the fixed DR7 parameters, use fixed catalog/preload paths, and call the fixed preloader and learner. Its defaults enable catalog construction and disable preloading, learning, sampling, and processing. Most of its previous downstream processing steps are commented out.
- `build_catalog_dr7.m` now reads the DR7 no-BAL FITS catalog instead of a DR16 file, processes the full catalog instead of a 1,000-QSO subset, and saves to `data/dr7/processed/catalog.mat`.
- `process_qsos_dr12_par.m` adds a DR12 processing implementation.

These scripts contain hard-coded paths under `data/` and some MEX commands contain the HPC include/library path `/rhome/rmona003/bin/`. Those settings affect where the new workflows read/write data and where they expect dependencies.

### Voigt transition implementation

The branch adds `voigt_iP_transition.c` and `voigt_iP_fixed.c`, plus `test_voigt_iP_transition.m`. These provide separate transition-profile implementations and a MATLAB test/diagnostic script. The DR7 runner's MEX setup is changed to build `voigt_iP_fixed.c` when `voigtPrep` is enabled.

### Diagnostics and integrity auditing

- `diagnose_civ_pilot.m` and `diagnose_civ_single_run.m` add pilot and saved-output numerical checks.
- `hpcc_output_integrity_audit/audit_chunks.py` adds a Python/HDF5 audit for posterior-corrected chunk MAT files.
- `posterior_corrected_diagnostics_20261102/` adds an audit runner, report builder, corrected and original MATLAB diagnostic copies, a small patch file, and a MATLAB sweep launcher.
- The audit tooling checks stored array dimensions, posterior normalization and component-prefix structure, posterior-score algebra, MAP values and intervals, likelihood reconstruction, coverage rejection counts, REW-pair ratios, chunk identity, and whether files changed during an audit. The Python scripts require Python packages including `h5py` and `numpy`.
- The second commit (`2d72d40`) removes the dated per-chunk JSON summaries previously committed under `hpcc_output_integrity_audit/` and `posterior_corrected_diagnostics_20261102/`, along with `C13-full-hist2D.png`. Consequently, those JSON result snapshots are not present in the branch's final tree even though audit scripts remain.

The commit contents alone do not establish that the new code has been run successfully. No tests were run as part of this comparison.

## Data and diagnostic artifacts

### Captured HPCC output

The branch adds **54 `.out` files** under `CIV/miningDR16/outHPCC/`, totaling **9,186,945 inserted lines**. Most names follow `prc<start-index>.out` at five-thousand-row intervals from `prc1.out` through `prc240001.out`. Additional names are `prc245000.out`, `prc245001.out`, `prc250000.out`, `prc250001.out`, and `prc25000.out`; the last uses a different, shorter index format. These are captured outputs, not source files. Their presence alone does not verify that every associated processing job completed successfully.

### Plots

- Adds **20** diagnostic images under `plt-bad-DR/`, with filenames identifying QSO indices `230`, `393`, `492`, `710`, `1074`, `3187`, `3246`, `3387`, and `4616` (some indices have multiple component plots).
- Deletes `C13-full-hist2D.png`.

### Other path change

The `CIV/miningDR16/CDDF-CIV-DR16` submodule entry is deleted. On `master` it points to submodule commit `5c2ea7a5987380acd7c014e1a6ddaac05e9ffa78`; this branch no longer records that gitlink.

## Changed-file inventory

The source additions and edits are grouped below so a teammate can jump to the relevant workflow. Paths are relative to `CIV/miningDR16/`.

**New MATLAB files:**

`diagnose_civ_pilot.m`, `diagnose_civ_single_run.m`, `learn_qso_model_fixed.m`, `mean_posterior_continuum.m`, `mergeProcessedDR16.m`, `objective_fixed.m`, `plot_learned_covariance_dr7.m`, `pltQSOv1.m`, `pltQSOv2.m`, `preload_qsos_dr16_fixed.m`, `preload_qsos_dr7_fixed.m`, `process_qsos_dr12_par.m`, `process_qsos_dr16_par_fixed.m`, `process_qsos_dr16_par_fixed_precoverage.m`, `process_qsos_dr16_par_nofine.m`, `read_spec_dr7_fixed.m`, `reporter.m`, `run_script_dr16_fixed.m`, `set_parameters_dr7_fixed.m`, and `test_voigt_iP_transition.m`.

**Edited MATLAB files:**

`EW_dr16_voigt.m`, `build_catalog_dr16.m`, `build_catalog_dr7.m`, `generate_c4_samples.m`, `preload_qsos_dr16.m`, `process_qsos_dr16_par.m`, `run_script_dr16.m`, `run_script_dr7.m`, and `set_parameters_dr16.m`.

**New C/MEX and Python files:**

`voigt_iP_fixed.c`, `voigt_iP_transition.c`, `hpcc_output_integrity_audit/audit_chunks.py`, `posterior_corrected_diagnostics_20261102/audit_chunks.py`, and `posterior_corrected_diagnostics_20261102/build_report.py`.

**New diagnostic support files:**

`posterior_corrected_diagnostics_20261102/diagnose_civ_single_run.corrected.m`, `posterior_corrected_diagnostics_20261102/diagnose_civ_single_run.original.m`, `posterior_corrected_diagnostics_20261102/diagnostic_fix.diff`, and `posterior_corrected_diagnostics_20261102/run_matlab_diagnostics.m`.

## Review points

- The branch contains parallel versions of the DR16 processor (`process_qsos_dr16_par.m`, `..._fixed.m`, `..._fixed_precoverage.m`, and `..._nofine.m`). Check which entry point a run script selects before comparing output products.
- Several runners and catalog builders now use explicit `data/dr16/...` or `data/dr7/...` paths. Existing deployments that relied on `processed_directory(...)` may need corresponding files or path updates.
- The posterior-corrected audit JSON snapshots and the old histogram plot were removed by the latest commit; recover them from earlier branch history if needed.
- The 54 checked-in output logs account for almost all of the branch's line-count increase. Review whether these run artifacts belong in source control separately from the MATLAB and C changes.
- A whitespace scan (`git diff --check master...HEAD`) reports trailing whitespace and extra blank lines in several added or edited files, including captured logs. These were not changed as part of preparing this report.
- No runtime validation is claimed here. The comparison is based on commit ancestry, diff statistics, filenames, and source inspection.
