% Matched baseline for process_qsos_dr16_par_fixed.
%
% Run this script in place of process_qsos_dr16_par_fixed after the usual
% DR16 setup. It preserves the pre-coverage-gate likelihood while retaining
% all other current fixed-processor changes. Output files are tagged
% `_precoverage` and therefore cannot overwrite a coverage-gated run.

doublet_coverage_enabled = false;
process_qsos_dr16_par_fixed
