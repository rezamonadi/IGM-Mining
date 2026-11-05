function run_matlab_diagnostics
% Run the actual corrected MATLAB diagnostic on every root-level chunk.
% From the repository: addpath('posterior_corrected_diagnostics_20261102'); run_matlab_diagnostics
outdir = fileparts(mfilename('fullpath'));
root = fileparts(outdir);
addpath(root);
files = dir(fullfile(root, 'processed_qsos_tst_posterior_corrected_S_*_E_*.mat'));
starts = arrayfun(@(f) str2double(regexp(f.name, '(?<=_S_)\d+', 'match', 'once')), files);
[~,order] = sort(starts); files = files(order);
summary = table;
for i = 1:numel(files)
    name = files(i).name;
    [~,stem] = fileparts(name);
    logfile = fullfile(outdir, [stem '.matlab.log']);
    diary(logfile);
    fprintf('\nMATLAB diagnostic %d/%d: %s\n', i, numel(files), name);
    try
        report = diagnose_civ_single_run(fullfile(root,name), ...
            'MakePlots',false, 'CheckSampleLikelihoods',true, 'BlockQSOs',1000);
        save(fullfile(outdir,[stem '.matlab_report.mat']), 'report');
        writetable(report.rew.by_p_c4, fullfile(outdir,[stem '.matlab_ratios.csv']));
        row = table(string(name), report.shape(1), report.processed_qsos, ...
            report.active_cells, report.pass_numerical, report.pass_rew_ratio, ...
            report.sample_likelihoods.mismatch, string(''), ...
            'VariableNames', {'file','qsos','processed_qsos','active_cells', ...
            'pass_numerical','pass_rew_ratio','likelihood_mismatches','error'});
    catch err
        fprintf('%s\n',getReport(err,'extended','hyperlinks','off'));
        row = table(string(name), NaN, NaN, NaN, false, false, NaN, ...
            string(err.message), 'VariableNames', {'file','qsos','processed_qsos', ...
            'active_cells','pass_numerical','pass_rew_ratio','likelihood_mismatches','error'});
    end
    diary off;
    summary = [summary;row]; %#ok<AGROW>
    writetable(summary,fullfile(outdir,'matlab_summary.csv'));
end
disp(summary);
end
