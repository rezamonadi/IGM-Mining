from pathlib import Path
import csv, json, datetime, re
import h5py, numpy as np
D=Path(__file__).resolve().parent
R=D.parent
read=lambda n:list(csv.DictReader((D/n).open()))
s=read('summary.csv');t=read('ratios_total.csv');o=read('ratio_offenders.csv');c=read('ratios_by_chunk.csv')
files=sorted(R.glob('processed_qsos_tst_posterior_corrected_S_*_E_*.mat'),key=lambda p:int(re.search(r'_S_(\d+)',p.name)[1]))
assert len(s)==len(files), 'Sweep incomplete: do not publish report.'
for a,b in zip(s,s[1:]):assert int(a['end'])+1==int(b['start']), 'Noncontiguous coverage'
failed=[r for r in s if r['pass_numerical']!='True']
with h5py.File(files[0]) as f:
    eligible=int(np.count_nonzero(f['test_ind'][:]))
def save_rows(name,rows,fields=None):
    fields=fields or list(rows[0])
    with (D/name).open('w',newline='') as fp:
        w=csv.DictWriter(fp,fieldnames=fields);w.writeheader();w.writerows(rows)
for threshold in [.5,.75,.85,.95,.99]:
    rows=[r for r in o if float(r['p_c4'])>=threshold]
    rows.sort(key=lambda r:(-int(r['ratio_lt_1_draws'])-int(r['ratio_gt_2_draws']),int(r['filtered_catalog_index']),int(r['component'])))
    save_rows('ratio_offenders_p'+str(threshold).replace('.','')+'.csv',rows,list(o[0]))
if not (D/'numerical_failures.csv').exists():
    save_rows('numerical_failures.csv',[],['start','result_qso_row','filtered_catalog_index','component','check'])
manifest=[]
for p in files:
    st=p.stat();manifest.append({'file':p.name,'size_bytes':st.st_size,'mtime_ns':st.st_mtime_ns})
(D/'source_manifest.json').write_text(json.dumps(manifest,indent=2)+'\n')
lines=[
'Posterior-corrected C IV chunk audit',
'Generated '+datetime.datetime.now().astimezone().isoformat(),
'',
'EXECUTION AND SCOPE',
'MATLAB launch was attempted with -batch, clean MATLAB_PREFDIR, and -nojvm/-nodisplay.',
'Each attempt exited with code 1 before executing commands or producing a MATLAB log.',
'Therefore these results are an independent Python/h5py audit, NOT a completed MATLAB diagnose_civ_single_run run.',
'The corrected MATLAB diagnostic and a sweep runner are supplied for execution in a working MATLAB environment.',
'',
f"Files: {len(s)}; saved QSO rows: {sum(int(r['qsos']) for r in s):,}; coverage: {s[0]['start']} through {s[-1]['end']} (contiguous, no overlaps).",
f"Processed QSOs: {sum(int(r['processed_qsos']) for r in s):,}; active QSO/component cells: {sum(int(r['active_cells']) for r in s):,}.",
f"Numerical passes: {len(s)-len(failed)}/{len(s)}. Posterior draws per active cell: {sorted(set(r['draws_per_cell'] for r in s))}.",
'',
f'Saved test_ind selects {eligible:,} eligible QSOs; audited chunks cover {100*sum(int(r["qsos"]) for r in s)/eligible:.2f}% of that saved selection.',
'',
'NUMERICAL CHECKS',
'Checks include saved shapes and filename row counts; complete QSO coverage; probability bounds and normalization;',
'component-prefix completion; joint-score algebra and probability reconstruction; finite MAP values and redshift bounds;',
'all five joint draw arrays and their 16th/84th percentile intervals; full L1 and L2 proposal-to-marginal likelihood reconstruction;',
'L2 rejection counts; shared test_ind identity across chunks; and unchanged source size/mtime during each audit.',
'Percentiles use NumPy method=hazen to match the MATLAB percentile convention used by the saved intervals.',
]
checks=['shape_failures','probability_invalid','prefix_violations','algebra_invalid_cells','map_invalid_cells','interval_invalid_cells','L1_likelihood_mismatches','L2_likelihood_mismatches','L1_invalid_sample_cells','L2_invalid_sample_cells','coverage_count_mismatches','negative_rew_pairs']
for k in checks:lines.append(f"  {k}: {sum(int(r[k]) for r in s):,}")
for k in ['L1_max_abs_error','L2_max_abs_error','probability_sum_max_error']:lines.append(f"  {k}: {max(float(r[k]) for r in s):.6g}")
lines += ['', 'DIAGNOSTIC CORRECTION',
'The original check removed -Inf proposals before taking the mean. Coverage-rejected proposals have zero likelihood,',
'and must remain in the original proposal denominator. The corrected check retains them and explicitly rejects NaN/+Inf.',
f"The old L2 check would flag {sum(int(r['L2_legacy_diagnostic_mismatches']) for r in s):,} cells; corrected reconstruction mismatches: {sum(int(r['L2_likelihood_mismatches']) for r in s):,}.",
f"Cells containing L2 coverage-rejected proposals: {sum(int(r['coverage_rejected_cells']) for r in s):,}.",
'Original and corrected diagnostic sources and diagnostic_fix.diff are preserved alongside this report.',
'', 'REW RATIOS',
'Ratios are W1548/W1550, with finite pairs and W1550 > 1e-8.',
'Flags use the existing diagnostic tolerances: ratio < 1-1e-5 or ratio > 2+1e-3.',
'Counts are posterior draws, not independent absorbers. Selection is p_c4 only, without a p_c4 > p_c4L1 cut.',
'Physical ratio warnings are separate from numerical integrity.',
'', 'p_c4 minimum | cells | ratio draws | below 1 | above 2 | outside | percent | affected cells | affected QSOs']
for r in t:
    lines.append(f"{r['p_c4_threshold']:>12} | {int(r['selected_cells']):,} | {int(r['ratio_pairs']):,} | {int(r['ratio_lt_1']):,} | {int(r['ratio_gt_2']):,} | {int(r['ratio_outside']):,} | {100*float(r['ratio_outside_fraction']):.6f}% | {int(r['affected_cells']):,} | {int(r['affected_qsos']):,}")
lines += ['', 'PER-CHUNK RESULTS', 'start-end | processed QSOs | active cells | numerical pass | all flagged ratio draws | p>=0.95 flags/draws | p>=0.99 flags/draws']
for r in s:
    rr={float(x['p_c4_threshold']):x for x in c if x['start']==r['start']}
    lines.append(f"{r['start']}-{r['end']} | {r['processed_qsos']} | {r['active_cells']} | {r['pass_numerical']} | {r['ratio_outside_draws']} | {rr[.95]['ratio_outside']}/{rr[.95]['ratio_pairs']} | {rr[.99]['ratio_outside']}/{rr[.99]['ratio_pairs']}")
hi=[r for r in o if float(r['p_c4'])>=.99];hi.sort(key=lambda r:-int(r['ratio_lt_1_draws'])-int(r['ratio_gt_2_draws']))
lines += ['', 'HIGH-CONFIDENCE FLAGGED CELLS (p_c4 >= 0.99)', 'filtered catalog index | chunk start | local row | component | p_c4 | below 1 draws | above 2 draws | ratio min | ratio max']
for r in hi:lines.append(' | '.join(r[k] for k in ['filtered_catalog_index','start','result_qso_row','component','p_c4','ratio_lt_1_draws','ratio_gt_2_draws','ratio_min','ratio_max']))
lines += ['', 'IDENTIFIERS AND LIMITS',
'filtered_catalog_index = chunk start + result_qso_row - 1. It indexes the filtered processing set, not the raw DR16Q catalog.',
'All posterior draws are retained in these counts; no QSO, including 3246, has been excluded.',
'No source MAT file was modified. No spectra were refitted, and ratio warnings have not been physically adjudicated.',
'Catalog-scale completeness beyond the available chunks is not established by this sweep.',
'', 'RUN THE ACTUAL MATLAB SWEEP',
"From the repository in MATLAB: addpath('posterior_corrected_diagnostics_20261102'); run_matlab_diagnostics",
'This creates per-chunk .matlab.log and .matlab_report.mat files plus matlab_summary.csv.',
'Independent audit outputs: summary.csv, ratios_by_chunk.csv, ratios_total.csv, ratio_offenders.csv,',
'ratio_offenders_p*.csv, numerical_failures.csv, source_manifest.json, and per-chunk JSON files.'
]
(D/'report.txt').write_text('\n'.join(lines)+'\n')
print('\n'.join(lines[:30]))
print('HIGH_CONFIDENCE_CELLS',len(hi))
