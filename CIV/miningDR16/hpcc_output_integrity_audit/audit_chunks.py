"""Independent HDF5 audit mirroring diagnose_civ_single_run; no MATLAB execution."""
from pathlib import Path
import csv, json, hashlib, re, time
import numpy as np
import h5py
ROOT=Path(__file__).resolve().parent.parent
OUT=Path(__file__).resolve().parent
THRESH=[0,.5,.75,.85,.95,.99]
files=sorted(ROOT.glob('processed_qsos_tst_posterior_corrected_S_*_E_*.mat'),key=lambda p:int(re.search(r'_S_(\d+)',p.name)[1]))
summary=[]; thresholds=[]; offenders=[]; failures=[]; reference_test=None

def write_csv(name,rows):
    if rows:
        with (OUT/name).open('w',newline='') as fp:
            w=csv.DictWriter(fp,fieldnames=list(rows[0]));w.writeheader();w.writerows(rows)

def record_bad(mask,label,start):
    for c,q in np.argwhere(mask):
        failures.append(dict(start=start,result_qso_row=int(q+1),filtered_catalog_index=int(start+q),component=int(c+1),check=label))

for path in files:
    ts=time.time(); st=path.stat(); start,end=map(int,re.search(r'_S_(\d+)_E_(\d+)',path.name).groups())
    print('START',path.name,flush=True)
    try:
        with h5py.File(path,'r') as f:
            p=np.stack([f[k][:] for k in ['p_no_c4','p_c4L1','p_c4']]); nc,nq=p.shape[1:]; active=np.all(np.isfinite(p),axis=0); any_active=np.any(np.isfinite(p),axis=0)
            s=dict(file=path.name,start=start,end=end,qsos=nq,processed_qsos=int(np.count_nonzero(active.any(axis=0))),active_cells=int(active.sum()),draws_per_cell=int(f['num_posterior_draws'][0,0]),active_by_component=' '.join(map(str,active.sum(axis=1))),shape_failures=0)
            nd=s['draws_per_cell']
            model_keys=['no_c4','c4L1','c4L2']
            for kind in ['log_likelihoods_','log_posteriors_']:
                for k in model_keys:
                    if f[kind+k].shape!=(nc,nq): s['shape_failures']+=1
            if nq!=end-start+1:s['shape_failures']+=1
            s['probability_invalid']=int(np.count_nonzero(any_active&~active)+np.count_nonzero(active&((p.min(axis=0)<-1e-12)|(p.max(axis=0)>1+1e-12)|(np.abs(p.sum(axis=0)-1)>1e-10))))
            s['prefix_violations']=int(np.count_nonzero(np.any(active[1:]&~active[:-1],axis=0)))
            post=np.stack([f['log_posteriors_'+k][:] for k in model_keys]); ll=np.stack([f['log_likelihoods_'+k][:] for k in model_keys]); pri=np.stack([f['log_priors_no_c4'][:],f['log_priors_c4'][:],f['log_priors_c4'][:]])
            joint=ll+pri
            with np.errstate(invalid='ignore',over='ignore'):
                e=np.exp(joint-joint.max(axis=0)); prob=e/e.sum(axis=0)
            algebra=active&(~np.all(np.isfinite(joint)&np.isfinite(post),axis=0)|np.any(np.abs(joint-post)>1e-9*np.maximum(1,np.maximum(np.abs(joint),np.abs(post))),axis=0)|np.any(np.abs(prob-p)>1e-9,axis=0))
            s['algebra_invalid_cells']=int(algebra.sum());record_bad(algebra,'algebra',start)
            s['probability_sum_max_error']=float(np.max(np.abs(p.sum(axis=0)[active]-1)))
            zmin=f['min_z_c4s'][:].reshape(-1);zmax=f['max_z_c4s'][:].reshape(-1)
            zmap=f['map_z_c4L2'][:];nmap=f['map_N_c4L2'][:];smap=f['map_sigma_c4L2'][:]
            badmap=active&(~np.isfinite(zmap)|~np.isfinite(nmap)|~np.isfinite(smap)|(smap<=0)|(zmap<zmin)|(zmap>zmax))
            s['map_invalid_cells']=int(badmap.sum());record_bad(badmap,'MAP',start)
            s['interval_invalid_cells']=0
            for name in ['z_c4L2','logN_c4L2','sigma_c4L2','REW_1548','REW_1550']:
                d=f['posterior_'+name+'_draws'][:];lo=f['posterior_'+name+'_16'][:];hi=f['posterior_'+name+'_84'][:]
                if d.shape!=(nc,nd,nq) or lo.shape!=(nc,nq) or hi.shape!=(nc,nq):s['shape_failures']+=1
                pct=np.quantile(d,[.16,.84],axis=1,method='hazen')
                tol_lo=1e-9*np.maximum(1,np.maximum(np.abs(lo),np.abs(pct[0])));tol_hi=1e-9*np.maximum(1,np.maximum(np.abs(hi),np.abs(pct[1])))
                bad=active&(~np.all(np.isfinite(d),axis=1)|~np.isfinite(lo)|~np.isfinite(hi)|(lo>hi)|(np.abs(lo-pct[0])>tol_lo)|(np.abs(hi-pct[1])>tol_hi))
                s['interval_invalid_cells']+=int(bad.sum());record_bad(bad,'interval_'+name,start)
                del d,pct
            for model in ['L1','L2']:
                # Read the cube once; HDF5 chunks span all components and 1000 QSOs.
                cube=f['sample_log_likelihoods_c4'+model][:]; stored=f['log_likelihoods_c4'+model][:]
                bad_count=0; legacy_count=0; maxerr=0.; invalid_samples=0; rejection_mismatch=0; rejected_cells=0
                for c in range(nc):
                    qidx=np.flatnonzero(active[c]);x=cube[c][:,qidx]
                    invalid=np.any(np.isnan(x)|np.isposinf(x),axis=0);invalid_samples+=int(invalid.sum())
                    mx=np.max(x,axis=0)
                    with np.errstate(invalid='ignore',divide='ignore'):
                        total=np.exp(x-mx).sum(axis=0);expected=mx+np.log(total/x.shape[0]);legacy=mx+np.log(total/np.isfinite(x).sum(axis=0))
                    actual=stored[c,qidx];delta=np.abs(actual-expected);tol=1e-9*np.maximum(1,np.maximum(np.abs(actual),np.abs(expected)))
                    bad=invalid|~np.isfinite(expected)|(delta>tol);bad_count+=int(bad.sum())
                    legacy_count+=int(np.count_nonzero(np.abs(actual-legacy)>1e-9*np.maximum(1,np.maximum(np.abs(actual),np.abs(legacy)))))
                    if len(delta):maxerr=max(maxerr,float(np.nanmax(delta)))
                    for q in qidx[bad]:failures.append(dict(start=start,result_qso_row=int(q+1),filtered_catalog_index=int(start+q),component=c+1,check='likelihood_'+model))
                    if model=='L2':
                        cnt=np.isneginf(x).sum(axis=0);rejected_cells+=int(np.count_nonzero(cnt))
                        rejection_mismatch+=int(np.count_nonzero(cnt!=f['num_doublet_coverage_rejected'][c,:][qidx]))
                    del x
                s[model+'_likelihood_mismatches']=bad_count;s[model+'_max_abs_error']=maxerr;s[model+'_invalid_sample_cells']=invalid_samples;s[model+'_legacy_diagnostic_mismatches']=legacy_count
                if model=='L2':s['coverage_rejected_cells']=rejected_cells;s['coverage_count_mismatches']=rejection_mismatch
                del cube
            a=f['posterior_REW_1548_draws'][:];b=f['posterior_REW_1550_draws'][:]
            good=active[:,None,:]&np.isfinite(a)&np.isfinite(b); valid=good&(b>1e-8)
            ratio=np.divide(a,b,out=np.full_like(a,np.nan),where=b>1e-8);low=valid&(ratio<1-1e-5);high=valid&(ratio>2+1e-3)
            s['finite_rew_pairs']=int(good.sum());s['negative_rew_pairs']=int(np.count_nonzero(good&((a<-1e-10)|(b<-1e-10))));s['ratio_outside_draws']=int(np.count_nonzero(low|high))
            s['coverage_enabled']=bool(f['doublet_coverage_enabled'][0,0]);s['coverage_min_pixels']=int(f['doublet_coverage_min_pixels'][0,0]);s['coverage_velocity_kms']=float(f['doublet_coverage_velocity_kms'][0,0])
            for t in THRESH:
                cell=active&(p[2]>=t);sel=valid&cell[:,None,:];out=(low|high)&sel
                thresholds.append(dict(start=start,end=end,p_c4_threshold=t,selected_cells=int(cell.sum()),ratio_pairs=int(sel.sum()),ratio_lt_1=int((low&sel).sum()),ratio_gt_2=int((high&sel).sum()),ratio_outside=int(out.sum()),affected_cells=int(np.count_nonzero(out.any(axis=1))),affected_qsos=int(np.count_nonzero(out.any(axis=(0,1))))))
            for c,q in np.argwhere((low|high).any(axis=1)):
                rr=ratio[c,:,q];offenders.append(dict(start=start,result_qso_row=int(q+1),filtered_catalog_index=int(start+q),component=int(c+1),p_c4=float(p[2,c,q]),p_c4L1=float(p[1,c,q]),ratio_lt_1_draws=int(low[c,:,q].sum()),ratio_gt_2_draws=int(high[c,:,q].sum()),ratio_min=float(np.nanmin(rr)),ratio_max=float(np.nanmax(rr))))
            digest=hashlib.sha256(f['test_ind'][:].tobytes()).hexdigest()
            if reference_test is None:reference_test=digest
            s['test_ind_matches_first']=digest==reference_test
            after=path.stat();s['unchanged_during_audit']=st.st_size==after.st_size and st.st_mtime_ns==after.st_mtime_ns
            checks=['shape_failures','probability_invalid','prefix_violations','algebra_invalid_cells','map_invalid_cells','interval_invalid_cells','L1_likelihood_mismatches','L2_likelihood_mismatches','coverage_count_mismatches']
            s['pass_numerical']=all(s[k]==0 for k in checks) and s['processed_qsos']==nq and s['unchanged_during_audit'] and s['test_ind_matches_first']
            s['elapsed_seconds']=round(time.time()-ts,2);summary.append(s)
            (OUT/(path.stem+'.json')).write_text(json.dumps(s,indent=2)+'\n')
            print('DONE',start,'pass',s['pass_numerical'],'ratio flags',s['ratio_outside_draws'],'seconds',s['elapsed_seconds'],flush=True)
    except Exception as e:
        print('ERROR',path.name,repr(e),flush=True)
        (OUT/(path.stem+'.error.txt')).write_text(repr(e)+'\n')
    write_csv('summary.csv',summary);write_csv('ratios_by_chunk.csv',thresholds);write_csv('ratio_offenders.csv',offenders);write_csv('numerical_failures.csv',failures)

aggregated=[]
for t in THRESH:
    rows=[r for r in thresholds if r['p_c4_threshold']==t]
    row={'p_c4_threshold':t}
    for k in ['selected_cells','ratio_pairs','ratio_lt_1','ratio_gt_2','ratio_outside','affected_cells','affected_qsos']:row[k]=sum(r[k] for r in rows)
    row['ratio_outside_fraction']=row['ratio_outside']/max(1,row['ratio_pairs']);aggregated.append(row)
write_csv('ratios_total.csv',aggregated)
print('TOTAL',len(summary),'of',len(files),'files','numerical passes',sum(s['pass_numerical'] for s in summary),flush=True)
for row in aggregated:print(row,flush=True)
