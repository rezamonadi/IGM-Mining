for f=[1:5000:250001]
    if f ==245001
        continue;
    end
    fid = sprintf("processed_qsos_tst_posterior_corrected_S_%d_E_%d.mat",f, f+4999);
    report = diagnose_civ_single_run(fid);
    save(sprintf('rep-%d.mat', f), "report")
end