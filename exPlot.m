variables_to_load = {'all_wavelengths', 'all_flux', 'all_noise_variance', ...
     'all_pixel_mask', 'all_sigma_pixel'};
load(sprintf('%s/preloaded_qsos_%s', processed_directory('dr12'), ...
        training_set_name), variables_to_load{:});



pltQSO(this_flux, this_wavelengths, c4_muL2, c4_muL1, ind_zoomL2, ind_zoomL1, ttl, fid)
