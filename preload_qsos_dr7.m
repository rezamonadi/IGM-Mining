% preload_qsos: loads spectra from SDSS FITS files, applies further
% filters, and applies some basic preprocessing such as normalization
% and truncation to the region of interest

% load QSO catalog
tic;

num_quasars = numel(all_zqso_dr7);
all_wavelengths    =  cell(num_quasars, 1);
all_flux           =  cell(num_quasars, 1);
all_noise_variance =  cell(num_quasars, 1);
all_pixel_mask     =  cell(num_quasars, 1);
all_sigma_pixel     =  cell(num_quasars, 1);
all_normalizers    = zeros(num_quasars, 1);


% to track reasons for filtering out QSOs
filter_flags = zeros(num_quasars, 1, 'uint8');
% 
% filtering bit 0: z_QSO < 1.7
ind = (all_zqso_dr7 < z_qso_cut);
filter_flags(ind) = bitset(filter_flags(ind), 1, true);



for i = 1:num_quasars


  if (filter_flags(i)~=0)
    continue;
  end
  
  %------dr7-----------
  [this_wavelengths, this_flux, this_noise_variance, this_pixel_mask, this_sigma_pixel] ...
      = file_loader(all_mjd_dr7(i), all_plate_dr7(i),all_fiber_dr7(i));


  
  this_sigma_pixel = this_sigma_pixel';
	% Here file_loader uses dr7 spectrum reader function and given mpf to read
	% spectrum 
  
  % bit 1: the spectrum with median SN<4 should be filtered
  if (nanmedian(this_flux./sqrt(this_noise_variance))<4)
    filter_flags(i) = bitset(filter_flags(i), 2, true);
    continue;
  end
  % CIV removal region
  % if (civ_rmove_in_train =true)
    
  %   maskCIV = 
  %   this_pixel_mask()
  % % Masking Sky lines 
  if SkyLine==1
    this_pixel_mask((abs(this_wavelengths-5579)<5) & (abs(this_wavelengths-6302)<5))=1;
    this_pixel_mask((this_wavelengths>6868) & (this_wavelengths<6932))=1;
    this_pixel_mask((this_wavelengths>7594) & (this_wavelengths<7700))=1;
  end

  % masking pixels with SN<4
  this_pixel_mask(this_flux./sqrt(this_noise_variance)<4)=1;
  
  this_rest_wavelengths = emitted_wavelengths(this_wavelengths, all_zqso_dr7(i));
  % normalizing here
 
  ind = (this_rest_wavelengths >= normalization_min_lambda) & ...
          (this_rest_wavelengths <= normalization_max_lambda) & ...
          (~this_pixel_mask);
 
  this_median = nanmedian(this_flux(ind));
  


  % bit 2: cannot normalize (all normalizing pixels are masked)
  if (isnan(this_median))
    filter_flags(i) = bitset(filter_flags(i), 3, true);
    continue;
  end
  
  ind = (this_rest_wavelengths >= min_lambda) & ...
        (this_rest_wavelengths <= max_lambda) & ...
        (~this_pixel_mask);
        
  % bit 3: not enough pixels available
  if (nnz(ind) < min_num_pixels)
    filter_flags(i) = bitset(filter_flags(i), 4, true);
    continue;
  end

  all_normalizers(i) = this_median;
 
  this_flux           = this_flux           / this_median;
  this_noise_variance = this_noise_variance / this_median^2;
 

  % ind = (this_rest_wavelengths >= loading_min_lambda) & ...
  %       (this_rest_wavelengths <= loading_max_lambda);

  % add one pixel on either side
  available_ind = find(~ind & ~this_pixel_mask);
  ind(min(available_ind(available_ind > find(ind, 1, 'last' )))) = true;
  ind(max(available_ind(available_ind < find(ind, 1, 'first')))) = true;

  all_wavelengths{i}    =    this_wavelengths(ind);%no longer (ind)
  all_flux{i}           =           this_flux(ind);
  all_noise_variance{i} = this_noise_variance(ind);
  all_pixel_mask{i}     =     this_pixel_mask(ind);
  all_sigma_pixel{i}    =     this_sigma_pixel(ind); 

  fprintf('loaded quasar %i of %i (%i/%i/%04i)\n', ...
          i, num_quasars, all_plate_dr7(i), all_mjd_dr7(i), all_fiber_dr7(i));
  
end

variables_to_save = {'normalization_min_lambda', 'normalization_max_lambda', ...
                     'min_num_pixels', 'all_wavelengths', 'all_flux', ...
                     'all_noise_variance', 'all_pixel_mask', ...
                     'all_normalizers', 'all_sigma_pixel', 'filter_flags'};
save(sprintf('%s/preloaded_qsos-mnp-%d', processed_directory(release),...
               min_num_pixels), ...
     variables_to_save{:});

% % write new filter flags to catalog
% save(sprintf('%s/filter_flags', processed_directory(release)), ...
%      'filter_flags');

toc;

