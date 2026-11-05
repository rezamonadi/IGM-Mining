% preload_qsos: loads spectra from SDSS FITS files, applies further
% filters, and applies some basic preprocessing such as normalization
% and truncation to the region of interest

% load QSO catalog
tic;

% ----------dr7----------------------------------------
variables_to_load = {'all_plate_dr7', 'all_mjd_dr7', 'all_fiber_dr7',...
'all_QSO_ID','all_zqso'};



load('data/dr7/processed/catalog.mat',  variables_to_load{:});

num_quasars = numel(all_zqso);

all_wavelengths    =  cell(num_quasars, 1);
all_flux           =  cell(num_quasars, 1);
all_noise_variance =  cell(num_quasars, 1);
all_pixel_mask     =  cell(num_quasars, 1);
all_sigma_pixel     =  cell(num_quasars, 1);
all_normalizers    = zeros(num_quasars, 1);

% Diagnostics for the sky-line masking and subsequent quality cuts.
num_sky_5579_masked = 0;
num_sky_6302_masked = 0;


% to track reasons for filtering out QSOs
filter_flags = zeros(num_quasars, 1, 'uint8');
% 
% filtering bit 0: z_QSO < 1.5
ind = (all_zqso < z_qso_cut);
filter_flags(ind) = bitset(filter_flags(ind), 1, true);

% This is the population reaching the normalization and usable-pixel tests.
num_before_quality_cuts = nnz(filter_flags == 0);



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
  
  % % % Masking Sky lines  (Cooksey+2013)

  sky_5579_ind = abs(this_wavelengths - 5579) < 5;
  sky_6302_ind = abs(this_wavelengths - 6302) < 5;
  num_sky_5579_masked = num_sky_5579_masked + nnz(sky_5579_ind & ~this_pixel_mask);
  num_sky_6302_masked = num_sky_6302_masked + nnz(sky_6302_ind & ~this_pixel_mask);
  this_pixel_mask(sky_5579_ind | sky_6302_ind) = 1;
  this_rest_wavelengths = emitted_wavelengths(this_wavelengths, all_zqso(i));
  % normalizing here
 
  ind = (this_rest_wavelengths >= normalization_min_lambda) & ...
          (this_rest_wavelengths <= normalization_max_lambda) & ...
          (~this_pixel_mask);
 
  this_median = nanmedian(this_flux(ind));
  
  % bit 2: cannot normalize (all normalizing pixels are masked)
  if (~isfinite(this_median) || this_median <= 0)
    filter_flags(i) = bitset(filter_flags(i), 3, true);
    continue;
  end
  
  ind = (this_rest_wavelengths >= min_lambda) & ...
        (this_rest_wavelengths <= max_lambda) & ...
        (~this_pixel_mask) & (this_sigma_pixel>0);
        
  % bit 3: not enough pixels available
  if (nnz(ind) < min_num_pixels)
    filter_flags(i) = bitset(filter_flags(i), 4, true);
    continue;
  end

  all_normalizers(i) = this_median;

  this_flux           = this_flux           / this_median;
  this_noise_variance = this_noise_variance / this_median^2;
 
  % add one pixel on either side
  available_ind = find(~ind & ~this_pixel_mask);
  first_ind = find(ind, 1, 'first');
  last_ind = find(ind, 1, 'last');
  right_ind = available_ind(available_ind > last_ind);
  left_ind = available_ind(available_ind < first_ind);
  if ~isempty(right_ind)
      ind(min(right_ind)) = true;
  end
  if ~isempty(left_ind)
      ind(max(left_ind)) = true;
  end

  all_wavelengths{i}    =    this_wavelengths(ind);%no longer (ind)
  all_flux{i}           =           this_flux(ind);
  all_noise_variance{i} = this_noise_variance(ind);
  all_pixel_mask{i}     =     this_pixel_mask(ind);
  all_sigma_pixel{i}    =     this_sigma_pixel(ind); 

  fprintf('loaded quasar %i of %i (%i/%i/%04i)\n', ...
          i, num_quasars, all_plate_dr7(i), all_mjd_dr7(i), all_fiber_dr7(i));
  
end

variables_to_save = {'loading_min_lambda', 'loading_max_lambda', ...
                     'normalization_min_lambda', 'normalization_max_lambda', ...
                     'min_num_pixels', 'all_wavelengths', 'all_flux', ...
                     'all_noise_variance', 'all_pixel_mask', ...
                     'all_normalizers', 'all_sigma_pixel', 'filter_flags', ...
                     'num_sky_5579_masked', 'num_sky_6302_masked'};
save('data/dr7/processed/preloaded_qsos_dr7.mat', variables_to_save{:}, '-v7.3');

% write new filter flags to catalog
save(sprintf('%s/filter_flags', processed_directory(release)), ...
     'filter_flags');

num_retained = nnz(cellfun(@(x) ~isempty(x), all_wavelengths));
fprintf('sky-line pixels newly masked: 5579 A window = %d, 6302 A window = %d\n', ...
        num_sky_5579_masked, num_sky_6302_masked);
fprintf('spectra before normalization/pixel quality cuts: %d; retained: %d\n', ...
        num_before_quality_cuts, num_retained);

toc;


