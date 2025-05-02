 % preload_qsos: loads spectra from SDSS FITS files, applies further
% filters, and applies some basic preprocessing such as normalization
% and truncation to the region of interest

% load QSO catalog
tic;

% ----------dr7----------------------------------------
variables_to_load = {'all_plate_dr7', 'all_mjd_dr7', 'all_fiber_dr7',...
'all_QSO_ID','all_zqso'};



load(sprintf('%s/catalog', processed_directory(release)), ...
    variables_to_load{:});

num_quasars = numel(all_zqso);

all_wavelengths    =  cell(num_quasars, 1);
all_flux           =  cell(num_quasars, 1);
all_noise_variance =  cell(num_quasars, 1);
all_pixel_mask     =  cell(num_quasars, 1);
all_sigma_pixel    =  cell(num_quasars, 1);
all_normalizers    = zeros(num_quasars, 1);


% to track reasons for filtering out QSOs
filter_flags = zeros(num_quasars, 1, 'uint8');
% 
% filtering bit 0: z_QSO < 0.36
ind = (all_zqso < z_qso_cut);
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
  
  % % % Masking Sky lines  (Cooksey+2013)

  this_pixel_mask((abs(this_wavelengths-5579)<5) | (abs(this_wavelengths-6302)<5))=1;
 
  this_rest_wavelengths = emitted_wavelengths(this_wavelengths, all_zqso(i));
  % normalizing here

  % if (all_zqso(i) < 0.6)
  % 
  %  ind = (this_rest_wavelengths >= normalization_min_lambda_1) & ...
  %          (this_rest_wavelengths <= normalization_max_lambda_1) & ...
  %          (~this_pixel_mask);
  % end
  % 
  % if (all_zqso(i)>= 0.6 && all_zqso(i) < 1.0)
  % 
  %  ind = (this_rest_wavelengths >= normalization_min_lambda_2) & ...
  %          (this_rest_wavelengths <= normalization_max_lambda_2) & ...
  %          (~this_pixel_mask);
  % end
  
  % if (all_zqso(i)>= 1.0 && all_zqso(i)< 2.5)
 
   ind = (this_rest_wavelengths >= normalization_min_lambda_3) & ...
           (this_rest_wavelengths <= normalization_max_lambda_3) & ...
           (~this_pixel_mask);
  % end
  
  % % if (all_zqso(i)>= 2.5 && all_zqso(i) > 4.7)
  % if (all_zqso(i)>= 2.5)
  % 
  %  ind = (this_rest_wavelengths >= normalization_min_lambda_4) & ...
  %          (this_rest_wavelengths <= normalization_max_lambda_4) & ...
  %          (~this_pixel_mask);

  % end

  this_median = nanmedian(this_flux(ind));
  
  % bit 2: cannot normalize (all normalizing pixels are masked)
  if (isnan(this_median))
    filter_flags(i) = bitset(filter_flags(i), 3, true);
    continue;
  end

  
  % bit 3: not enough pixels available

  ind = (this_rest_wavelengths >= min_lambda) & ...
        (this_rest_wavelengths <= max_lambda) & ...
        (this_sigma_pixel>0) & ...
        (~this_pixel_mask);
  
  if (nnz(ind) < min_num_pixels)
    filter_flags(i) = bitset(filter_flags(i), 4, true);
    continue;
  end

  all_normalizers(i) = this_median;

  this_flux           = this_flux           / this_median;
  this_noise_variance = this_noise_variance / this_median^2;
 
  % add one pixel on either side
  available_ind = find(~ind & ~this_pixel_mask);
  ind(min(available_ind(available_ind > find(ind, 1, 'last' )))) = true;
  ind(max(available_ind(available_ind < find(ind, 1, 'first')))) = true;

  all_wavelengths{i}    =    this_wavelengths(ind); %no longer (ind)
  all_flux{i}           =           this_flux(ind);
  all_noise_variance{i} = this_noise_variance(ind);
  all_pixel_mask{i}     =     this_pixel_mask(ind);
  all_sigma_pixel{i}    =     this_sigma_pixel(ind); 

  fprintf('loaded quasar %i of %i (%i/%i/%04i) %i\n', ...
          i, num_quasars, all_plate_dr7(i), all_mjd_dr7(i), all_fiber_dr7(i), nnz(ind));
  
end

variables_to_save = {'loading_min_lambda',...
                     'loading_max_lambda', ...
                     'min_num_pixels',...
                     'all_wavelengths',...
                     'all_flux', ...
                     'all_noise_variance',...
                     'all_pixel_mask', ...
                     'all_normalizers',...
                     'all_sigma_pixel',...
                     'filter_flags'};
save(sprintf('%s/preloaded_qsos_%s.mat', processed_directory(release), training_set_name), ...
     variables_to_save{:});

% write new filter flags to catalog
save(sprintf('%s/filter_flags', processed_directory(release)), ...
     'filter_flags');

toc;

