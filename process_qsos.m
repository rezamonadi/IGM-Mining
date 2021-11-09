% process_qsos: run DLA detection algorithm on specified objects
% load C4 catalog
Full_catalog = ...
    load(sprintf('%s/catalog', processed_directory(training_release)));
% train_ind -> those LOSs without any problem (filter_flag==0) that does not have
% Civ and useful for training null model (a model without Civ absorption line)
% prior_ind -> those LOSs with Civ absorption and in the half part of test
% test_ind -> second half without any filter flag for testing null and absorption model
% on the LOSs that we know have Civ or not. So, we asses our algorithm in this way.

if (ischar(prior_ind))
    prior_ind = eval(prior_ind);
end

% My prior_ind here is already those OK sight of lines that have CIV
prior.z_qsos  = Full_catalog.all_zqso(prior_ind);
prior.c4_ind = prior_ind;
prior.z_c4 = Full_catalog.all_z_civ(prior_ind,:);

% filter out CIVs from prior catalog corresponding to region of spectrum below
% Ly-alpha QSO rest. In Roman's code, for detecting DLAs, instead of
% Ly-alpha they have Lyman-limit
for i = size(prior.z_c4)
    for j=1:17
        if (observed_wavelengths(civ_1548_wavelength, prior.z_c4(i,j)) < ...
                observed_wavelengths(min_lambda, prior.z_qsos(i)))
            % make a list of all z_c4s and make a list of their z_zqsos
            % prior.z_qsos has the same length as all of the c4 systems
            prior.c4_ind(i) = false;
        end
    end
end
prior = rmfield(prior, 'z_c4');

% enable processing specific QSOs via setting to_test_ind
if (ischar(test_ind))
    test_ind = eval(test_ind);
end

all_wavelengths    =    all_wavelengths(test_ind);
all_flux           =           all_flux(test_ind);
all_noise_variance = all_noise_variance(test_ind);
all_pixel_mask     =     all_pixel_mask(test_ind);
all_sigma_pixel    =    all_sigma_pixel(test_ind);
z_qsos             =   Full_catalog.all_zqso(test_ind);
num_quasars        =                numel(z_qsos);
% preprocess model interpolants
% griddedInterpolant does an interpolation and gives a function handle
% based on the grided data. If the data is 2xD, {x,y} are the the same size as
% row  and columns. M is like M=f(x,y) and is like a matrix with each element
% M(i,j) = f(x(i), y(j))
mu_interpolator = ...
    griddedInterpolant(rest_wavelengths,        mu,        'linear');
M_interpolator = ...
    griddedInterpolant({rest_wavelengths, 1:k}, M,         'linear');


sample_sigma_c4 = min_sigma + (max_sigma-min_sigma)*offset_sigma_samples;

% initialize results with nan
min_z_c4s                   = nan(num_quasars, 1);
max_z_c4s                   = nan(num_quasars, 1);
log_priors_no_c4            = nan(num_quasars, 1);
log_priors_c4               = nan(num_quasars, 1);
log_likelihoods_no_c4       = nan(num_quasars, 1);
sample_log_likelihoods_c4   = nan(num_quasars, num_C4_samples);
log_likelihoods_c4          = nan(num_quasars, 1);
log_posteriors_no_c4        = nan(num_quasars, 1);
log_posteriors_c4           = nan(num_quasars, 1);
map_N_c4                  = nan(num_quasars, 1);
map_z_c4                 = nan(num_quasars, 1);
map_sigma_c4              = nan(num_quasars, 1);
Dz                          = nan(num_quasars, 1);
all_residual0               = nan(num_quasars, 1);
all_residual2               = nan(num_quasars, 1);
all_offsetl2                = nan(num_quasars, 1);
DoubletRatio                = nan(num_quasars, 1);
EqW1                        = nan(num_quasars, 1);
EqW2                        = nan(num_quasars, 1);
DLikelihood                 = nan(num_quasars, 1);
fw1_100                     = nan(num_quasars, 1);
fw1_200                     = nan(num_quasars, 1);
fw2_100                     = nan(num_quasars, 1);
fw2_200                     = nan(num_quasars, 1);
min_ratio                   = nan(num_quasars, 1);
N_pts1                      = nan(num_quasars, 1);
N_pts2                      = nan(num_quasars, 1);


for quasar_ind = 1:num_quasars
    tic;
    z_qso = z_qsos(quasar_ind);
    fprintf('processing quasar %i/%i (z_QSO = %0.4f) ...', ...
                              quasar_ind, num_quasars, z_qso);
    
    this_wavelengths    =    all_wavelengths{quasar_ind};
    this_wavelengths    =              this_wavelengths';
    this_flux           =           all_flux{quasar_ind}; 
    this_flux           =                     this_flux';
    this_noise_variance = all_noise_variance{quasar_ind};
    this_noise_variance =           this_noise_variance';
    this_pixel_mask     =     all_pixel_mask{quasar_ind};
    this_pixel_mask     =               this_pixel_mask';
    this_sigma_pixel    =     all_sigma_pixel{quasar_ind};
    this_sigma_pixel    =               this_sigma_pixel';
        
    %-------------interpolating observed flux ------------------\
    flux_interpolant = griddedInterpolant(this_wavelengths, this_flux);
    noise_variance_interpolant = griddedInterpolant(this_wavelengths, this_noise_variance);

    % building a finer wavelength array by adding the mean of ith and ith +1 element
    fine_wavelength = zeros(2*length(this_wavelengths)-1,1);
    fine_pixel_mask = zeros(2*length(this_wavelengths)-1,1);
    for i=1:length(this_wavelengths)-1
        fine_wavelength(2*i) = (this_wavelengths(i) + this_wavelengths(i+1))/2;
        fine_wavelength(2*i-1) = this_wavelengths(i);
        fine_pixel_mask(2*i) = (this_pixel_mask(i) + this_pixel_mask(i+1))/2;
        fine_pixel_mask(2*i-1) = this_pixel_mask(i);
        

    end
    fine_wavelength(2*length(this_wavelengths)-1) = this_wavelengths(end);
    fine_pixel_mask(2*length(this_wavelengths)-1) = this_pixel_mask(end);
    fine_pixel_mask(fine_pixel_mask==0.5)=1;
    
    this_wavelengths= fine_wavelength;
    this_flux = flux_interpolant(this_wavelengths);
    this_noise_variance = noise_variance_interpolant(this_wavelengths);
    % convert to QSO rest frame
    this_rest_wavelengths = emitted_wavelengths(this_wavelengths, z_qso);
    
    unmasked_ind = (this_rest_wavelengths >= min_lambda) & ...
        (this_rest_wavelengths <= max_lambda); %& (this_sigma_pixel>0);
    % keep complete copy of equally spaced wavelengths for absorption
    % computation
    this_unmasked_wavelengths = this_wavelengths(unmasked_ind);

    % interpolate 
    
    % [mask_ind] remove flux pixels with pixel_mask; pixel_mask is defined
    % in read_spec_dr7.m
    ind = unmasked_ind & (~fine_pixel_mask);% & (this_sigma_pixel>0);
    
    this_wavelengths      =      this_wavelengths(ind);
    this_rest_wavelengths = this_rest_wavelengths(ind);
    this_flux             =             this_flux(ind);
    this_noise_variance = this_noise_variance(ind);
    % this_noise_variance   =   this_noise_variance(ind);
    % this_sigma_pixel      =      this_sigma_pixel(ind);
    %   this_lya_zs = ...
    %       (this_wavelengths - lya_wavelength) / ...
    %       lya_wavelength;
    
    % c4 existence prior
    less_ind = (prior.z_qsos < (z_qso + prior_z_qso_increase));
    this_num_c4    = nnz(prior.c4_ind(less_ind));
    this_num_quasars = nnz(less_ind);
    this_p_c4 = this_num_c4 / this_num_quasars;
    log_priors_c4(quasar_ind) = ...
        log(                   this_num_c4) - log(this_num_quasars);
    log_priors_no_c4(quasar_ind) = ...
        log(this_num_quasars - this_num_c4) - log(this_num_quasars);
    
    % ind_dla = ismember(dla_QSO_ID, all_QSO_ID(quasar_ind));
    % if nnz(ind_dla)>0
                           
    %     log_posteior_dla = log_posteriors_dla(ind_dla);
                              
    %     log_posteior_no_dla = log_posteriors_no_dla(ind_dla);
    %     log_priors_c4(quasar_ind) = log_priors_c4(quasar_ind)+ ...
    %                                 log(b)+log_posteior_dla;
    %     log_priors_no_c4(quasar_ind) = log_priors_no_c4(quasar_ind)+ ...
    %                                 log(b)+log_posteior_no_dla;                                    
    % end

    fprintf('\n');
    fprintf(' ...     p(   CIV | z_QSO)  mvn      : %0.3f\n',     this_p_c4);
    fprintf(' ...     p(no CIV | z_QSO)        : %0.3f\n', 1 - this_p_c4);
    
    % interpolate model onto given wavelengths
    this_mu = mu_interpolator( this_rest_wavelengths);
    this_M  =  M_interpolator({this_rest_wavelengths, 1:k});
    
    %   this_log_omega = log_omega_interpolator(this_rest_wavelengths);
    %   this_omega2 = exp(2 * this_log_omega);
    
    %   this_scaling_factor = 1 - exp(-tau_0 .* (1 + this_lya_zs).^beta) + c_0;
    
    %   this_omega2 = this_omega2 .* thiscaling_factor.^2;
    
    % baseline: probability of no DLA model
    %   disp(size(this_M))
    %   disp(size(this_flux))
    %   disp(size(this_noise_variance))
    log_likelihoods_no_c4(quasar_ind) = ...
        log_mvnpdf_low_rank(this_flux, this_mu, this_M, ...
        this_noise_variance);
    
    log_posteriors_no_c4(quasar_ind) = ...
        log_priors_no_c4(quasar_ind) + log_likelihoods_no_c4(quasar_ind);
    
    fprintf(' ... log p(D | z_QSO, no CIV)     : %0.2f\n', ...
        log_likelihoods_no_c4(quasar_ind));
    
    min_z_c4s(quasar_ind) = min_z_c4(this_wavelengths, z_qso);
    max_z_c4s(quasar_ind) = max_z_c4(z_qso, max_z_cut);
    
    sample_z_c4 = ...
        min_z_c4s(quasar_ind) +  ...
        (max_z_c4s(quasar_ind) - min_z_c4s(quasar_ind)) * offset_z_samples;


    % Temperature samples
    
    % ensure enough pixels are on either side for convolving with
    % instrument profile
    padded_wavelengths = ...
        [logspace(log10(min(this_unmasked_wavelengths)) - width * pixel_spacing, ...
        log10(min(this_unmasked_wavelengths)) - pixel_spacing,         ...
        width)';                                                       ...
        this_unmasked_wavelengths;                                              ...
        logspace(log10(max(this_unmasked_wavelengths)) + pixel_spacing,         ...
        log10(max(this_unmasked_wavelengths)) + width * pixel_spacing, ...
        width)'                                                        ...
        ];
    
    % [mask_ind] to retain only unmasked pixels from computed absorption profile
    % this has to be done by using the unmasked_ind which has not yet
    % been applied this_pixel_mask.
    ind = (~fine_pixel_mask(unmasked_ind));
    
    % compute probabilities under DLA model for each of the sampled
    % (normalized offset, log(N HI)) pairs
    parfor i = 1:num_C4_samples
  
        
        % Limitting red-shift in the samples
        % this_z_c4 = (this_wavelengths / 1549) - 1;
        % absorption corresponding to this sample with one absorption line as a noise model 

        %     dla_omega2 = this_omega2 .* absorption.^2;
        
        % absorption corresponding to this sample with two absorption lines as a doublet 
        num_lines=2;
        % absorption = voigt(padded_wavelengths, sample_z_c4(i), ...
        % nciv_samples(i),num_lines , sample_sigma_c4(i), this_sigma_pixel);
        absorption = voigt0(padded_wavelengths, sample_z_c4(i), ...
        nciv_samples(i),num_lines , sample_sigma_c4(i));

        % precomputing voigt function and putting the average
        %  of the voigt instead of the actual value

        wavelengths_fine = linspace(min(this_unmasked_wavelengths),...
        max(this_unmasked_wavelengths), 100*length(this_unmasked_wavelengths));
        wavelengths_fine = wavelengths_fine';
        dl = wavelengths_fine(2) - wavelengths_fine(1);
        
        % voigt_fine= ones(length(wavelengths_fine), 1);
        % mask_not_ones_voigt = (wavelengths_fine>=1538*(1+sample_z_c4(i))) ...
        %                 & (wavelengths_fine<=1560*(1+sample_z_c4(i)));
        
        % padded_wavelengths_fine = [...
        % linspace(min(wavelengths_fine(mask_not_ones_voigt)) - width * dl, ...
        % min(wavelengths_fine(mask_not_ones_voigt)) - dl, width)';                                                       ...
        % wavelengths_fine(mask_not_ones_voigt);                                              ...
        % linspace(min(wavelengths_fine(mask_not_ones_voigt)) + dl, ...
        % min(wavelengths_fine(mask_not_ones_voigt)) + width*dl, width)'];
        % voigt_fine(mask_not_ones_voigt) = ...
        %             voigt0(padded_wavelengths_fine, sample_z_c4(i), ...
        %     nciv_samples(i),num_lines , sample_sigma_c4(i));
        % plot(this_unmasked_wavelengths, absorption, 'c', 'Marker', 'o')
        % hold on 
        % plot(wavelengths_fine, voigt_fine, 'r'); hold on 
        % plot(wavelengths_fine(mask_not_ones_voigt), voigt_fine(mask_not_ones_voigt), 'b', 'Marker', 'o')
        % hold off
        % for j=1:length(this_unmasked_wavelengths)
        %     mask_j = abs(this_unmasked_wavelengths(j) - wavelengths_fine)<5*dl;
        %     % sum(mask_j)
        %     absorption(j)= mean(voigt_fine(mask_j));
        %     % fprintf(' mean_voigt:%d\n', mean(voigt_fine(mask_j)));
        end
      
        absorption = absorption(ind);
        c4_mu     = this_mu     .* absorption;
        c4_M      = this_M      .* absorption;
        %     dla_omega2 = this_omega2 .* absorption.^2;
        sample_log_likelihoods_c4(quasar_ind, i) = ...
          log_mvnpdf_low_rank(this_flux, c4_mu, c4_M, ...
          this_noise_variance);
        
    end
       
    % compute sample probabilities and log likelihood of DLA model in
    % numerically safe manner for one line
    
       
    % compute sample probabilities and log likelihood of DLA model in
    % numerically safe manner for  doublet 
    max_log_likelihood = max(sample_log_likelihoods_c4(quasar_ind, :));
    sample_probabilities = ...
        exp(sample_log_likelihoods_c4(quasar_ind, :) - ...
        max_log_likelihood);
    log_likelihoods_c4(quasar_ind) = ...
        max_log_likelihood + log(mean(sample_probabilities));
    
    log_posteriors_c4(quasar_ind) = ...
        log_priors_c4(quasar_ind) + log_likelihoods_c4(quasar_ind);
    
    fprintf(' ... log p(D | z_QSO,    CIV)     : %0.2f\n', ...
        log_likelihoods_c4(quasar_ind));
    fprintf(' ... log p(CIV | D, z_QSO)        : %0.2f\n', ...
        log_posteriors_c4(quasar_ind));
     
    
    

    [max2, maxind] = nanmax(sample_log_likelihoods_c4(quasar_ind, :));
    [min2, ~] = nanmin(sample_log_likelihoods_c4(quasar_ind, :));
    DLikelihood(quasar_ind,1) = max2 -min2;

    map_z_c4(quasar_ind)    = sample_z_c4(maxind);        
    Dz(quasar_ind) = (z_qso - map_z_c4(quasar_ind))/(1+z_qso);

    
        
    

    map_N_c4(quasar_ind)  = log_nciv_samples(maxind);
    map_sigma_c4(quasar_ind)  = sample_sigma_c4(maxind);
    fprintf('\nmap(N): %.2f, map(z_c4): %.2f, map(b/1e5): %.2f\n',map_N_c4(quasar_ind),...
    map_z_c4(quasar_ind), map_sigma_c4(quasar_ind)/1e5);
    fprintf(' took %0.3fs.\n', toc);


    % % Goodness of fit evaluation
    % all_residual0(quasar_ind,1) = mean(abs(this_flux - this_mu));
    
    
    % % absorptionL2 = voigt(padded_wavelengths, map_z_c4L2(quasar_ind), ...
    % %                      10^map_N_c4L2(quasar_ind), num_lines,...
    % %                      map_sigma_c4L2(quasar_ind), this_sigma_pixel);
    % % c4_muL2     = this_mu     .* absorptionL2;
    % % all_residual2(quasar_ind,1) = mean(abs(this_flux - c4_muL2));
    % % all_offsetl2(quasar_ind,1)  = mean(this_flux - c4_muL2);
    % % fprintf('Residual-> 0: %.5f, 1: %.5f, 2: %.5f\n', all_residual0(quasar_ind),...
    %                     %    all_residual1(quasar_ind), all_residual2(quasar_ind));
    % % dw = 0.01;
    % % w1 = 1548-5:dw:1548+5;
    % % w2 = 1550-5:dw:1550+5;
    % % padded_wavelengths1 = [linspace(1548-5, 1548-5+width*dw, width)';
    % %                       w1'; linspace(1548+5, 1548+width*dw, width)']; 
    % % padded_wavelengths2 = [linspace(1550-5, 1550-5+width*dw, width)';
    % %                       w2'; linspace(1550+5, 1550+width*dw, width)'];
    % L1 = voigt(padded_wavelengths,...
    %             map_z_c4(quasar_ind), 10^map_N_c4(quasar_ind), 1,...
    %                       map_sigma_c4(quasar_ind), this_sigma_pixel);
    % L2 = voigt_dr(padded_wavelengths,...
    %             map_z_c4(quasar_ind), 10^map_N_c4(quasar_ind), 1,...
    %                       map_sigma_c4(quasar_ind), this_sigma_pixel);
 
    % EqW1(quasar_ind,1) =trapz(this_unmasked_wavelengths, 1-L1);
    % EqW2(quasar_ind,1) =trapz(this_unmasked_wavelengths, 1-L2);
    
    % % min ratio 
    % ind_min1 = find(L1 == min(L1));
    % ind_min2 = find(L2 == min(L2));
    % min_ratio(quasar_ind,1) = (max(L1) - L1(ind_min1))/(max(L2) - L2(ind_min2));

    
    
    % % fw_100 width at 0.01 belwo the max(L)
    % H1 = max(L1) - min(L1);
    % H2 = max(L2) - min(L2);

    % ind_sides1 = find(max(L1) - L1 >=0.01*H1);
    % ind_sides2 = find(max(L2) - L2 >=0.01*H2);
    % Dips1 = sort(L1(ind_sides1), 'descend');
    % Dips2 = sort(L2(ind_sides2), 'descend');
    % fw1_100(quasar_ind,1) = abs(this_unmasked_wavelengths(L1==Dips1(1))-...
    %                           this_unmasked_wavelengths(L1==Dips1(2)));
    % fw2_100(quasar_ind,1) = abs(this_unmasked_wavelengths(L2==Dips2(1))-...
    %                           this_unmasked_wavelengths(L2==Dips2(2)));
    % % fw_200 width at 0.02 belwo the ma(L)
    % ind_sides1 = find(max(L1) - L1 >=0.02*H1);
    % ind_sides2 = find(max(L2) - L2 >=0.02*H2);
    % Dips1 = sort(L1(ind_sides1), 'descend');
    % Dips2 = sort(L2(ind_sides2), 'descend');
    % fw1_200(quasar_ind,1) = abs(this_unmasked_wavelengths(L1==Dips1(1))-...
    %                           this_unmasked_wavelengths(L1==Dips1(2)));
    % fw2_200(quasar_ind,1) = abs(this_unmasked_wavelengths(L2==Dips2(1))-...
    %                           this_unmasked_wavelengths(L2==Dips2(2)));                              


    % % plot(this_unmasked_wavelengths,L1, '-'); 
    % % hold on; 
    % % mask  = (L1 == Dips1(1)) | (L1 == Dips1(2));
    % % plot(this_unmasked_wavelengths(mask), L1(mask), 'X');
    % % hold on 
    % % plot(this_unmasked_wavelengths,L2, '-'); 
    % % hold on; 
    % % mask  = (L2 == Dips2(1)) | (L2 == Dips2(2));
    % % plot(this_unmasked_wavelengths(mask), L2(mask), 'X');
    % % break
    % DoubletRatio(quasar_ind,1)=EqW1(quasar_ind,1)/EqW2(quasar_ind,1);
    % fprintf('EW1: %.2e, EW2:%.2e, DR: %.2e, DL:%.2e\nfw1_100:%.3f, fw2_100:%.3f\n fw1_200:%.3f fw2_200:%.3f',...
    %  EqW1(quasar_ind,1),EqW2(quasar_ind,1), DoubletRatio(quasar_ind,1),...
    %  DLikelihood(quasar_ind,1), fw1_100(quasar_ind,1), fw2_100(quasar_ind,1),...
    % fw1_200(quasar_ind,1), fw2_200(quasar_ind,1));



    % % number of points around absorption
    % l1 = 1550.7810;
    % l2 = 1548.2040;
    % mask_pts1 = (this_unmasked_wavelengths>(1+map_z_c4(quasar_ind))*(l1-0.5)) & ...
    %             (this_unmasked_wavelengths<(1+map_z_c4(quasar_ind))*(l1+0.5));
    % mask_pts2 = (this_unmasked_wavelengths>(1+map_z_c4(quasar_ind))*(l2-0.5)) & ...
    %             (this_unmasked_wavelengths<(1+map_z_c4(quasar_ind))*(l2+0.5));

    % N_pts1(quasar_ind,1) = sum(mask_pts1);
    % N_pts2(quasar_ind,1) = sum(mask_pts2);
   
end



% compute model posteriors in numerically safe manner
max_log_posteriors = ...
    max([log_posteriors_no_c4, log_posteriors_c4], [], 2);

model_posteriors = ...
    exp([log_posteriors_no_c4, log_posteriors_c4] - max_log_posteriors);

model_posteriors = model_posteriors ./ sum(model_posteriors, 2);

p_no_c4 = model_posteriors(:, 1);
p_c4    = 1 - p_no_c4;

% save results
variables_to_save = {'training_release', 'training_set_name', ...
    'c4_catalog_name', 'prior_ind', 'release', ...
    'test_ind', 'prior_z_qso_increase', ...
    'max_z_cut', 'min_z_c4s', 'max_z_c4s', ...
    'log_priors_no_c4', 'log_priors_c4', ...
    'log_likelihoods_no_c4',  ...
    'sample_log_likelihoods_c4', 'log_likelihoods_c4'...
    'log_posteriors_no_c4',  'log_posteriors_c4',...
    'model_posteriors', 'p_no_c4', 'p_c4' ...
    'map_z_c4', 'map_N_c4', 'map_sigma_c4', 'Dz',...
     'all_residual0',  'all_residual2', ...
     'DoubletRatio', 'all_offsetl2', 'EqW1', 'EqW2', 'DLikelihood','fw1_100','fw2_100',...
     'fw1_200', 'fw2_200', 'min_ratio', 'N_pts1', 'N_pts2'};

filename = sprintf('%s/processed_qsos_R%s', ...
    processed_directory(release), ...
    testing_set_name);

save(filename, variables_to_save{:}, '-v7.3');
