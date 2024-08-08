% process_qsos: run CIV detection algorithm on specified objects

% Averaging     -> yes
% Single model  -> yes
% Multi CIV     -> yes
% Multi Singlet -> no
%
% removing CIV found region after each run 
% also removing the map Z_civ
% load C4 catalog
Full_catalog = ...
    load(sprintf('%s/catalog', processed_directory(release)));
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
prior.z_c4 = Full_catalog.all_z_civ(prior_ind);

% filter out CIVs from prior catalog corresponding to region of spectrum below
% Ly-alpha QSO rest. In Roman's code, for detecting DLAs, instead of
% Ly-alpha they have Lyman-limit
for i = size(prior.z_c4)
    if (observed_wavelengths(civ_1548_wavelength , prior.z_c4(i)) < ...
            observed_wavelengths(min_lambda, prior.z_qsos(i)))
        prior.c4_ind(i) = false;
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
z_qsos             =           all_zqso(test_ind);
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

% initialize results with nan
min_z_c4s                   = nan(num_quasars, 1);
max_z_c4s                   = nan(num_quasars, 1);
log_priors_no_c4            = nan(num_quasars, 1);
log_priors_c4               = nan(num_quasars, max_civ);
log_likelihoods_no_c4       = nan(num_quasars, max_civ);
sample_log_likelihoods_c4L1 = nan(num_quasars, num_C4_samples);
sample_log_likelihoods_c4L2 = nan(num_quasars, num_C4_samples, max_civ);
log_likelihoods_c4L1        = nan(num_quasars, max_civ);
log_likelihoods_c4L2        = nan(num_quasars, max_civ);
log_posteriors_no_c4        = nan(num_quasars, max_civ);
log_posteriors_c4L1         = nan(num_quasars, max_civ);
log_posteriors_c4L2         = nan(num_quasars, max_civ);
map_N_c4L1                  = nan(num_quasars, max_civ);
map_N_c4L2                  = nan(num_quasars, max_civ);
map_z_c4L1                  = nan(num_quasars, max_civ);
map_z_c4L2                  = nan(num_quasars, max_civ);
map_sigma_c4L1              = nan(num_quasars, max_civ);
map_sigma_c4L2              = nan(num_quasars, max_civ);
p_c4                        = nan(num_quasars, max_civ);
p_c4L1                        = nan(num_quasars, max_civ);
p_no_c4                        = nan(num_quasars, max_civ);
EW                          = nan(num_quasars, max_civ);
num_pixel_civ               = nan(num_quasars, max_civ, 2);
sample_sigma_c4 = min_sigma + (max_sigma-min_sigma)*offset_sigma_samples;
% sample_sigma_c4 = sigma_samples;
ID = all_QSO_ID(test_ind);
% plt_count=0;
% FN_IDs = importdata('FN-list.csv');
% in_Kathy_FN_list = ismember(ID, FN_IDs);
N_civ_test = all_N_civ(test_ind,:);
z_PM_test = all_z_civ(test_ind,:);
z_PM_prior = all_z_civ(prior_ind,:);
% ind_E = 49 + ind_S;
% if (ind_S==1201)
    % ind_E = 1301;
% end
% for quasar_ind = ind_S:ind_E
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
    
    % convert to QSO rest frame
    this_rest_wavelengths = emitted_wavelengths(this_wavelengths, z_qso);
    
    unmasked_ind = (this_rest_wavelengths >= min_lambda) & ...
        (this_rest_wavelengths <= max_lambda);% & (this_sigma_pixel>0);
    % keep complete copy of equally spaced wavelengths for absorption
    % computation
    this_unmasked_wavelengths = this_wavelengths(unmasked_ind);
    
    % [mask_ind] remove flux pixels with pixel_mask; pixel_mask is defined
    % in read_spec_dr7.m
    ind = unmasked_ind & (~this_pixel_mask);
    
    this_wavelengths      =      this_wavelengths(ind);
    this_rest_wavelengths = this_rest_wavelengths(ind);
    this_flux             =             this_flux(ind);
    this_noise_variance   =   this_noise_variance(ind);
    this_sigma_pixel      =      this_sigma_pixel(ind);
    %   this_lya_zs = ...
    %       (this_wavelengths - lya_wavelength) / ...
    %       lya_wavelength;
    
    % c4 existence prior
    less_ind = (prior.z_qsos < (z_qso + prior_z_qso_increase));
    less_systems = z_PM_prior(less_ind,:);
    this_num_quasars = nnz(less_ind);
    for i=1:max_civ
        this_num_c4 = nnz(less_systems(:,i)>0 & less_systems(:,i+1)==-1);
        this_p_c4(i)   =   (this_num_c4 / this_num_quasars);
    end

    % the prior for having exactly k CIVs is (M / N)^k - (M / N)^(k + 1) 
    % for i =1:(max_civ-1)
        % this_p_c4(i) = this_p_c4(i) - this_p_c4(i+1);
    % end

    log_priors_c4(quasar_ind,:) = log(this_p_c4(:));
    log_priors_no_c4(quasar_ind) = ...
        log(this_num_quasars - this_num_c4) - log(this_num_quasars);
    
    fprintf('\n');
    for i = 1:max_civ
        fprintf(' ...     p(%i  CIVs | z_QSO)       : %0.3f\n', i, this_p_c4(i));
    end
    fprintf(' ...     p(no CIV  | z_QSO)       : %0.3f\n', exp(log_priors_no_c4(quasar_ind)) );
    
    % interpolate model onto given wavelengths
    this_mu = mu_interpolator( this_rest_wavelengths);
    this_M  =  M_interpolator({this_rest_wavelengths, 1:k});
    
    
    
    min_z_c4s(quasar_ind) = min_z_c4(this_wavelengths, z_qso);
    % instead of this_wavelengths I puting 1310A where is the lower limit of CIV search in C13
    min_z_c4s(quasar_ind) = min_z_c4(1310, z_qso);
    max_z_c4s(quasar_ind) = max_z_c4(z_qso, max_z_cut);
    
    sample_z_c4 = ...
        min_z_c4s(quasar_ind) +  ...
        (max_z_c4s(quasar_ind) - min_z_c4s(quasar_ind)) * offset_z_samples;


    % Temperature samples
    
    % ensure enough pixels are on either side for convolving with
    % instrument profile

   % building a finer wavelength and mask arrays 
   % by adding the mean of ith and ith +1 element
    
    
    padded_wavelengths = ...
        [logspace(log10(min(this_unmasked_wavelengths)) - width * pixel_spacing, ...
        log10(min(this_unmasked_wavelengths)) - pixel_spacing,...
        width)';...
        this_unmasked_wavelengths;...
        logspace(log10(max(this_unmasked_wavelengths)) + pixel_spacing,...
        log10(max(this_unmasked_wavelengths)) + width * pixel_spacing,...
        width)'...
        ];


        % % when broadening is off
        % padded_wavelengths = this_unmasked_wavelengths;
    
    % [mask_ind] to retain only unmasked pixels from computed absorption profile
    % this has to be done by using the unmasked_ind which has not yet
    % been applied this_pixel_mask.
    ind = (~this_pixel_mask(unmasked_ind));
    
    % compute probabilities under DLA model for each of the sampled
    % (normalized offset, log(N HI)) pairs
    lenW_unmasked = length(this_unmasked_wavelengths);
    ind_not_remove = true(size(this_flux));
    for num_c4=1:max_civ
        fprintf('num_civ:%d\n',num_c4);
        this_z_1548 = (this_wavelengths / 1548.2040) - 1;
        this_z_1550 = (this_wavelengths / 1550.7810) - 1;
        if(num_c4>1)
            if((p_c4(quasar_ind, num_c4-1)>p_c4L1(quasar_ind, num_c4-1)) & ...
               (p_c4(quasar_ind, num_c4-1)>p_no_c4(quasar_ind, num_c4-1)))
                ind_not_remove = ind_not_remove  & ...
                    (abs(this_z_1548 - map_z_c4L2(quasar_ind, num_c4-1))>kms_to_z(dv_mask)*(1+z_qso)) & ...
                    (abs(this_z_1550 - map_z_c4L2(quasar_ind, num_c4-1))>kms_to_z(dv_mask)*(1+z_qso));
            end

            if((p_c4L1(quasar_ind, num_c4-1)>p_c4(quasar_ind, num_c4-1)) & ...
                        (p_c4L1(quasar_ind, num_c4-1)>p_no_c4(quasar_ind, num_c4-1)))
                ind_not_remove = ind_not_remove  & ...
                (abs(this_z_1548 - map_z_c4L1(quasar_ind, num_c4-1))>kms_to_z(dv_mask)*(1+z_qso));
            end

            if((p_no_c4(quasar_ind, num_c4-1)>p_c4(quasar_ind, num_c4-1)) & ...
                (p_no_c4(quasar_ind, num_c4-1)>p_c4L1(quasar_ind, num_c4-1)))
                
                fprintf('No more than %d CIVs in this spectrum.', num_c4-1)
                break;
            end

        end

        log_likelihoods_no_c4(quasar_ind, num_c4) = ...
        log_mvnpdf_low_rank(this_flux(ind_not_remove), this_mu(ind_not_remove),...
        this_M(ind_not_remove), this_noise_variance(ind_not_remove));

        log_posteriors_no_c4(quasar_ind, num_c4) = ...
            log_priors_no_c4(quasar_ind) + log_likelihoods_no_c4(quasar_ind, num_c4);

        fprintf(' ... log p(D | z_QSO, no CIV)     : %0.2f\n', ...
        log_likelihoods_no_c4(quasar_ind));
        parfor i = 1:num_C4_samples
            % Limitting red-shift in the samples
            % absorption corresponding to this sample with one absorption line as a noise model 

            % absorption corresponding to this sample with two absorption lines as a doublet 
            % compute fine absorption 
            num_lines=2;
            absorptionL2_fine = voigt_iP(finer(padded_wavelengths, nAVG), sample_z_c4(i), ...
            nciv_samples(i),num_lines , sample_sigma_c4(i), finer(this_sigma_pixel, nAVG));
            
            % average fine absorption and shrink it to the size of original array
            % as large as the unmasked_wavelengths
            
            absorptionL2 = Averager(absorptionL2_fine, nAVG, lenW_unmasked);
            absorptionL2 = absorptionL2(ind);
            c4_muL2     = this_mu     .* absorptionL2;
            c4_ML2      = this_M      .* absorptionL2;
            %     dla_omega2 = this_omega2 .* absorption.^2;
            
            sample_log_likelihoods_c4L2(quasar_ind, i, num_c4) = ...
            log_mvnpdf_low_rank(this_flux(ind_not_remove),...
                                c4_muL2(ind_not_remove),...
                                c4_ML2(ind_not_remove), ...
                                this_noise_variance(ind_not_remove));
            if SingleLineModel==1
                num_lines=1;
                absorptionL1_fine = voigt_iP(finer(padded_wavelengths, nAVG), sample_z_c4(i), ...
                nciv_samples(i),num_lines , sample_sigma_c4(i), finer(this_sigma_pixel, nAVG));
                
                % average fine absorption and shrink it to the size of original array
                % as large as the unmasked_wavelengths
                
                absorptionL1 = Averager(absorptionL1_fine, nAVG, lenW_unmasked);
                absorptionL1 = absorptionL1(ind);
                c4_muL1     = this_mu     .* absorptionL1;
                c4_ML1      = this_M      .* absorptionL1;
                %     dla_omega2 = this_omega2 .* absorption.^2;
                sample_log_likelihoods_c4L1(quasar_ind, i) = ...
                log_mvnpdf_low_rank(this_flux(ind_not_remove),...
                c4_muL1(ind_not_remove), c4_ML1(ind_not_remove), ...
                this_noise_variance(ind_not_remove));
            
            end
        end
            
        % compute sample probabilities and log likelihood of DLA model in
        % numerically safe manner for one line
        if SingleLineModel==1
            max_log_likelihoodL1 = max(sample_log_likelihoods_c4L1(quasar_ind, :));
            sample_probabilitiesL1 = ...
                exp(sample_log_likelihoods_c4L1(quasar_ind, :) - ...
                max_log_likelihoodL1);
            log_likelihoods_c4L1(quasar_ind, num_c4) = ...
                max_log_likelihoodL1 + log(mean(sample_probabilitiesL1)) ...
                - log(num_C4_samples)*(num_c4-1);;
              
            log_posteriors_c4L1(quasar_ind, num_c4) = ...
            log_priors_c4(quasar_ind, num_c4) + log_likelihoods_c4L1(quasar_ind, num_c4);
            
            fprintf(' ... log p(D | z_QSO,    L1)     : %0.2f\n', ...
                log_likelihoods_c4L1(quasar_ind, num_c4));
            fprintf(' ... log p(L1 | D, z_QSO)        : %0.2f\n', ...
                log_posteriors_c4L1(quasar_ind, num_c4));
        end
        
        % compute sample probabilities and log likelihood of DLA model in
        % numerically safe manner for  doublet 
        max_log_likelihoodL2 = max(sample_log_likelihoods_c4L2(quasar_ind, :, num_c4));
        sample_probabilitiesL2 = ...
            exp(sample_log_likelihoods_c4L2(quasar_ind, :, num_c4)  ... 
            - max_log_likelihoodL2);
        log_likelihoods_c4L2(quasar_ind, num_c4) = ...
            max_log_likelihoodL2 + log(mean(sample_probabilitiesL2))...
            - log(num_C4_samples)*(num_c4-1);
            
        
        log_posteriors_c4L2(quasar_ind, num_c4) = ...
            log_priors_c4(quasar_ind, num_c4) + log_likelihoods_c4L2(quasar_ind, num_c4);
        
        fprintf(' ... log p(D | z_QSO,    CIV)     : %0.2f\n', ...
            log_likelihoods_c4L2(quasar_ind, num_c4));
        fprintf(' ... log p(CIV | D, z_QSO)        : %0.2f\n', ...
            log_posteriors_c4L2(quasar_ind, num_c4));
        %  fprintf(' ... Num_CIV                      : %d\n ', ...
        %    Full_catalog.all_Num_c4_sys(quasar_ind))
        % fprintf('... FilterFlag                    : %d\n ', filter)
        if SingleLineModel==1
            [~, maxindL1] = nanmax(sample_log_likelihoods_c4L1(quasar_ind, :));
            map_z_c4L1(quasar_ind, num_c4 )    = sample_z_c4(maxindL1);        
            map_N_c4L1(quasar_ind, num_c4)  = log_nciv_samples(maxindL1);
            map_sigma_c4L1(quasar_ind, num_c4)  = sample_sigma_c4(maxindL1);
            fprintf('L1\nmap(N): %.2f, map(z_c4): %.2f, map(b/1e5): %.2f\n',map_N_c4L1(quasar_ind, num_c4),...
             map_z_c4L1(quasar_ind, num_c4), map_sigma_c4L1(quasar_ind, num_c4)/1e5);
        end

        [~, maxindL2] = nanmax(sample_log_likelihoods_c4L2(quasar_ind, :, num_c4));
        map_z_c4L2(quasar_ind, num_c4)    = sample_z_c4(maxindL2);        
        map_N_c4L2(quasar_ind, num_c4)  = log_nciv_samples(maxindL2);
        map_sigma_c4L2(quasar_ind, num_c4)  = sample_sigma_c4(maxindL2);
        fprintf('L2\nmap(N): %.2f, map(z_c4): %.2f, map(b/1e5): %.2f\n',...
        map_N_c4L2(quasar_ind, num_c4), map_z_c4L2(quasar_ind, num_c4),...
        map_sigma_c4L2(quasar_ind, num_c4)/1e5);

        max_log_posteriors = max([log_posteriors_no_c4(quasar_ind, num_c4), log_posteriors_c4L1(quasar_ind, num_c4), log_posteriors_c4L2(quasar_ind,num_c4)], [], 2);

        model_posteriors = ...
                exp(bsxfun(@minus, ...           
                [log_posteriors_no_c4(quasar_ind, num_c4), log_posteriors_c4L1(quasar_ind, num_c4), log_posteriors_c4L2(quasar_ind, num_c4)], ...
                max_log_posteriors));
        model_posteriors = ...
        bsxfun(@times, model_posteriors, 1 ./ sum(model_posteriors, 2));
    
        p_no_c4(quasar_ind, num_c4) = model_posteriors(1);
        p_c4L1(quasar_ind, num_c4)  = model_posteriors(2);
        p_c4(quasar_ind, num_c4)    = 1 - p_no_c4(quasar_ind, num_c4) -...
                                    p_c4L1(quasar_ind, num_c4);

        c4_pixel_ind1 = abs(this_wavelengths - (1+map_z_c4L2(quasar_ind, num_c4))*1548.2040)<3;
        c4_pixel_ind2 = abs(this_wavelengths - (1+map_z_c4L2(quasar_ind, num_c4))*1550.7810)<3;
        num_pixel_civ(quasar_ind, num_c4, 1) = nnz(c4_pixel_ind1);
        num_pixel_civ(quasar_ind, num_c4, 2) = nnz(c4_pixel_ind2);
        fprintf('CIV pixels:[%d, %d]\n', num_pixel_civ(quasar_ind, num_c4, :)); 
        if(plotting==1)
            % plotting
        
            this_ID = ID{quasar_ind};
            max_log_posteriors = max([log_posteriors_no_c4(quasar_ind, num_c4), log_posteriors_c4L1(quasar_ind, num_c4), log_posteriors_c4L2(quasar_ind,num_c4)], [], 2);
            num_lines=1;
            absorptionL1_fine= voigt_iP(finer(padded_wavelengths, nAVG),...
                            map_z_c4L1(quasar_ind, num_c4), 10^map_N_c4L1(quasar_ind, num_c4),...
                            num_lines, map_sigma_c4L1(quasar_ind, num_c4), finer(this_sigma_pixel, nAVG));
            absorptionL1 = Averager(absorptionL1_fine, nAVG, lenW_unmasked);
            absorptionL1 = absorptionL1(ind);
            c4_muL1    = this_mu     .* absorptionL1;

            num_lines=2;
            absorptionL2_fine= voigt_iP(finer(padded_wavelengths, nAVG),...
                map_z_c4L2(quasar_ind,  num_c4), 10^map_N_c4L2(quasar_ind, num_c4),...
                num_lines, map_sigma_c4L2(quasar_ind, num_c4), finer(this_sigma_pixel, nAVG));
            absorptionL2 = Averager(absorptionL2_fine, nAVG, lenW_unmasked);
            absorptionL2 = absorptionL2(ind);
            c4_muL2    = this_mu     .* absorptionL2;
            
            % Equivalent width calculation 
            EW(quasar_ind, num_c4) = trapz(this_wavelengths(ind), 1-absorptionL2);
            

            ttl = sprintf('ID:%s, zQSO:%.2f\n P(CIV)=%.2f, z_{CIV}=%.3f, P(S)=%.2f, z_S =%.3f\nz_{PM}=[%.3f, %.3f, %.3f, %.3f]\nEW=%.3f, N=%.2f, nPix=[%d, %d]', ...
                this_ID, z_qso, p_c4(quasar_ind, num_c4), map_z_c4L2(quasar_ind, num_c4), ...
                p_c4L1(quasar_ind, num_c4), map_z_c4L1(quasar_ind, num_c4), z_PM_test(quasar_ind,1:4),...
                EW(quasar_ind, num_c4), map_N_c4L2(quasar_ind, num_c4),...
                num_pixel_civ(quasar_ind, num_c4, 1), num_pixel_civ(quasar_ind, num_c4, 2));
            fid = sprintf('plt-vCut-3000-dv-mask-350/ind-%d-c4-%d.png', quasar_ind, num_c4);
            ind_zoomL2 = (abs(this_z_1548-map_z_c4L2(quasar_ind, num_c4))<2.5*kms_to_z(dv_mask)*(1+z_qso));
            ind_zoomL1 = (abs(this_z_1548-map_z_c4L1(quasar_ind, num_c4))<2.5*kms_to_z(dv_mask)*(1+z_qso));
            pltQSO(this_flux, this_wavelengths, c4_muL2, c4_muL1, ind_zoomL2, ind_zoomL1, ttl, fid)
        end
        

      
    end

    fprintf(' took %0.3fs.\n', toc);
   
end
% compute model posteriors in numerically safe manner



% save results
variables_to_save = {'release', 'training_set_name', ...
     'prior_ind', 'release', ...
    'test_ind', 'prior_z_qso_increase', ...
    'max_z_cut', 'min_z_c4s', 'max_z_c4s', ...
    'log_priors_no_c4', 'log_priors_c4', ...
    'log_likelihoods_no_c4',  ...
    'sample_log_likelihoods_c4L2', 'log_likelihoods_c4L2'...
    'log_posteriors_no_c4', 'log_posteriors_c4L1', 'log_posteriors_c4L2',...
    'model_posteriors', 'p_no_c4', 'p_c4L1' ...
    'map_z_c4L2', 'map_N_c4L2', 'map_sigma_c4L2', 'p_c4', 'EW'};

% filename = sprintf('%s/processed_qsos_tst_%s_S_%d_E_%d.mat', ...
%     processed_directory(release), ...
%     testing_set_name, ind_S, ind_E);
filename = sprintf('%s/processed_qsos_tst_6%s.mat', ...
    processed_directory(release), ...
    testing_set_name);

save(filename, variables_to_save{:}, '-v7.3');