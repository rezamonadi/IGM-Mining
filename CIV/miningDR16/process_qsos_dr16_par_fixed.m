% process_qsos: run CIV detection algorithm on specified objects

% Averaging     -> yes
% Single model  -> yes
% Multi CIV     -> yes
% Multi Singlet -> no
%
% removing CIV found region after each run
% also removing the map Z_civ
% load C4 catalog
PM_catalog = ...
    load('data/dr7/processed/catalog');
% train_ind -> those LOSs without any problem (filter_flag==0) that does not have
% Civ and useful for training null model (a model without Civ absorption line)
% prior_ind -> those LOSs with Civ absorption and in the half part of test
% test_ind -> second half without any filter flag for testing null and absorption model
% on the LOSs that we know have Civ or not. So, we asses our algorithm in this way.
all_QSO_ID_dr16_test = all_QSO_ID_dr16(test_ind);
if (ischar(prior_ind))
    prior_ind = eval(prior_ind);
end

% My prior_ind here is already those OK sight of lines that have CIV
prior.z_qsos  = PM_catalog.all_zqso(prior_ind);
prior.c4_ind = prior_ind;
prior.z_c4 = PM_catalog.all_z_civ3(prior_ind);

% filter out CIVs from prior catalog corresponding to region of spectrum below
% Ly-alpha QSO rest. In Roman's code, for detecting DLAs, instead of
% Ly-alpha they have Lyman-limit
for i = 1: numel(prior.z_c4)
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
z_qsos             =     all_zqso_dr16(test_ind);
all_num_quasars    =                numel(z_qsos);




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
log_priors_no_c4            = nan(num_quasars, max_civ);
log_priors_c4               = nan(num_quasars, max_civ);
log_likelihoods_no_c4       = nan(num_quasars, max_civ);
sample_log_likelihoods_c4L1 = nan(num_quasars, num_C4_samples, max_civ);
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
map_REW_1548_c4L1           = nan(num_quasars, max_civ);
map_REW_1548_c4L2           = nan(num_quasars, max_civ);
map_REW_1550_c4L2           = nan(num_quasars, max_civ);
p_c4                        = nan(num_quasars, max_civ);
p_c4L1                      = nan(num_quasars, max_civ);
p_no_c4                     = nan(num_quasars, max_civ);
REW_1548_dr16_samples       = nan(1, num_C4_samples, max_civ);
REW_1550_dr16_samples       = nan(1, num_C4_samples, max_civ);
% Compact L2-posterior representation for catalog construction.  Draws at
% the same position in each array refer to the same CIV parameter sample.
num_posterior_draws          = 256;
posterior_z_c4L2_draws       = nan(num_quasars, num_posterior_draws, max_civ);
posterior_logN_c4L2_draws    = nan(num_quasars, num_posterior_draws, max_civ);
posterior_sigma_c4L2_draws   = nan(num_quasars, num_posterior_draws, max_civ);
posterior_REW_1548_draws     = nan(num_quasars, num_posterior_draws, max_civ);
posterior_REW_1550_draws     = nan(num_quasars, num_posterior_draws, max_civ);
posterior_z_c4L2_16          = nan(num_quasars, max_civ);
posterior_z_c4L2_84          = nan(num_quasars, max_civ);
posterior_logN_c4L2_16       = nan(num_quasars, max_civ);
posterior_logN_c4L2_84       = nan(num_quasars, max_civ);
posterior_sigma_c4L2_16      = nan(num_quasars, max_civ);
posterior_sigma_c4L2_84      = nan(num_quasars, max_civ);
posterior_REW_1548_16        = nan(num_quasars, max_civ);
posterior_REW_1548_84        = nan(num_quasars, max_civ);
posterior_REW_1550_16        = nan(num_quasars, max_civ);
posterior_REW_1550_84        = nan(num_quasars, max_civ);
num_pixel_civ               = nan(num_quasars, max_civ, 2);
% A C IV proposal is evaluated only when its required transitions have
% enough surviving native pixels in the likelihood: L1 needs 1548 support;
% L2 needs support around both 1548 and 1550. This is a coverage criterion,
% not a cut on the measured REW doublet ratio. The default is enabled;
% process_qsos_dr16_par_fixed_precoverage.m provides the matched baseline
% without either coverage gate.
if ~exist('doublet_coverage_enabled', 'var')
    doublet_coverage_enabled = true;
end
doublet_coverage_velocity_kms = 300;
doublet_coverage_min_pixels = 3;
num_doublet_coverage_rejected = nan(num_quasars, max_civ);
num_singlet_coverage_rejected = nan(num_quasars, max_civ);
% min_sigma = 5e5;
sigma_civ_samples = min_sigma + (max_sigma-min_sigma)*offset_sigma_samples;

    % `finer(x, nAVG)` makes a fine grid for the Voigt/WDISP convolution.
    % The convolved profile is sampled at the native SDSS wavelengths below;
    % it is not box-averaged across a pixel.
fine_factor = nAVG + 1;
fine_width = width * fine_factor;


z_PM_prior = all_z_civ_C13(prior_ind,:);

this_quasar_ind = 0;
ind_E = num_quasars - 1 + ind_S;
if (ind_S==250001)
    ind_E = all_num_quasars;
    num_quasars = ind_E - ind_S + 1
end

for all_quasar_ind = ind_S:ind_E
    % if ~ismember(all_quasar_ind, plotted_qso_ind)
    %     continue;
    % end

    this_quasar_ind = this_quasar_ind + 1;

    tic;
    z_qso = z_qsos(all_quasar_ind);
    fprintf('processing quasar %i/%i (z_QSO = %0.4f) ...\n', ...
        all_quasar_ind, num_quasars, z_qso);

    this_wavelengths    =    all_wavelengths{all_quasar_ind};
    % this_wavelengths    =              this_wavelengths';
    this_flux           =           all_flux{all_quasar_ind};
    % this_flux           =                     this_flux';
    this_noise_variance = all_noise_variance{all_quasar_ind};
    % this_noise_variance =           this_noise_variance';
    this_pixel_mask     =     all_pixel_mask{all_quasar_ind};
    % this_pixel_mask     =  this_pixel_mask';
    this_sigma_pixel   = all_sigma_pixel{all_quasar_ind};
    % this_sigma_pixel    =               this_sigma_pixel';


    % convert to QSO rest frame
    this_rest_wavelengths = emitted_wavelengths(this_wavelengths, z_qso);
    unmasked_ind = (this_rest_wavelengths >= min_lambda) & ...
        (this_rest_wavelengths <= max_lambda) & (this_sigma_pixel>0);



    % keep complete copy of equally spaced wavelengths for absorption
    % computation
    this_unmasked_wavelengths = this_wavelengths(unmasked_ind);
    this_unmasked_sigma_pixel = this_sigma_pixel(unmasked_ind); % avoiding mismathed sizes for padded variavles and NaNs
    %     % [mask_ind] remove flux pixels with pixel_mask; pixel_mask is defined
    %     % in read_spec_dr7.m
    ind = unmasked_ind & (~this_pixel_mask);
    this_wavelengths      =      this_wavelengths(ind);
    this_rest_wavelengths = this_rest_wavelengths(ind);
    this_flux             =             this_flux(ind);
    this_noise_variance   =   this_noise_variance(ind);
    this_sigma_pixel      =      this_sigma_pixel(ind);

    % c4 existence prior
    less_ind = (prior.z_qsos < (z_qso + prior_z_qso_increase));
    less_systems = z_PM_prior(less_ind,:);

    this_num_quasars = nnz(less_ind);
    this_p_c4(1) = nnz(less_systems(:,1)>0 )/this_num_quasars;  % P(at least 1 CIV)
    for i=2:max_civ
        this_p_c4(i) = nnz(less_systems(:,i)>0)/nnz(less_systems(:,i-1)>0); % P(at least n CIV | at least n-1 CIV)
        if (this_p_c4(i)==0)
            this_p_c4(i) = 0.001;
        end
    end

    fprintf('\n');
    for i = 1:max_civ

        log_priors_no_c4(this_quasar_ind, i) = log(1 - this_p_c4(i));
        fprintf(' ...     p(%i  CIVs | z_QSO)       : %0.3f\n', i, this_p_c4(i));
        fprintf(' ...     p(no CIV  | z_QSO)       : %0.3f\n', exp(log_priors_no_c4(this_quasar_ind, i)) );
    end
    log_priors_c4(this_quasar_ind,:) = log(this_p_c4(:));

    % interpolate model onto given wavelengths
    this_mu = mu_interpolator( this_rest_wavelengths);
    this_M  =  M_interpolator({this_rest_wavelengths, 1:k});

    min_z_c4s(this_quasar_ind) = min_z_c4(this_wavelengths, z_qso);
    max_z_c4s(this_quasar_ind) = max_z_c4(z_qso, max_z_cut);

    sample_z_c4 = ...
        min_z_c4s(this_quasar_ind) +  ...
        (max_z_c4s(this_quasar_ind) - min_z_c4s(this_quasar_ind)) * offset_z_samples;



    % ensure enough pixels are on either side for convolving with
    % instrument profile and sigma_pixels

    % padded_wavelengths = ...
    %     [logspace(log10(min(this_unmasked_wavelengths)) - width * pixel_spacing, ...
    %     log10(min(this_unmasked_wavelengths)) - pixel_spacing,...
    %     width)';...
    %     this_unmasked_wavelengths;...
    %     logspace(log10(max(this_unmasked_wavelengths)) + pixel_spacing,...
    %     log10(max(this_unmasked_wavelengths)) + width * pixel_spacing,...
    %     width)'...
    %     ];


    % padded_sigma_pixel = ...
    %     [this_sigma_pixel(1)*ones(width,1);...
    %     this_unmasked_sigma_pixel;...
    %     this_sigma_pixel(end)*ones(width,1)];

    this_unmasked_wavelengths_fine = finer(this_unmasked_wavelengths, nAVG)';
    padded_wavelengths_fine = ...
        [logspace(log10(min(this_unmasked_wavelengths)) - width * pixel_spacing, ...
        log10(min(this_unmasked_wavelengths)) - pixel_spacing/fine_factor,...
        fine_width)';...
        this_unmasked_wavelengths_fine;...
        logspace(log10(max(this_unmasked_wavelengths)) + pixel_spacing/fine_factor,...
        log10(max(this_unmasked_wavelengths)) + width * pixel_spacing,...
        fine_width)'...
        ];

    % Keep this array aligned with this_unmasked_wavelengths.  The latter
    % retains masked pixels so that the absorption profile can subsequently
    % be reduced with `ind`; this_sigma_pixel no longer has those pixels.
    padded_sigma_pixels_fine = ...
        [this_unmasked_sigma_pixel(1)*ones(fine_width,1);...
        finer(this_unmasked_sigma_pixel, nAVG)';...
        this_unmasked_sigma_pixel(end)*ones(fine_width,1)];


    % [mask_ind] to retain only unmasked pixels from computed absorption profile
    % this has to be done by using the unmasked_ind which has not yet
    % been applied this_pixel_mask.
    ind = (~this_pixel_mask(unmasked_ind));

    % compute probabilities under DLA model for each of the sampled
    % (normalized offset, log(N HI)) pairs
    lenW_unmasked = length(this_unmasked_wavelengths);
    % The MEX trims the 3*fine_factor convolution-padding samples from both
    % ends, so its profile is on this_unmasked_wavelengths_fine. These indices
    % select the native SDSS wavelengths from that returned fine-grid profile.
    fine_native_ind = 1:fine_factor:(1 + (lenW_unmasked - 1) * fine_factor);
    ind_not_remove = true(size(this_flux));
    for num_c4=1:max_civ
        fprintf('num_civ:%d\n',num_c4);
        this_z_1548 = (this_wavelengths / civ_1548_wavelength) - 1;
        this_z_1550 = (this_wavelengths / civ_1550_wavelength) - 1;
        % if(num_c4>1)
        %     % if((p_c4(this_quasar_ind, num_c4-1)>p_c4L1(this_quasar_ind, num_c4-1)) & ...
        %     %    (p_c4(this_quasar_ind, num_c4-1)>p_no_c4(this_quasar_ind, num_c4-1)))
        %     %     ind_not_remove = ind_not_remove  & ...
        %     %         (abs(this_z_1548 - map_z_c4L2(this_quasar_ind, num_c4-1))>kms_to_z(dv_mask)*(1+map_z_c4L2(this_quasar_ind, num_c4-1))) & ...
        %     %         (abs(this_z_1550 - map_z_c4L2(this_quasar_ind, num_c4-1))>kms_to_z(dv_mask)*(1+map_z_c4L2(this_quasar_ind, num_c4-1)));
        %     % end
        %     if (p_c4(this_quasar_ind, num_c4-1)>0.75)
        %         ind_not_remove = ind_not_remove  & ...
        %             (abs(this_z_1548 - map_z_c4L2(this_quasar_ind, num_c4-1))>kms_to_z(dv_mask)*(1+map_z_c4L2(this_quasar_ind, num_c4-1))) & ...
        %             (abs(this_z_1550 - map_z_c4L2(this_quasar_ind, num_c4-1))>kms_to_z(dv_mask)*(1+map_z_c4L2(this_quasar_ind, num_c4-1)));
        %     end
        %     %
        %     % if((p_c4L1(this_quasar_ind, num_c4-1)>p_c4(this_quasar_ind, num_c4-1)) & ...
        %     %             (p_c4L1(this_quasar_ind, num_c4-1)>p_no_c4(this_quasar_ind, num_c4-1)))
        %     %     ind_not_remove = ind_not_remove  & ...
        %     %     (abs(this_z_1548 - map_z_c4L1(this_quasar_ind, num_c4-1))>kms_to_z(dv_mask)*(1+map_z_c4L1(this_quasar_ind, num_c4-1)));
        %     % end



        %      if (p_c4L1(this_quasar_ind, num_c4-1)>0.75)
        %         ind_not_remove = ind_not_remove  & ...
        %         (abs(this_z_1548 - map_z_c4L1(this_quasar_ind, num_c4-1))>kms_to_z(dv_mask)*(1+map_z_c4L1(this_quasar_ind, num_c4-1)));
        %     end


        %     if((p_no_c4(this_quasar_ind, num_c4-1)>p_c4(this_quasar_ind, num_c4-1)) & ...
        %         (p_no_c4(this_quasar_ind, num_c4-1)>p_c4L1(this_quasar_ind, num_c4-1)))

        %         fprintf('No more than %d CIVs in this spectrum.\n', num_c4-1)
        %         break;
        %     end

        % end

        if num_c4 > 1
            previous = num_c4 - 1;

            pDoublet = p_c4(this_quasar_ind, previous);
            pSinglet = p_c4L1(this_quasar_ind, previous);
            pNull    = p_no_c4(this_quasar_ind, previous);

            % Stop if the null model is preferred.
            if pNull >= max(pDoublet, pSinglet)
                fprintf('Stopping: null model preferred at search %d.\n', ...
                    previous);
                break;
            end

            oldMask = ind_not_remove;

            if pDoublet > 0.75
                zPrevious = map_z_c4L2(this_quasar_ind, previous);
                dz = kms_to_z(dv_mask) * (1 + zPrevious);

                ind_not_remove = ind_not_remove & ...
                    (abs(this_z_1548 - zPrevious) > dz) & ...
                    (abs(this_z_1550 - zPrevious) > dz);

            elseif pSinglet > 0.75
                zPrevious = map_z_c4L1(this_quasar_ind, previous);
                dz = kms_to_z(dv_mask) * (1 + zPrevious);

                ind_not_remove = ind_not_remove & ...
                    (abs(this_z_1548 - zPrevious) > dz);

            else
                % No confident feature to mask: do not repeat the same search.
                fprintf(['Stopping: neither absorber model exceeds 0.75 ' ...
                    'at search %d.\n'], previous);
                break;
            end

            % Ensure the next search actually uses different data.
            if isequal(ind_not_remove, oldMask)
                fprintf('Stopping: masking removed no additional pixels.\n');
                break;
            end

            if ~any(ind_not_remove)
                fprintf('Stopping: no usable pixels remain.\n');
                break;
            end
        end

        log_likelihoods_no_c4(this_quasar_ind, num_c4) = ...
            log_mvnpdf_low_rank(this_flux(ind_not_remove), this_mu(ind_not_remove),...
            this_M(ind_not_remove, :), this_noise_variance(ind_not_remove));

        log_posteriors_no_c4(this_quasar_ind, num_c4) = ...
            log_priors_no_c4(this_quasar_ind, num_c4) + log_likelihoods_no_c4(this_quasar_ind, num_c4);

        fprintf(' ... log p(D | z_QSO, no CIV)     : %0.2f\n', ...
            log_likelihoods_no_c4(this_quasar_ind, num_c4));
        fprintf(' ... log p(no CIV | D, z_QSO)     : %0.2f\n', ...
            log_posteriors_no_c4(this_quasar_ind, num_c4));

        % `this_wavelengths` contains only pixels which survived the
        % preload mask. `ind_not_remove` additionally removes pixels
        % assigned to an earlier confident absorber on later searches.
        % Require local support for both lines in the data used here.
        doublet_coverage_valid = false(1, num_C4_samples);
        singlet_coverage_valid = false(1, num_C4_samples);

        parfor i = 1:num_C4_samples
            % Limitting red-shift in the samples
            % absorption corresponding to this sample with one absorption line as a noise model

            % absorption corresponding to this sample with two absorption lines as a doublet
            % compute fine absorption
            if doublet_coverage_enabled
                coverage_dz = kms_to_z(doublet_coverage_velocity_kms) * ...
                    (1 + sample_z_c4(i));
                has_usable_1548_pixels = nnz(ind_not_remove & ...
                    abs(this_z_1548 - sample_z_c4(i)) <= coverage_dz) >= ...
                    doublet_coverage_min_pixels;
                has_usable_1550_pixels = nnz(ind_not_remove & ...
                    abs(this_z_1550 - sample_z_c4(i)) <= coverage_dz) >= ...
                    doublet_coverage_min_pixels;
                doublet_coverage_valid(i) = has_usable_1548_pixels && ...
                    has_usable_1550_pixels;
                singlet_coverage_valid(i) = has_usable_1548_pixels;
            else
                doublet_coverage_valid(i) = true;
                singlet_coverage_valid(i) = true;
            end
            if doublet_coverage_valid(i)
                absorptionL2_fine = voigt_iP_transition(padded_wavelengths_fine, sample_z_c4(i), ...
                    nciv_samples(i), 3 , sigma_civ_samples(i), ...
                    padded_sigma_pixels_fine, fine_factor);

                % Sample the fine convolved profile at native wavelengths. Do not
                % apply an additional pixel box average after the WDISP convolution.
                absorptionL2 = absorptionL2_fine(fine_native_ind);
                absorptionL2 = absorptionL2(ind);
                absorptionL2(~isfinite(absorptionL2)) = 1;
                c4_muL2     = this_mu     .* absorptionL2;
                c4_ML2      = this_M      .* absorptionL2;

                sample_log_likelihoods_c4L2(this_quasar_ind, i, num_c4) = ...
                    log_mvnpdf_low_rank(this_flux(ind_not_remove),...
                    c4_muL2(ind_not_remove),...
                    c4_ML2(ind_not_remove, :), ...
                    this_noise_variance(ind_not_remove));


                 %REW 1550 
            
                absorption_1550_fine = voigt_iP_transition(padded_wavelengths_fine, sample_z_c4(i), ...
                    nciv_samples(i), 2 , sigma_civ_samples(i), ...
                    padded_sigma_pixels_fine, fine_factor);
                REW_1550_dr16_samples(1, i, num_c4) = ...
                    trapz(this_unmasked_wavelengths_fine, 1-absorption_1550_fine) / ...
                    (1+sample_z_c4(i));

            else
                sample_log_likelihoods_c4L2(this_quasar_ind, i, num_c4) = -Inf;
            end

            if singlet_coverage_valid(i)
                absorptionL1_fine = voigt_iP_transition(padded_wavelengths_fine, sample_z_c4(i), ...
                nciv_samples(i), 1 , sigma_civ_samples(i), ...
                padded_sigma_pixels_fine, fine_factor);

            


                REW_1548_dr16_samples(1, i, num_c4) = ...
                    trapz(this_unmasked_wavelengths_fine, 1-absorptionL1_fine) / ...
                    (1+sample_z_c4(i));
                absorptionL1 = absorptionL1_fine(fine_native_ind);
                absorptionL1 = absorptionL1(ind);
                absorptionL1(~isfinite(absorptionL1))=1;
                c4_muL1     = this_mu     .* absorptionL1;
                c4_ML1      = this_M      .* absorptionL1;



                sample_log_likelihoods_c4L1(this_quasar_ind, i, num_c4) = ...
                    log_mvnpdf_low_rank(this_flux(ind_not_remove),...
                    c4_muL1(ind_not_remove), c4_ML1(ind_not_remove, :), ...
                    this_noise_variance(ind_not_remove));
            else

                sample_log_likelihoods_c4L1(this_quasar_ind, i, num_c4) = -Inf;
            end

           
        end

        num_doublet_coverage_rejected(this_quasar_ind, num_c4) = ...
            nnz(~doublet_coverage_valid);


        % compute sample probabilities and log likelihood of DLA model in
        % numerically safe manner for one line

        max_log_likelihoodL1 = max(sample_log_likelihoods_c4L1(this_quasar_ind, :, num_c4));
        sample_probabilitiesL1 = ...
            exp(sample_log_likelihoods_c4L1(this_quasar_ind, :, num_c4) - ...
            max_log_likelihoodL1);
        log_likelihoods_c4L1(this_quasar_ind, num_c4) = ...
            max_log_likelihoodL1 + log(mean(sample_probabilitiesL1));

        log_posteriors_c4L1(this_quasar_ind, num_c4) = ...
            log_priors_c4(this_quasar_ind, num_c4) + log_likelihoods_c4L1(this_quasar_ind, num_c4);

        fprintf(' ... log p(D | z_QSO,    L1)     : %0.2f\n', ...
            log_likelihoods_c4L1(this_quasar_ind, num_c4));
        fprintf(' ... log p(L1 | D, z_QSO)        : %0.2f\n', ...
            log_posteriors_c4L1(this_quasar_ind, num_c4));


        % compute sample probabilities and log likelihood of DLA model in
        % numerically safe manner for  doublet
        max_log_likelihoodL2 = max(sample_log_likelihoods_c4L2(this_quasar_ind, :, num_c4));
        sample_probabilitiesL2 = ...
            exp(sample_log_likelihoods_c4L2(this_quasar_ind, :, num_c4)  ...
            - max_log_likelihoodL2);
        log_likelihoods_c4L2(this_quasar_ind, num_c4) = ...
            max_log_likelihoodL2 + log(mean(sample_probabilitiesL2));



        log_posteriors_c4L2(this_quasar_ind, num_c4) = ...
            log_priors_c4(this_quasar_ind, num_c4) + log_likelihoods_c4L2(this_quasar_ind, num_c4);

        fprintf(' ... log p(D | z_QSO,    CIV)     : %0.2f\n', ...
            log_likelihoods_c4L2(this_quasar_ind, num_c4));
        fprintf(' ... log p(CIV | D, z_QSO)        : %0.2f\n', ...
            log_posteriors_c4L2(this_quasar_ind, num_c4));
        %  fprintf(' ... Num_CIV                      : %d\n ', ...
        %    Full_catalog.all_Num_c4_sys(this_quasar_ind))
        % fprintf('... FilterFlag                    : %d\n ', filter)
        [~, maxindL1] = max(sample_log_likelihoods_c4L1(this_quasar_ind, :, num_c4), [], 'omitnan');
        map_z_c4L1(this_quasar_ind, num_c4 )    = sample_z_c4(maxindL1);
        map_N_c4L1(this_quasar_ind, num_c4)  = log_nciv_samples(maxindL1);
        map_sigma_c4L1(this_quasar_ind, num_c4)  = sigma_civ_samples(maxindL1);
        map_REW_1548_c4L1(this_quasar_ind, num_c4) = ...
            REW_1548_dr16_samples(1, maxindL1, num_c4);
        fprintf('L1\nmap(N): %.2f, map(z_c4): %.2f, map(b/1e5): %.2f\n',map_N_c4L1(this_quasar_ind, num_c4),...
            map_z_c4L1(this_quasar_ind, num_c4), map_sigma_c4L1(this_quasar_ind, num_c4)/1e5);

        [~, maxindL2] = max(sample_log_likelihoods_c4L2(this_quasar_ind, :, num_c4),[], 'omitnan');
        map_z_c4L2(this_quasar_ind, num_c4)    = sample_z_c4(maxindL2);
        map_N_c4L2(this_quasar_ind, num_c4)  = log_nciv_samples(maxindL2);
        map_sigma_c4L2(this_quasar_ind, num_c4)  = sigma_civ_samples(maxindL2);
        map_REW_1548_c4L2(this_quasar_ind, num_c4) = ...
            REW_1548_dr16_samples(1, maxindL2, num_c4);
        map_REW_1550_c4L2(this_quasar_ind, num_c4) = ...
            REW_1550_dr16_samples(1, maxindL2, num_c4);

        fprintf('L2\nmap(N): %.2f, map(z_c4): %.2f, map(b/1e5): %.2f\n',...
            map_N_c4L2(this_quasar_ind, num_c4), map_z_c4L2(this_quasar_ind, num_c4),...
            map_sigma_c4L2(this_quasar_ind, num_c4)/1e5);

        max_log_posteriors = max([log_posteriors_no_c4(this_quasar_ind, num_c4), log_posteriors_c4L1(this_quasar_ind, num_c4), log_posteriors_c4L2(this_quasar_ind,num_c4)], [], 2);

        model_posteriors = ...
            exp(bsxfun(@minus, ...
            [log_posteriors_no_c4(this_quasar_ind, num_c4), log_posteriors_c4L1(this_quasar_ind, num_c4), log_posteriors_c4L2(this_quasar_ind, num_c4)], ...
            max_log_posteriors));
        model_posteriors = ...
            bsxfun(@times, model_posteriors, 1 ./ sum(model_posteriors, 2));

        p_no_c4(this_quasar_ind, num_c4) = model_posteriors(1);
        p_c4L1(this_quasar_ind, num_c4)  = model_posteriors(2);
        p_c4(this_quasar_ind, num_c4)    = 1 - p_no_c4(this_quasar_ind, num_c4) -...
            p_c4L1(this_quasar_ind, num_c4);

        % Draw a compact, joint representation of the L2 posterior.  The
        % likelihood-only weights are appropriate because the CIV proposals
        % are drawn from the prior.  Keeping matched indices preserves the
        % z, logN, sigma, REW1548, and REW1550 correlations.
        posterior_weightsL2 = sample_probabilitiesL2(:);
        posterior_weightsL2 = posterior_weightsL2 / sum(posterior_weightsL2);
        posterior_draw_ind = randsample(num_C4_samples, num_posterior_draws, ...
            true, posterior_weightsL2);

        posterior_z_draws = reshape(sample_z_c4(posterior_draw_ind), 1, []);
        posterior_logN_draws = reshape(log_nciv_samples(posterior_draw_ind), 1, []);
        posterior_sigma_draws = reshape(sigma_civ_samples(posterior_draw_ind), 1, []);
        posterior_REW_1548 = reshape(REW_1548_dr16_samples(1, ...
            posterior_draw_ind, num_c4), 1, []);
        posterior_REW_1550 = reshape(REW_1550_dr16_samples(1, ...
            posterior_draw_ind, num_c4), 1, []);

        posterior_z_c4L2_draws(this_quasar_ind, :, num_c4) = posterior_z_draws;
        posterior_logN_c4L2_draws(this_quasar_ind, :, num_c4) = posterior_logN_draws;
        posterior_sigma_c4L2_draws(this_quasar_ind, :, num_c4) = posterior_sigma_draws;
        posterior_REW_1548_draws(this_quasar_ind, :, num_c4) = posterior_REW_1548;
        posterior_REW_1550_draws(this_quasar_ind, :, num_c4) = posterior_REW_1550;

        posterior_z_c4L2_16(this_quasar_ind, num_c4) = prctile(posterior_z_draws, 16);
        posterior_z_c4L2_84(this_quasar_ind, num_c4) = prctile(posterior_z_draws, 84);
        posterior_logN_c4L2_16(this_quasar_ind, num_c4) = prctile(posterior_logN_draws, 16);
        posterior_logN_c4L2_84(this_quasar_ind, num_c4) = prctile(posterior_logN_draws, 84);
        posterior_sigma_c4L2_16(this_quasar_ind, num_c4) = prctile(posterior_sigma_draws, 16);
        posterior_sigma_c4L2_84(this_quasar_ind, num_c4) = prctile(posterior_sigma_draws, 84);
        posterior_REW_1548_16(this_quasar_ind, num_c4) = prctile(posterior_REW_1548, 16);
        posterior_REW_1548_84(this_quasar_ind, num_c4) = prctile(posterior_REW_1548, 84);
        posterior_REW_1550_16(this_quasar_ind, num_c4) = prctile(posterior_REW_1550, 16);
        posterior_REW_1550_84(this_quasar_ind, num_c4) = prctile(posterior_REW_1550, 84);

        % c4_pixel_ind1 = abs(this_wavelengths - (1+map_z_c4L2(this_quasar_ind, num_c4))*1548.2040)<3;
        % c4_pixel_ind2 = abs(this_wavelengths - (1+map_z_c4L2(this_quasar_ind, num_c4))*1550.7810)<3;
        % num_pixel_civ(this_quasar_ind, num_c4, 1) = nnz(c4_pixel_ind1);
        % num_pixel_civ(this_quasar_ind, num_c4, 2) = nnz(c4_pixel_ind2);
        % fprintf('CIV pixels:[%d, %d]\n', num_pixel_civ(this_quasar_ind, num_c4, :));


    
        % fprintf('REW(%d,%d)=%e\n', this_quasar_ind, num_c4, REW_1548_dr16_samples(this_quasar_ind, num_c4));

        % plotting
        if(plotting==1)

            max_log_posteriors = max([log_posteriors_no_c4(this_quasar_ind, num_c4),...
                log_posteriors_c4L1(this_quasar_ind, num_c4),...
                log_posteriors_c4L2(this_quasar_ind,num_c4)], [], 2);
            
            absorptionL1_fine= voigt_iP_transition(padded_wavelengths_fine,...
                map_z_c4L1(this_quasar_ind, num_c4), 10^map_N_c4L1(this_quasar_ind, num_c4),...
                1, map_sigma_c4L1(this_quasar_ind, num_c4), ...
                padded_sigma_pixels_fine, fine_factor);
                
            absorptionL1 = absorptionL1_fine(fine_native_ind);
            absorptionL1 = absorptionL1(ind);
            c4_muL1    = this_mu     .* absorptionL1;
            c4_ML1  = this_M  .* absorptionL1;

                        
            absorptionL2_fine= voigt_iP_transition(padded_wavelengths_fine,...
                map_z_c4L2(this_quasar_ind,  num_c4), 10^map_N_c4L2(this_quasar_ind, num_c4),...
                3, map_sigma_c4L2(this_quasar_ind, num_c4), ...
                padded_sigma_pixels_fine, fine_factor);



            absorptionL2 = absorptionL2_fine(fine_native_ind);
            absorptionL2 = absorptionL2(ind);
            c4_muL2    = this_mu     .* absorptionL2;
            c4_ML2  = this_M  .* absorptionL2;
            % posterior mean GP model flux model, conditioned  on the current observed flux  
            model_post_L2 = mean_posterior_continuum( ...
                    this_flux(ind_not_remove), ...
                    c4_muL2(ind_not_remove), ...
                    c4_ML2(ind_not_remove,:), ...
                    this_noise_variance(ind_not_remove));

            model_post_L1 = mean_posterior_continuum( ...
                    this_flux(ind_not_remove), ...
                    c4_muL1(ind_not_remove), ...
                    c4_ML1(ind_not_remove,:), ...
                    this_noise_variance(ind_not_remove));

            ttl = sprintf('ID:%s, zQSO:%.2f\n P(CIV)=%.2f, z_{L2}=%.5f, P(S)=%.2f, z_{L1} = %.5f, REW=%.3f, SN=%.2f', ...
                all_QSO_ID_dr16_test{all_quasar_ind}, z_qso, p_c4(this_quasar_ind, num_c4), ...
                map_z_c4L2(this_quasar_ind, num_c4),p_c4L1(this_quasar_ind, num_c4), map_z_c4L1(this_quasar_ind, num_c4) ,map_REW_1548_c4L2(this_quasar_ind, num_c4),...
                median(this_flux./sqrt(this_noise_variance)))

            fid = sprintf('plt-bad-DR/ind-%d-c4-%d.png',  all_quasar_ind, num_c4);
            % ind_zoomL2 = (abs(this_z_1548-map_z_c4L2(this_quasar_ind, num_c4))<5*kms_to_z(dv_mask)*(1+map_z_c4L2(all_quasar_ind, num_c4)));
            % ind_zoomL1 = (abs(this_z_1548-map_z_c4L1(this_quasar_ind, num_c4))<5*kms_to_z(dv_mask)*(1+map_z_c4L1(all_quasar_ind, num_c4)));
            fprintf('min(Z):%.3f, max(Z):%.3f\n', min(sample_z_c4), max(sample_z_c4))
            indMAP = (abs(sigma_civ_samples - map_sigma_c4L2(this_quasar_ind, num_c4))<25e5);
            fprintf('min(Z(indMAP)):%.3f, max(Z(indMAP)):%.3f\n', min(sample_z_c4(indMAP)), max(sample_z_c4(indMAP)))

            pltQSOv2(this_flux(ind_not_remove), this_wavelengths(ind_not_remove), model_post_L2, model_post_L1, ttl, fid)


        end



    end

    fprintf(' took %0.3fs per qso.\n', toc/num_c4);

end
% % compute model posteriors in numerically safe manner


%  save results
if saving==1
    variables_to_save = {'training_set_name', ...
        'prior_ind',  ...
        'test_ind', 'prior_z_qso_increase', ...
        'max_z_cut', 'min_z_c4s', 'max_z_c4s', ...
        'log_priors_no_c4', 'log_priors_c4', ...
        'log_likelihoods_no_c4',  ...
        'sample_log_likelihoods_c4L2', ...
        'sample_log_likelihoods_c4L1',...
        'log_likelihoods_c4L2',...
        'log_likelihoods_c4L1',...
        'log_posteriors_no_c4', 'log_posteriors_c4L1', 'log_posteriors_c4L2',...
        'p_no_c4', 'p_c4L1', ...
        'map_z_c4L2', 'map_N_c4L2', 'map_sigma_c4L2', 'p_c4',...
        'map_REW_1548_c4L1', 'map_REW_1548_c4L2', 'map_REW_1550_c4L2', ...
        'map_z_c4L1', 'map_N_c4L1', 'map_sigma_c4L1', ...
        'num_posterior_draws', ...
        'posterior_z_c4L2_draws', 'posterior_logN_c4L2_draws', 'posterior_sigma_c4L2_draws', ...
        'posterior_REW_1548_draws', 'posterior_REW_1550_draws', ...
        'posterior_z_c4L2_16', 'posterior_z_c4L2_84', ...
        'posterior_logN_c4L2_16', 'posterior_logN_c4L2_84', ...
        'posterior_sigma_c4L2_16', 'posterior_sigma_c4L2_84', ...
        'posterior_REW_1548_16', 'posterior_REW_1548_84', ...
        'posterior_REW_1550_16', 'posterior_REW_1550_84', ...
        'doublet_coverage_enabled', 'doublet_coverage_velocity_kms', 'doublet_coverage_min_pixels', ...
        'num_doublet_coverage_rejected'};

    if doublet_coverage_enabled
        coverage_output_tag = '';
    else
        coverage_output_tag = '_precoverage';
    end
    filename = sprintf(['data/dr16/processed/processed_qsos_tst_%s%s_S_%d_E_%d.mat'], ...
                testing_set_name, coverage_output_tag, ind_S, ind_E);

    save(filename, variables_to_save{:}, '-v7.3');
end

