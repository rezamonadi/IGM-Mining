
% Comparing column density and redshift 
% from our MAP analysis to C13 measurements
clc 
clear
set_parameters_dr7;
build_catalog_dr7;
filename ='data/dr7/processed/processed_qsos_tst_mask-1-prior-1-OccamRazor-1-nC4-30000-plt-0-MaskinP-0-fixedPr.mat';
variables_to_load = {'prior_ind', 'release', ...
   'test_ind', 'prior_z_qso_increase', ...
   'max_z_cut', 'min_z_c4s', 'max_z_c4s', ...
   'log_priors_no_c4', 'log_priors_c4', ...
    'log_likelihoods_no_c4',  ...
    'sample_log_likelihoods_c4L2', 'log_likelihoods_c4L2'...
    'log_posteriors_no_c4', 'log_posteriors_c4L1', 'log_posteriors_c4L2',...
    'model_posteriors', 'p_no_c4', ...
    'map_z_c4L2', 'map_N_c4L2', 'map_sigma_c4L2', 'p_c4'};
load(filename, variables_to_load{:});
load('data/C4_catalogs/Cooksey_C4_cat/processed/CIV-cat.mat','c4_QSO_ID','Z_c4','NCIV');

num_quasars = sum(test_ind);
ID = all_QSO_ID(test_ind);
NCIV_Cooksey_compare = [0];
jj=0;
for quasar_ind=1:num_quasars
    all_c4_NCIV_test=all_N_civ(test_ind,:);
    all_c4_Z_test=all_z_civ(test_ind, :);
    this_c4s = all_c4_NCIV_test(quasar_ind,:); % all of the Ns for 17 possible CIV systems in C13 for the current test QSO
    this_Zs  = all_c4_Z_test(quasar_ind,:);    % all of the Zs for 17 possible CIV systems in C13 for the current test QSO
    
    matched_Z_c4= this_Zs(abs(this_Zs -map_z_c4L2(quasar_ind))==min(abs(this_Zs-map_z_c4L2(quasar_ind))));
    % Redshift of the closest absorber in the C13 to our found absorber in redshift space
    matched_N_c4= this_c4s(abs(this_Zs -map_z_c4L2(quasar_ind))==min(abs(this_Zs-map_z_c4L2(quasar_ind))));
    % N of the closest absorber in the C13 to our found absorber in redshift space
    this_ID = ID{quasar_ind};
    num_this_systems = nnz(all_c4_NCIV_test(quasar_ind,:)>0); % total number of detected absorbers in C13 for the current test QSO
    
    if(num_this_systems>0)
        jj=jj+1;
        NCIV_Cooksey_compare(jj) = matched_N_c4(1) - map_N_c4L2(quasar_ind); % N comparision
        Z_Cooksey_compare(jj) = matched_Z_c4(1) - map_z_c4L2(quasar_ind);    % z comparison
        
    end
end
fig=figure();
h=histogram(NCIV_Cooksey_compare, 20);
set(get(gca, 'XLabel'), 'String', 'N^{CIV}_{C13} - N^{CIV}_{MAP}');
set(get(gca, 'YLabel'), 'String', 'Frequency');
set(gca, 'FontSize',20)
h.FaceColor = [0,191/255,1]; % changing the Face Color of histogram
h.EdgeColor = [0,0,0];       % changing the Edge Color of histogram
exportgraphics(fig, sprintf('%s-N(C13)-N(MAP).pdf',...
                            training_set_name),...
                            'ContentType','vector')
fig=figure();
h = histogram(Z_Cooksey_compare, 15);
set(get(gca, 'XLabel'), 'String', 'z^{CVI}_{C13}-z^{CIV}_{MAP}');
set(get(gca, 'YLabel'), 'String', 'Frequency');
set(gca, 'FontSize',20)
h.FaceColor = [0,191/255,1]; % changing the Face Color of histogram
h.EdgeColor = [0,0,0];       % changing the Edge Color of histogram


exportgraphics(fig, sprintf('%s-z(C13)-z(MAP).pdf',...
                            training_set_name),'ContentType','vector')

