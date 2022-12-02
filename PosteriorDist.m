
filename = 'data/dr12/processed/processed_dr12.mat';
fprintf('loading...\n');
load(filename);
% loading full 3D poster
sample_log_likelihoods_c4L2    = savingCat.all_sample_log_likelihoods_c4L2;
map_z_c4L2     = savingCat.all_map_z_c4L2;
map_N_c4L2     = savingCat.all_map_N_c4L2;
map_sigma_c4L2 = savingCat.all_map_sigma_c4L2;  
p_c4           = savingCat.all_p_c4; 
z_qsos             =     all_zqso_dr12(test_ind);
test_ind       = savingCat.test_ind;
num_quasars = nnz(test_ind)
for this_quasar_ind=1:30
    fprintf('QSO-%d\n',this_quasar_ind);
    z_qso = z_qsos(this_quasar_ind);
    min_z_c4s(this_quasar_ind) = min_z_c4(1310, z_qso);
    max_z_c4s(this_quasar_ind) = max_z_c4(z_qso, max_z_cut);
    binZ = (max_z_c4s(this_quasar_ind) - min_z_c4s(this_quasar_ind))/500;
    sample_z_c4 = ...
        min_z_c4s(this_quasar_ind) +  ...
        (max_z_c4s(this_quasar_ind) - min_z_c4s(this_quasar_ind)) * offset_z_samples;
    for num_c4=1:7
        if p_c4(this_quasar_ind, num_c4)>0.85
            fig = figure();
            weight = sample_log_likelihoods_c4L2(this_quasar_ind, :, num_c4);
            % weight = weight - min(weight);
            % weight = weight/sum(weight);

            for ii=1:500
                this_min_z = min_z_c4s(this_quasar_ind) + (ii-1)*binZ;
                this_max_z = min_z_c4s(this_quasar_ind) + ii*binZ;
                indMean = (sample_z_c4>= this_min_z) & (sample_z_c4<this_max_z);
                mean_weight(ii) = mean(weight(indMean));
                mean_sample_z_c4(ii)  =  (this_max_z+this_min_z)*0.5;
            end
            % [x_sort, ind_x_sort] = sort(mean_sample_z_c4);
            % y_sort = mean_weight(ind_x_sort);
            % plot (x_sort, y_sort)
            % hold on
            % xline(map_z_c4L2(this_quasar_ind, num_c4))
            % hold on

            indPlot = abs(mean_sample_z_c4-map_z_c4L2(this_quasar_ind, num_c4))<0.01;
            x = mean_sample_z_c4(indPlot);
            y = mean_weight(indPlot); 

            [x_sort, ind_x_sort] = sort(x);
            y_sort = y(ind_x_sort);
            % axes('position', [0.2 0.68 0.18 0.15]);
            % box on 
            p = plot(x_sort, y_sort);
            p.Color = [0.8500 0.3250 0.0980, 0.6]; 
            p.LineWidth=1;
            hold on
            xline(map_z_c4L2(this_quasar_ind, num_c4))
            exportgraphics(fig, sprintf('ConfIntvl/Z-Lik-QSO-%d-c4-%dmean500_nobx.png', this_quasar_ind, num_c4), 'Resolution', 800)
            close all
        end
    end
end


















% % 2D histograms 

% all_z_civ     = nan(num_quasars*7,1);
% all_N_civ     = nan(num_quasars*7,1);
% all_sigma_civ = nan(num_quasars*7,1);
% j=0;
% for this_quasar_ind=1:num_quasars
%     if mod(this_quasar_ind,100)==0
%         fprintf('QSO-%d\n',this_quasar_ind);
%     end
%     for num_c4=1:7
%         if p_c4(this_quasar_ind, num_c4)>0.85
%             j=j+1;
%             all_z_civ(j,1)=map_z_c4L2(this_quasar_ind,num_c4);
%             all_N_civ(j,1)=map_N_c4L2(this_quasar_ind,num_c4);
%             all_sigma_civ(j,1)=map_sigma_c4L2(this_quasar_ind,num_c4);
%         end
%     end
% end


% fig = figure();
% X = [all_z_civ,  all_N_civ];
% hist3(X,'CdataMode','auto', 'Nbins', [50,50])
% ylabel('$N_{CIV}$', 'interpreter', 'latex')
% xlabel('$Z_{CIV}$', 'interpreter', 'latex')
% colorbar
% view(2)
% exportgraphics(fig, 'Z-N-hist2D-50.png', 'resolution',800)
% close all

% fig = figure();
% X = [all_sigma_civ, all_N_civ];
% hist3(X,'CdataMode','auto', 'Nbins', [50,50])
% ylabel('$N_{CIV}$', 'interpreter', 'latex')
% xlabel('$\sigma_{CIV}$', 'interpreter', 'latex')
% colorbar
% view(2)
% exportgraphics(fig, 'sigma-N-hist2D-50.png', 'resolution',800)
