
%%%%%%%%-------------C13 Npts-----------%%%%

% clc
% clear
% set_parameters;
% build_catalog;
% variables_to_load= {'all_plate_dr7', 'all_mjd_dr7', 'all_fiber_dr7', ...
% 'all_QSO_ID', 'all_RA', 'all_DEC', 'all_zqso', 'EW1', 'EW2',...
% 'all_N_civ','all_z_civ', 'all_RATING', 'dla_QSO_ID','log_posteriors_dla',...
%  'log_posteriors_no_dla' , 'c4_QSO_ID'};
    
%  load(sprintf('%s/catalog', processed_directory(release)), ...
%  variables_to_load{:});
%  load(sprintf('%s/preloaded_qsos', processed_directory(release)));
% ind_has_c4 =ismember(all_QSO_ID, c4_QSO_ID);

% % checking how many Npts1 and Npts2 are there in C13's absorbers
% num_quasars = length(all_wavelengths);
% nSys =17; % total number of possible absorbers in a spectrum 
% N_pts1 = nan(num_quasars,nSys);
% N_pts2 = nan(num_quasars,nSys);
% for quasar_ind = 1:num_quasars    
%     z_qso            =          all_zqso(quasar_ind);
%     this_wavelengths = all_wavelengths{quasar_ind};
%     this_wavelengths =           this_wavelengths';
%     this_pixel_mask  =  all_pixel_mask{quasar_ind};
%     this_pixel_mask  =            this_pixel_mask';
%     this_sigma_pixel = all_sigma_pixel{quasar_ind};
%     this_sigma_pixel =           this_sigma_pixel';

%     % convert to QSO rest frame
%     this_rest_wavelengths = emitted_wavelengths(this_wavelengths, z_qso);

%     unmasked_ind = (this_rest_wavelengths >= min_lambda) & ...
%         (this_rest_wavelengths <= max_lambda) & (this_sigma_pixel>0);
%     % keep complete copy of equally spaced wavelengths for absorption
%     % computation
%     this_unmasked_wavelengths = this_wavelengths(unmasked_ind);

%     l1 = 1550.7810;
%     l2 = 1548.2040;
%     h=(l1-l2)*0.25;
%     for Sys=1:nSys
%         if(all_z_civ(quasar_ind,Sys)>0)
%             mask_pts1 = (this_unmasked_wavelengths>(1+all_z_civ(quasar_ind,Sys))*(l1-h)) & ...
%                         (this_unmasked_wavelengths<(1+all_z_civ(quasar_ind,Sys))*(l1+h));
%             mask_pts2 = (this_unmasked_wavelengths>(1+all_z_civ(quasar_ind,Sys))*(l2-h)) & ...
%                         (this_unmasked_wavelengths<(1+all_z_civ(quasar_ind,Sys))*(l2+h));

%             N_pts1(quasar_ind,Sys) = sum(mask_pts1);
%             N_pts2(quasar_ind,Sys) = sum(mask_pts2);
%         end
%     end
% end
% N_pts1 = reshape(N_pts1 ,num_quasars*nSys,1);
% N_pts2 = reshape(N_pts2 ,num_quasars*nSys,1);
% size(N_pts1)
% N_pts1 = N_pts1(isnan(N_pts1)==false);
% N_pts2 = N_pts2(isnan(N_pts2)==false);
% size(N_pts1)
% maxP = max(max(N_pts1), max(N_pts2));
% P=zeros(maxP,maxP);
% for i=0:maxP
%     for j=0:maxP
%         P(i+1,j+1) = sum(N_pts1==i & N_pts2==j);
%     end
% end
% P = P/sum(sum(P));



testing_set_name = 'N-pts-voigt-05A-p1p2'
filename = sprintf('data/dr7/processed/processed_qsos_R%s.mat', testing_set_name);   
% filename ='data/dr7/processed/processed_qsos_Rprior-fixed.mat';
load(filename, 'test_ind','p_no_c4', 'p_c4', 'map_z_c4', ...
                'map_N_c4', 'max_z_c4s', 'Dz','map_sigma_c4',...
                'N_pts1', 'N_pts2', 'log_posteriors_no_c4',  'log_posteriors_c4',...
                'log_posteriors_no_c4',  'log_posteriors_c4');



cooksey_catalog = ...
fitsread('data/dr7/distfiles/dr7qso_CIV_noBAL.fit', 'binarytable');
all_zqso          = cooksey_catalog{4};
ind_has_c4 =ismember(all_QSO_ID, c4_QSO_ID);
true_label = ind_has_c4(test_ind);
num_quasars = length(N_pts1);
tr = 0.90;

% ind_FP = p_c4>tr & ~ind_has_c4(test_ind);
% sum(ind_FP)

true_label = ~(all_RATING(test_ind,1)<2 & all_RATING(test_ind,2)<2 & ...
all_RATING(test_ind,3)<2 & all_RATING(test_ind,4)<2 & ...
all_RATING(test_ind,6)<2 & all_RATING(test_ind,5)<2 & ...
all_RATING(test_ind,7)<2 & all_RATING(test_ind,8)<2 & ...
all_RATING(test_ind,9)<2 & all_RATING(test_ind,10)<2 & ...
all_RATING(test_ind,11)<2 & all_RATING(test_ind,12)<2 &  ...
all_RATING(test_ind,14)<2 & all_RATING(test_ind,13)<2 & ...
all_RATING(test_ind,15)<2 & all_RATING(test_ind,16)<2 & ...
all_RATING(test_ind,17)<2);


% sum(ind_FP)

% auc_max = -100;
% prior_np = [0,0,0,0];
% y_true = ind_has_c4(test_ind);
[X,Y,T,AUC] =perfcurve(true_label, p_c4, 'true');
% auc_0 = AUC
% for j=1:maxP1
%     prior_np(j) = nnz(true_label(N_pts2==j-1))/nnz(true_label);
%     prior_no_np(j) = nnz(~true_label(N_pts2==j-1))/nnz(true_label);
% end
% % r = prior_np./prior_no_np; 
% % for a=[0, logspace(-5,0,10)]
% %     for i=1:num_quasars
% %         for j=1:4
% %             if(N_pts2(i)==j)
% %                 p_c4_test(i,1) = p_c4(i) - a*(1-r(j));

% %             end
% %         end
% %     end

% %     y_score = p_c4_test;
% %     y_true = ind_has_c4(test_ind);
% %     [X,Y,T,AUC] =perfcurve(y_true, y_score, 'true');

% %     if(AUC>auc_max)
% %         auc_max = AUC;
% %         a_best=a;
% %     endfor i=1:num_quasars
% %     for j=1:4
% %         if(N_pts2(i)==j)
% %             p_c4_test(i,1) = p_c4(i) - a*(1-r(j));

% %         end
% %     end






















auc_max = -100;
NN=20;
ii=0;
for a=[0, logspace(-10,0,NN)]
  ii=ii+1;
    for i=1:num_quasars
  
        if(N_pts2(i)>=maxP-1 & N_pts1(i)>=maxP-1)
            p_c4_test(i,1)=p_c4(i)*(1+a);
        else
            p_c4_test(i,1)=p_c4(i);
        end
    end
    y_score = p_c4_test;
    y_true = ind_has_c4(test_ind);
    [X,Y,T,AUC] =perfcurve(y_true, y_score, 'true');

    if(AUC>auc_max) 
        auc_max = AUC
        a_best = a;
    end
end
% % % for i=1:num_quasars
% % %     if(N_pts2(i)==0)
% % %         p_c4_test(i,1)=p_c4(i)/(1+ a_best);
% % %     end
% % %     if(N_pts2(i)==1)
% % %         p_c4_test(i,1)=p_c4(i)/(1+b_best);
% % %     end
% % %     if (N_pts2(i)==2)
% % %         p_c4_test(i,1)=p_c4(i)*(1+c_best);
% % %     end
% % %     if(N_pts2(i)==3);
% % %         p_c4_test(i,1)=p_c4(i)*(1+d_best);
% % %     end
% % % end
% a=1e-15
% for i=1:num_quasars
%     for j=1:4
%         if(N_pts2(i)<=1)
%             p_c4_test(i,1) = p_c4(i)-a;
%         end
%         if(N_pts2(i)>=2)
%             p_c4_test(i,1) = p_c4(i)+a;
%             if(p_c4_test(i,1)>=1)
%                 p_c4_test(i,1)=1;
%             end
%         end
%     end
% end
% ind_TP = (ind_has_c4(test_ind) & p_c4_test>tr);
% TP = nnz(ind_TP);
% TN = nnz(~ind_has_c4(test_ind) & p_c4_test<tr);
% FN = nnz(ind_has_c4(test_ind) & p_c4_test<tr);
% ind_FP = ~ind_has_c4(test_ind) & p_c4_test>tr;
% FP = nnz(ind_FP);
% P = nnz(ind_has_c4(test_ind) );
% N = nnz(~ind_has_c4(test_ind));
% confusion_matrix=[TP/P, FN/P; FP/N, TN/N]
% Accuracy = (TP+TN)/(P+N);
% ErrorRate = (FP+FN)/(P+N);
% % % fprintf('a_best:%e\n', a_best)

% % fig= figure();
% y_score = p_c4_test;
% y_true = ind_has_c4(test_ind);
% [X,Y,T,AUC] =perfcurve(y_true, y_score, 'true');
% AUC

% y_score = p_c4;
% y_true = ind_has_c4(test_ind);
% [X,Y,T,AUC0] =perfcurve(y_true, y_score, 'true');
% AUC0
% % plot(X,Y)
% % legend(sprintf('AUC=%.5f',AUC))
% % set(get(gca, 'YLabel'), 'String', 'TPR');
% % set(get(gca, 'XLabel'), 'String', 'FPR');
% % set(get(gca, 'Title'), 'String', sprintf('p:%.2f, FP:%d\nCM:[%.4f, %.4f; %.4f, %.4f]\nAccuracy:%.4f, Error Rate:%.4f\n',...
% % tr, FP, TP/P, FN/P, FP/N,TN/N,Accuracy,ErrorRate));
% % exportgraphics(fig, sprintf('Penalize-FP/ROC-Penalized.png'))
% p1=zeros(4);
% p2=zeros(4);
% p3=zeros(4);
% p4=zeros(4);
% ind_FP = (~ind_has_c4(test_ind) & p_c4>tr);
% ind_TP = (ind_has_c4(test_ind) & p_c4>tr);

% for i=0:3
%     for j=0:3
%         p1(i+1,j+1) = nnz(N_pts1==i & N_pts2==j);
%         p2(i+1,j+1) = nnz( ind_FP & N_pts1==i & N_pts2==j)/nnz(ind_FP);
%         p3(i+1,j+1) = nnz( ind_TP & N_pts1==i & N_pts2==j)/nnz(ind_TP);
%         p4(i+1,j+1) = nnz( ind_FP & N_pts1==i & N_pts2==j)./nnz(N_pts1==i & N_pts2==j);
        
%     end
% end
% p1
% p2

% p = zeros(4,1);
% for i=0:3
%     p(i+1) = nnz(N_pts1==i & ind_FP)/nnz(ind_FP);
% end
% p