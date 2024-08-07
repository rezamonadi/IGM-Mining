% Confusion matrix builder 

clear
set_parameters_dr7;
build_catalog_dr7;
% testing_set_name = 'IntPol';
testing_set_name = 'norm_2p';
% 
filename = sprintf('data/dr7/processed/processed_qsos_R%s.mat', testing_set_name);   
% filename ='data/dr7/processed/processed_qsos_Rprior-fixed.mat';
variables_to_load = { 'training_set_name', ...
    'c4_catalog_name', 'prior_ind', 'release', ...
    'test_ind', 'prior_z_qso_increase', ...
    'max_z_cut', 'min_z_c4s', 'max_z_c4s', ...
    'log_priors_no_c4', 'log_priors_c4', ...
    'log_likelihoods_no_c4',  ...
    'sample_log_likelihoods_c4L2', 'log_likelihoods_c4L2'...
    'log_posteriors_no_c4', 'log_posteriors_c4L1', 'log_posteriors_c4L2',...
    'model_posteriors', 'p_no_c4', ...
    'map_z_c4L2', 'map_N_c4L2', 'map_sigma_c4L2', 'p_c4'};
load(filename, variables_to_load{:});
cooksey_catalog = ...
fitsread('data/dr7/distfiles/dr7qso_CIV_noBAL.fit', 'binarytable');
load(sprintf('%s/preloaded_qsos', processed_directory(release)));
all_zqso          = cooksey_catalog{4};
ind_has_c4 =ismember(all_QSO_ID, c4_QSO_ID);

tr=0.85;
predicted_c4s =  p_c4>tr;
% predicted_labels={0};
% true_labels={0};
ind_has_c4_tst = ind_has_c4(test_ind);
% for i=1:nnz(test_ind)
%     if(predicted_c4s(i)==0)
%         predicted_label{i}='No CIV';
%     else
%         predicted_label{i}='CIV';
%     end

%     if(ind_has_c4_tst(i)==0)
%         true_labels{i}='CIV';
%     else
%         true_labels{i} = 'No CIV';
%     end
% end

figure();
confusionchart(ind_has_c4_tst, predicted_c4s)

% plotconfusion(ind_has_c4(test_ind), predicted_c4s)
% test_ind = test_ind & (DoubletRatio>1 | DoubletRatio <2);
ind_TP = (ind_has_c4(test_ind) & p_c4>tr);
TP = nnz(ind_TP);
TN = nnz(~ind_has_c4(test_ind) & p_c4<tr);
FN = nnz(ind_has_c4(test_ind) & p_c4<tr);
ind_FP = ~ind_has_c4(test_ind) & p_c4>tr;
FP = nnz(ind_FP);
P = nnz(ind_has_c4(test_ind) );
N = nnz(~ind_has_c4(test_ind));
confusion_matrix=[TP/P, FN/P; FP/N, TN/N];
Accuracy = (TP+TN)/(P+N);
ErrorRate = (FP+FN)/(P+N);
Sensitivity = TP/P;
Specificity = TN/N;
fprintf('p:%.2f\nFP:%d\nCM:[%.4f, %.4f; %.4f, %.4f]\nAccuracy:%.4f\nError Rate:%.4f\n'...
,tr, FP, TP/P, FN/P, FP/N,TN/N,Accuracy,ErrorRate);



fig= figure();
y_score = p_c4;
y_true = ind_has_c4(test_ind);
[X,Y,T,AUC] =perfcurve(y_true, y_score, 'true');
plot(X,Y)
legend(sprintf('AUC=%.5f',AUC))
set(get(gca, 'YLabel'), 'String', 'TPR');
set(get(gca, 'XLabel'), 'String', 'FPR');
set(get(gca, 'Title'), 'String', sprintf('%s\np:%.2f, FP:%d\nCM:[%.4f, %.4f; %.4f, %.4f]\nAccuracy:%.4f, Error Rate:%.4f\n',...
testing_set_name ,tr, FP, TP/P, FN/P, FP/N,TN/N,Accuracy,ErrorRate));
exportgraphics(fig, sprintf('ROC-%s.pdf', testing_set_name)...
                ,'ContentType','vector')

% % % % Checking if our FPs are intrinsic
% % % % meaning that (Z_qso - MAP(Z_civ))/(1+Z_qso)>3000km/s/c

% % % % ind_intrinsic = (all_zqso(test_ind) - map_z_c4L2 )./(all_zqso(test_ind)+1)< kms_to_z(3000);
% % % % ind_intrinsic = Dz< kms_to_z(3000);

% % % % all_noise_variance = all_noise_variance(test_ind);
% % % % all_flux = all_flux(test_ind);
% % % % nnz(ind_intrinsic)
% % % % test_noise_med = zeros(nnz(test_ind),1);
% % % % for i=1:nnz(test_ind)
% % % %     sn_med(i) = 1/median(sqrt(all_noise_variance{i})./abs(all_flux{i}));
% % % % end
% % % % median(sn_med(ind_FP))
% % % % median(sn_med)

% % % % % % test_RATING = all_RATING(test_ind);
% % % % % % test_N_CIV = all_N_civ(test_ind);
% % % % % % histogram(test_N_CIV(ind_FP & test_RATING>0,1))
% % % % % % histogram(test_RATING(ind_FP,1))

% % % % % % nnz(test_RATING(ind_FP,1)==-1) 
% % % % % % nnz(test_RATING(ind_FP,1)==0) 
% % % % % % nnz(test_RATING(ind_FP,1)==1)
% % % % % % nnz(test_RATING(ind_FP,1)==2)
% % % % % % nnz(test_RATING(ind_FP,1)==3)

% % % % % histogram(all_residual2)
% % % % % hold on 
% % % % % histogram(all_residual2(ind_FP))
% % % % % legend('all', 'FP')
% % % % % set(get(gca, 'YLabel'), 'String', 'N');
% % % % % set(get(gca, 'XLabel'), 'String', 'mean|flux - \mu|');
% % % % % % fa = 0.189900;
% % % % % % fb  = 0.094750;
% % % % % % DR = sqrt(log(1548*map_N_c4L2*fa./sqrt(2)./map_sigma_c4L2)./...
% % % % % %           log(1550*map_N_c4L2*fb./sqrt(2)./map_sigma_c4L2));

% % % % % figure();
% % % % % histogram(DoubletRatio) 
% % % % % hold on         
% % % % % histogram(DoubletRatio(DoubletRatio>sqrt(2) & DoubletRatio<2))          
% % % % % legend('all', 'physical')
% % % % % saveas(gca, 'DR.png')

% % % % % figure();
% % % % % histogram(DLikelihood) 

% % % % % saveas(gca, 'DL.png')
% % % % % % % % voigt_dr(   )
% % % % % % minRes = linspace(min(all_residual2), max(all_residual2), 10000);
% % % % % % FP_res=zeros(10000,1);
% % % % % % for i=1:100
% % % % % %     FP_res(i) = nnz(ind_FP & all_residual2>minRes(i));
% % % % % % end
% % % % % % plot(minRes, FP_res)

% % % % % figure();
% % % % % minDR = linspace(min(DoubletRatio), max(DoubletRatio), 100);
% % % % % FP_DR=zeros(100,1);
% % % % % for i=1:100
% % % % %     FP_DR(i) = nnz(ind_FP & DoubletRatio>minDR(i));
% % % % % end
% % % % % plot(minDR, FP_DR)
% % % % % saveas(gca, 'minFPDoubletRatio.png')

% % % % % figure();
% % % % % histogram(all_offsetl2)
% % % % % hold on
% % % % % histogram(all_offsetl2(ind_FP))
% % % % % set(get(gca, 'XLabel'), 'String', 'mean(flux-\mu)')
% % % % % set(get(gca, 'YLabel'), 'String', 'N');
% % % % % legend('all', 'FP')
% % % % % saveas(gca, 'FPOffset.png')

% % % % % figure();
% % % % % minOffset = linspace(min(all_offsetl2), max(all_offsetl2), 100);
% % % % % FP_offset=zeros(100,1);
% % % % % for i=1:100
% % % % %     FP_offset(i) = nnz(ind_FP & all_offsetl2<minOffset(i));
% % % % % end
% % % % % plot(minOffset, FP_offset, 'LineWidth',3)

% % % % % set(get(gca, 'XLabel'), 'String', 'min(mean(flux-\mu))')
% % % % % set(get(gca, 'YLabel'), 'String', 'N(FP)');
% % % % % saveas(gca, 'FPminOffset.png')


% % % % % figure();
% % % % % minEqW1 = linspace(min(EqW1), max(EqW1), 100);
% % % % % FP_EqW1=zeros(100,1);
% % % % % for i=1:100
% % % % %     FP_EqW1(i) = nnz(ind_FP==1 & EqW1<minEqW1(i));
% % % % % end
% % % % % plot(minEqW1, FP_EqW1, 'LineWidth',3)
% % % % % hold on 
% % % % % minEqW2 = linspace(min(EqW2), max(EqW2), 100);
% % % % % FP_EqW2=zeros(100,1);
% % % % % for i=1:100
% % % % %     FP_EqW2(i) = nnz(ind_FP==1 & EqW2<minEqW2(i));
% % % % % end
% % % % % plot(minEqW2, FP_EqW2, 'LineWidth',3)
% % % % % legend('EW1', 'EW2')
% % % % % save(gca, 'Eqw-FP.png')


% % % % fig= figure();
% % % % FP_b=zeros(100,1);
% % % % minSigma = linspace(5e5, 40e5, 100);

% % % % for i=1:100
% % % %     FP_b(i) = nnz(ind_FP==1 & map_sigma_c4L2<minSigma(i));
% % % % end
% % % % plot(minSigma, FP_b, 'LineWidth',3)
% % % % % legend('EW1', 'EW2')
% % % % exportgraphics(fig,'Sigma-FP-sigma15-30.png', 'Resolution',400)


% % % fig2= figure();
% % % h1 = histogram(map_sigma_c4(ind_FP), 10);
% % % v1 = h1.Values;
% % % edgs = h1.BinEdges;
% % % h2 = histogram(map_sigma_c4, edgs);
% % % v2 = h2.Values;
% % % ratio = v1./v2;
% % % centers = 0.5*(edgs(2:end) + edgs(1:end-1));
% % % clf();
% % % plot(centers, ratio, '-*')
% % % set(get(gca, 'XLabel'), 'String', '\sigma');
% % % set(get(gca, 'YLabel'), 'String', 'Ratio(FP)');
% % % exportgraphics(fig2,sprintf('FP-sigsma-%s.png',testing_set_name),...
% % %  'Resolution',400)    



% % % fig3= figure();
% % % h1 = histogram(map_N_c4L2(ind_FP), 15);
% % % v1 = h1.Values;
% % % edgs = h1.BinEdges;
% % % h2 = histogram(map_N_c4L2, edgs);
% % % v2 = h2.Values;
% % % ratio = v1./v2;
% % % centers = 0.5*(edgs(2:end) + edgs(1:end-1));
% % % clf();
% % % plot(centers, ratio, '-*')
% % % set(get(gca, 'XLabel'), 'String', 'N');
% % % set(get(gca, 'YLabel'), 'String', 'Ratio(FP)');
% % % exportgraphics(fig3,'FP-N.png', 'Resolution',400)    

% % % fw1_100 = fw1_100./(1+map_z_c4);
% % % fw2_100 = fw2_100./(1+map_z_c4);
% % % fw1_200 = fw1_200./(1+map_z_c4);
% % % fw2_200 = fw2_200./(1+map_z_c4);

% % % fig4= figure();
% % % h1 = histogram(fw1_100(ind_FP), 5);
% % % v1 = h1.Values;
% % % edgs = h1.BinEdges;
% % % h2 = histogram(fw1_100, edgs);
% % % v2 = h2.Values;
% % % ratio = v1./v2;
% % % centers = 0.5*(edgs(2:end) + edgs(1:end-1));
% % % clf();
% % % plot(centers, ratio, '-*')
% % % set(get(gca, 'XLabel'), 'String', 'fw1(100)');
% % % set(get(gca, 'YLabel'), 'String', 'Ratio(FP)');
% % % exportgraphics(fig4,'FP-fw1_100.png', 'Resolution',800)    


% % % fig5= figure();
% % % h1 = histogram(fw1_200(ind_FP), 5);
% % % v1 = h1.Values;
% % % edgs = h1.BinEdges;
% % % h2 = histogram(fw1_200, edgs);
% % % v2 = h2.Values;
% % % ratio = v1./v2;
% % % centers = 0.5*(edgs(2:end) + edgs(1:end-1));
% % % clf();
% % % plot(centers, ratio, '-*')
% % % set(get(gca, 'XLabel'), 'String', 'fw1(200)');
% % % set(get(gca, 'YLabel'), 'String', 'Ratio(FP)');
% % % exportgraphics(fig5,'FP-fw1_200.png', 'Resolution',800)    




% % % fig6= figure();
% % % h1 = histogram(fw2_100(ind_FP), 5);
% % % v1 = h1.Values;
% % % edgs = h1.BinEdges;
% % % h2 = histogram(fw2_100, edgs);
% % % v2 = h2.Values;
% % % ratio = v1./v2;
% % % centers = 0.5*(edgs(2:end) + edgs(1:end-1));
% % % clf();
% % % plot(centers, ratio, '-*')
% % % set(get(gca, 'XLabel'), 'String', 'fw2(100)');
% % % set(get(gca, 'YLabel'), 'String', 'Ratio(FP)');
% % % exportgraphics(fig6,'FP-fw2_100.png', 'Resolution',800)    

% % % fig7= figure();
% % % h1 = histogram(fw2_200(ind_FP), 5);
% % % v1 = h1.Values;
% % % edgs = h1.BinEdges;
% % % h2 = histogram(fw2_200, edgs);
% % % v2 = h2.Values;
% % % ratio = v1./v2;
% % % centers = 0.5*(edgs(2:end) + edgs(1:end-1));
% % % clf();
% % % plot(centers, ratio, '-*')
% % % set(get(gca, 'XLabel'), 'String', 'fw2(200)');
% % % set(get(gca, 'YLabel'), 'String', 'Ratio(FP)');
% % % exportgraphics(fig7,'FP-fw2_200.png', 'Resolution',800)    



% % % % fig8= figure();
% % % % h1 = histogram(min_ratio(ind_FP), 10);
% % % % v1 = h1.Values;
% % % % edgs = h1.BinEdges;
% % % % h2 = histogram(min_ratio, edgs);
% % % % v2 = h2.Values;
% % % % ratio = v1./v2;
% % % % centers = 0.5*(edgs(2:end) + edgs(1:end-1));
% % % % clf();
% % % % plot(centers, ratio, '-*')
% % % % set(get(gca, 'XLabel'), 'String', 'minRatio');
% % % % set(get(gca, 'YLabel'), 'String', 'Ratio(FP)');
% % % % exportgraphics(fig8,'FP-min_ratio.png', 'Resolution',400)    

% % % histogram(DLikelihood)
% % % histogram(all_offsetl2)
% fig9= figure();
% v = [nnz(N_pts2==0), nnz(N_pts2==1), nnz(N_pts2==2), nnz(N_pts2==3)];
% vFP = [nnz(N_pts2(ind_FP)==0), nnz(N_pts2(ind_FP)==1),...
%       nnz(N_pts2(ind_FP)==2), nnz(N_pts2(ind_FP)==3)];
% vTP = [nnz(N_pts2(ind_TP)==0), nnz(N_pts2(ind_TP)==1),...
%       nnz(N_pts2(ind_TP)==2), nnz(N_pts2(ind_TP)==3)];
% ratioFP = vFP./v;
% ratioTP = vTP./v;
% plot([0,1,2,3], ratioFP, '-*')
% hold on
% plot([0,1,2,3], ratioTP, '-*')
% legend('FP', 'TP');
% set(get(gca, 'XLabel'), 'String', '#pts');
% set(get(gca, 'YLabel'), 'String', 'FP/TP');
% set(get(gca, 'Title'), 'String', '1548A\pm0.5A');
% exportgraphics(fig9,'Num_pts2-FP-TP-05AP1P2.png', 'Resolution',400)  
% % % histogram(all_offsetl2)
% fig10= figure();
% v = [nnz(N_pts1==0), nnz(N_pts1==1), nnz(N_pts1==2), nnz(N_pts1==3)];
% vFP = [nnz(N_pts1(ind_FP)==0), nnz(N_pts1(ind_FP)==1),...
%       nnz(N_pts1(ind_FP)==2), nnz(N_pts1(ind_FP)==3)];
% vTP = [nnz(N_pts1(ind_TP)==0), nnz(N_pts1(ind_TP)==1),...
%       nnz(N_pts1(ind_TP)==2), nnz(N_pts1(ind_TP)==3)];
% ratioFP = vFP./v;
% ratioTP = vTP./v;
% plot([0,1,2,3], ratioFP, '-*')
% hold on
% plot([0,1,2,3], ratioTP, '-*')
% legend('FP', 'TP');
% set(get(gca, 'XLabel'), 'String', '#pts');
% set(get(gca, 'YLabel'), 'String', 'FP/TP');
% set(get(gca, 'Title'), 'String', '1550A\pm0.5A');
% exportgraphics(fig10,'Num_pts1-FP-TP-05AP1P2.png', 'Resolution',400)       

% % % histogram(DLikelihood)
% % % histogram(all_offsetl2)
