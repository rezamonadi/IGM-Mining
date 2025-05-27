% set_parameters_dr7; 
% build_catalog_dr7;
% filename ='data/dr7/processed/processed_qsos_tst_N-1250-1610-S-35-115-civWVL.mat';
% % filename ='data/dr7/processed/processed_qsos_tst_mask-1-prior-1-OccamRazor-1-nC4-30000-plt-0-MaskinP-0-fixedPr.mat';
% load(filename);
% 
% load('EW/REW_DR7_sigma_width_4.mat')

% % Saving Kathy's data      
% variables_to_load= {'all_plate_dr7', 'all_mjd_dr7', 'all_fiber_dr7', ...
% 'all_QSO_ID', 'all_RA', 'all_DEC', 'all_zqso',...
% 'all_N_civ','all_z_civ1', 'all_z_civ2', 'all_z_civ3', 'all_RATING', 'c4_QSO_ID', 'all_EW1', 'all_errEW1', 'all_errEW2' };
% load(sprintf('%s/catalog', processed_directory(release)), ...
%     variables_to_load{:});          
% 
% load('data/dr7/processed/preloaded_qsos_L-1310-1555-mnp-400-normL-1420-1475-dlambda-0.50-k-20-mnv-0.25.mat')

load('data/dr7/processed/processed_qsos_tst_1_Masking_CIV-20A1%-Masking-20-CIV-samp-20k.mat');
num_quasars = nnz(test_ind);
all_wavelengths    =    all_wavelengths;   %G(test_ind);
all_flux           =           all_flux;   %G(test_ind);
% all_noise_variance = all_noise_variance(test_ind);
% all_pixel_mask     =     all_pixel_mask(test_ind);
% all_sigma_pixel    =    all_sigma_pixel(test_ind);
z_qsos             =      all_zqso(test_ind);
EW_C13_2796             =       all_EW1(test_ind,:);
errEW_C13_2796         =       all_errEW1(test_ind,:);
EW_C13_2803             =       all_EW2(test_ind,:);
errEW_C13_2803         =       all_errEW2(test_ind,:);

% % comparing Z and N
Z_C13_1 = all_z_MgII1(test_ind,:);
Z_C13_2 = all_z_MgII2(test_ind,:);
Z_C13_3 = all_z_MgII3(test_ind,:);
N_C13 = all_N_MgII(test_ind,:);
Rate_test = all_RATING(test_ind, :);
testID = all_QSO_ID(test_ind);
dv=10000000000000000000000000;
ind_low = zeros([num_quasars,1]);
ind_low_numMgII = zeros([num_quasars,1]);
ind_high = zeros([num_quasars,1]);
ind_high_numMgII = zeros([num_quasars,1]);
ind_0 = zeros([num_quasars,1]);
ind_0_numMgII = zeros([num_quasars,1]);

ind_dv_low = zeros([num_quasars,1]);
ind_dv_low_numMgII = zeros([num_quasars,1]);
ind_dv_high = zeros([num_quasars,1]);
ind_dv_high_numMgII = zeros([num_quasars,1]);
ind_dv_mid = zeros([num_quasars,1]);
ind_dv_mid_numMgII = zeros([num_quasars,1]);
ind_EW_large_PM1GP0 = nan([num_quasars, 1]);
DV_EW_large_PM1GP0 = nan([num_quasars, 1]);
ind_PM1GP0 = nan([num_quasars, 1]);
ii=0;
PM1GP0_QSO_ID = [];
for tr=0.85
    ii=ii+1;
    PM=0;
    GP=0;
    PM1GP1=0; 
    PM0GP1 = 0;
    PM1GP0 = 0;
    PMn1GP1 = 0; 
    PM1GP1dv0=0;
    PM1GP1highREW = 0; 
    lw =5;
    dNdZ_PM0GP1 = nan([num_quasars, max_MgII]);
    N_MgII_PM1GP1 = nan([num_quasars, max_MgII]);
    s_MgII_PM1GP1 = nan([num_quasars, max_MgII]);
    N_MgII_PM1GP1 = nan([num_quasars, max_MgII]);
    s_MgII_PM0GP1 = nan([num_quasars, max_MgII]);
    N_MgII_GP1 = nan([num_quasars, max_MgII]);
    s_MgII_GP1 = nan([num_quasars, max_MgII]);
    REW_2796_dr7_PM1GP1_PM = nan([num_quasars, max_MgII]);
    REW_2796_dr7_PM1GP1_PMerr = nan([num_quasars, max_MgII]);
    REW_2796_dr7_PM1GP1_GPflux = nan([num_quasars, max_MgII]);
    REW_2796_dr7_PM1GP1_GPfluxErr = nan([num_quasars, max_MgII]);
    REW_2796_dr7_PM1GP1_GPVoigt = nan([num_quasars, max_MgII]);
    REW_2803_dr7_PM1GP1_PM = nan([num_quasars, max_MgII]);
    REW_2803_dr7_PM1GP1_PMerr = nan([num_quasars, max_MgII]);
    REW_2803_dr7_PM1GP1_GPflux = nan([num_quasars, max_MgII]);
    REW_2803_dr7_PM1GP1_GPfluxErr = nan([num_quasars, max_MgII]);
    REW_2803_dr7_PM1GP1_GPVoigt = nan([num_quasars, max_MgII]);
    map_sigma_MgIIL2_PM1GP1 = nan([num_quasars,max_MgII]);
    map_N_MgIIL2_PM1GP1 = nan([num_quasars,max_MgII]);
    REW_2796_dr7_PM0GP1 = nan([num_quasars, max_MgII]);
    REW_2796_dr7_PM1GP0 = nan([num_quasars, 17]);
    DN = nan(num_quasars,max_MgII);
    rdz = nan(num_quasars,max_MgII);
    Dv_PM1GP1 = nan([num_quasars, max_MgII]);
    indPM1GP1_in_GP = nan(num_quasars,max_MgII);
    indPM0GP1_in_GP = nan(num_quasars,max_MgII);
    indPM1GP0_in_GP = nan(num_quasars,max_MgII);
    diffEW2796 = nan([num_quasars,max_MgII]);
    diffEW2803 = nan([num_quasars,max_MgII]);
    zMgII_PM1GP1_PM = nan([num_quasars, max_MgII]);
    nn=0;mm=0;
    clear highTPID indPM0GP1
 
    PM_ind = [];
    GP_ind = [];
    for quasar_ind=1:num_quasars
        for i=1:17
            if(Rate_test(quasar_ind,i)>=2)
                PM = PM+1;
                PM_ind(end+1) = quasar_ind;
            end
        end

       

       
        for j=1:max_MgII
            if (p_MgII(quasar_ind, j)>=tr)
                GP = GP+1;
                N_MgII_GP1(quasar_ind, j) = map_N_MgIIL2(quasar_ind, j);
%                 s_MgII_GP1(quasar_ind, j) = map_sigma_MgIIL2(quasar_ind, j);
                GP_ind(end+1) = quasar_ind;
            end
        end

        
       
        thisTP=0;
        for i=1:17
            % if(Rate_test(quasar_ind,i)>=2)
            if(Z_C13_1(quasar_ind, i)>0)
                

                DZ = Z_C13_1(quasar_ind, i) - map_z_MgIIL2(quasar_ind, :); %Z_seffert
                [DZ_min_abs, indMin] = min(abs(DZ));
                DZ_min = DZ(indMin);
                
                
                    j = indMin;
                    if (p_MgII(quasar_ind, j)>=tr)  % PM1GP1
                        %if DZ_min_abs<kms_to_z(dv)*(1+Z_C13_1(quasar_ind, i))
                            thisTP = thisTP+1;
                            % DZ_min = DZ_abs;
                            PM1GP1 = PM1GP1 + 1;
                            thisTP=thisTP+1;
                            N_MgII_PM1GP1(quasar_ind, j) = map_N_MgIIL2(quasar_ind, j);
    %                         s_c4_PM1GP1(quasar_ind, j) = map_sigma_c4L2(quasar_ind, j);
                            % 2796
%G                            REW_2796_dr7_PM1GP1_GPVoigt(quasar_ind, j) = REW_2796_DR7_voigt(quasar_ind, j);
%G                            REW_2796_dr7_PM1GP1_GPflux(quasar_ind, j) = REW_2796_DR7_flux(quasar_ind, j);
                            REW_2796_dr7_PM1GP1_PM(quasar_ind, j) = EW_C13_2796(quasar_ind, i);
                            REW_2796_dr7_PM1GP1_PMerr(quasar_ind, j) = errEW_C13_2796(quasar_ind,i);
  %g                          REW_2796_dr7_PM1GP1_GPfluxErr(quasar_ind, j) = ErrREW_2796_flux(quasar_ind,j);
                            map_sigma_MgIIL2_PM1GP1(quasar_ind,j) = map_sigma_MgIIL2(quasar_ind,j);
                            map_N_MgIIL2_PM1GP1(quasar_ind, j) = map_N_MgIIL2(quasar_ind, j);
                            zMgII_PM1GP1_PM(quasar_ind, j) = Z_C13_1(quasar_ind, i);
                            Dv_PM1GP1(quasar_ind, j) = DZ_min/(1+Z_C13_1(quasar_ind,i))*speed_of_light/1000;
                            indPM1GP1_in_GP(quasar_ind,j) = 1;
                            indPM1GP1(PM1GP1,2) = j;
                            DN(quasar_ind, j) = map_N_MgIIL2(quasar_ind, j) - all_N_MgII(quasar_ind,i);

                            errMax = max(REW_2796_dr7_PM1GP1_PMerr(quasar_ind, j), REW_2796_dr7_PM1GP1_GPfluxErr(quasar_ind, j));
                            diffEW2796(quasar_ind,j) = (REW_2796_dr7_PM1GP1_GPflux(quasar_ind, j) - REW_2796_dr7_PM1GP1_PM(quasar_ind, j))/errMax;
                            if diffEW2796(quasar_ind,j)>2
                                    ind_high(quasar_ind) = 1;
                                ind_high_numMgII(quasar_ind)=j;
                            end
                            if diffEW2796(quasar_ind,j)<-2
                                ind_low(quasar_ind) = 1;
                                ind_low_numMgII(quasar_ind)=j;
                            end
                            if abs(diffEW2796(quasar_ind,j))<0.02
                                ind_0(quasar_ind) = 1;
                                ind_0_numMgII(quasar_ind) = j;
                            end

                            % 2803
                %G            REW_2803_dr7_PM1GP1_GPVoigt(quasar_ind, j) = REW_2803_DR7_voigt(quasar_ind, j);
                 %G           REW_2803_dr7_PM1GP1_GPflux(quasar_ind, j) = REW_2803_DR7_flux(quasar_ind, j);
                            REW_2803_dr7_PM1GP1_PM(quasar_ind, j) = EW_C13_2803(quasar_ind, i);
                            REW_2803_dr7_PM1GP1_PMerr(quasar_ind, j) = errEW_C13_2803(quasar_ind,i);
                 %G           REW_2803_dr7_PM1GP1_GPfluxErr(quasar_ind, j) = ErrREW_2803_flux(quasar_ind,j);
                            errTotal = sqrt(REW_2803_dr7_PM1GP1_PMerr(quasar_ind, j)^2 +  REW_2803_dr7_PM1GP1_GPfluxErr(quasar_ind, j)^2);
                            diffEW2803(quasar_ind,j) = (REW_2803_dr7_PM1GP1_GPflux(quasar_ind, j) - REW_2803_dr7_PM1GP1_PM(quasar_ind, j))/errTotal;
                            


                            if (Dv_PM1GP1(quasar_ind, j)>50)
                                
                                ind_dv_high(quasar_ind) = 1;
                                ind_dv_high_numMgII(quasar_ind) = j;
                            end
                            if (Dv_PM1GP1(quasar_ind, j)<-150)
                                ind_dv_low(quasar_ind) =1;
                                ind_dv_low_numMgII(quasar_ind)=j;
                            end
                            if (Dv_PM1GP1(quasar_ind, j)>-50) & (Dv_PM1GP1(quasar_ind, j)<-48)
                                ind_dv_mid(quasar_ind) =1;
                                ind_dv_mid_numMgII(quasar_ind)=j;
                            end

                            if REW_2796_dr7_PM1GP1_GPflux(quasar_ind, i)> 1.2
                                PM1GP1highREW = PM1GP1highREW+1;
                            end
                            

    %                     else % --> this is a missed one by GP --> PM1GP0
    %                         PM1GP0 = PM1GP0 +1;
    %                         ind_PM1GP0(quasar_ind) = 1;
    %                         REW_2796_dr7_PM1GP0(quasar_ind, i) = EW_C13_2796(quasar_ind, i);
    %                         if REW_2796_dr7_PM1GP0(quasar_ind, i)> 1.2
    %                             ind_EW_large_PM1GP0(quasar_ind) =  1;
    %                             ind_EW_large_PM1GP0_numc4(quasar_ind) =  indMin;
    %                             DV_EW_large_PM1GP0(quasar_ind) =DZ_min/(1+z_qsos(quasar_ind))*speed_of_light/1000;
    %                            
    %                            
    %                         end
                        %else
                                PM1GP1dv0 =PM1GP1dv0+1;

                        %end
                    else
                        if (p_MgII(quasar_ind,i)>=0) % try >= 0 
                        PM1GP0 = PM1GP0 +1;
                       PM1GP0_QSO_ID(end+1) = quasar_ind;
                        PM1GP0_QSO_ID = testID(PM1GP0_QSO_ID) %%%%

                        %G   endmax_MgII
                        end
                     end
                   
            end

        end
        
        if thisTP>=5
            nn=nn+1;
            highTPID(nn) = testID(quasar_ind);
        end


        
        
        for j=1:max_MgII
            if (p_MgII(quasar_ind, j)>=tr) & (all(Rate_test(quasar_ind, :)<2))
                    PM0GP1 = PM0GP1 + 1;
                    dNdZ_PM0GP1(quasar_ind, j) =1/(max_z_MgIIs(quasar_ind) - min_z_MgIIs(quasar_ind));
                    N_MgII_PM0GP1(quasar_ind, j) = map_N_MgIIL2(quasar_ind, j);
%                     s_MgII_PM0GP1(quasar_ind, j) = map_sigma_MgIIL2(quasar_ind, j);
                    %g REW_2796_dr7_PM0GP1(quasar_ind, j)=REW_2796_DR7_flux(quasar_ind, j);  
                    mm = mm+1;
                    indPM0GP1(mm) = quasar_ind;
                    ind_PM0GP1_in_GP(quasar_ind, j) =1;
            end
        end

        
    %    

    end

    

    P(ii) = PM1GP1/GP 
    C(ii) = PM1GP1/PM
end
 PM_QSO_ID = testID(PM_ind);

commonQSO_IDs = intersect(PM_QSO_ID, GP_QSO_ID);
uncommonQSO_IDs = setxor(PM_QSO_ID, GP_QSO_ID);
size(commonQSO_IDs)
size(uncommonQSO_IDs)
% for ii2 = 1:length(testID)
%     for jj2 = 1:length(commonQSO_IDs)
%         if testID(ii2) == commonQSO_IDs(jj2)


% %% velocity  comparison for different zciv_PM
% fig = figure();
% x = reshape(Dv_PM_GP, [], 1);
% x = x(~isnan(x));
% p=histogram(x);
% p.NumBins = 50;
% p.Normalization = 'count';
% hold on 
% ax = gca;
% ax.FontSize = 14;
% ylabel('Count')

% exportgraphics(fig, 'dvPMGP.png', 'resolution',800);


% %%   scatter dv_{PM,GP} vs.. Z_PM
% fig = figure();
% y = reshape(Dv_PM1GP1, [], 1);
% y = y(~isnan(y));
% x = reshape(zciv_PM1GP1_PM, [], 1);
% x = x(~isnan(x));
% p=scatter(x,y, 10, 'filled')
% p.MarkerFaceAlpha = 0.7;
% hold on 
% xLine = linspace(min(x), max(x), 1000);
% p = plot(xLine, zeros([1000,1])+150, 'LineWidth',1, 'LineStyle','--')
% p.Color = [1,0,0,0.6];
% hold on 
% xLine = linspace(min(x), max(x), 1000);
% p = plot(xLine, zeros([1000,1])-150, 'LineWidth',1, 'LineStyle','--')
% p.Color = [1,0,0,0.6];
% hold on 
% xLine = linspace(min(x), max(x), 1000);
% p = plot(xLine, zeros([1000,1]), 'LineWidth',1.5, 'LineStyle','-')
% p.Color = [1,0,0,0.6];
% ylabel('$\delta v_{\rm PM,GP}~({\rm km.s^{-1}})$', 'interpreter', 'latex')
% xlabel('$z_{CIV}^{\rm PM}$', 'interpreter', 'latex')
% xlim([min(x)-0.05, max(x)+0.01])
% ylim([min(y)-5, max(y)+5])
% set(gca, 'fontsize', 12)
% exportgraphics(fig, 'dcatterDV_PM1GP1_ZPM.png', 'Resolution', 800)



% % Purity and completeness 
% fig = figure()
% % histogram(reshape(DN, length(DN)*7, 1), 15)
% % set(get(gca, 'XLabel'), 'String', 'MAP(N)-C13(N)');
% % set(get(gca, 'YLabel'), 'String', 'Frequency');
% % set(gca, 'FontSize', 15)

% % exportgraphics(fig, 'DN.pdf', 'ContentType', 'vector')
% trx = 0:0.001:1;
% p = plot(trx, C, 'LineWidth', 1.5, 'LineStyle','-');
% p.Color = [0.8, 0.2, 0.2, 0.5];
% hold on 
% p=plot(trx, P, 'LineWidth', 1.5, 'LineStyle','--');
% p.Color = [0.2, 0.2, 0.8, 0.5];

% hold on
% xlim([0, 1.02])
% legend({'Completeness', 'Purity'},'Location','southeast')
% xlabel('$P_n(\mathcal{M_D})$','Interpreter', 'Latex')
% ylabel('Purity/Completeness')
% exportgraphics(fig, 'PC-N-1250-1610-S-35-115-NoOcc-nc-10k.png', 'Resolution', 800)
 

% % REW of different categories 

% h = histogram(reshape(REW_2796_dr7_PM1GP1_PM(~isnan(REW_2796_dr7_PM1GP1_PM)),[],1), 'normalization', 'pdf');
% h.BinWidth =0.25;
% x3 = h.BinEdges;
% x3 = x3(2:end) - h.BinWidth/2;
% y3 = h.Values;


% fig = figure();
% clf();

% h = histogram(reshape(REW_2796_dr7_PM1GP0(~isnan(REW_2796_dr7_PM1GP0)),[],1), 'normalization', 'pdf');
% h.BinWidth =0.25;
% x1 = h.BinEdges;
% x1 = x1(2:end) - h.BinWidth/2;
% y1 = h.Values;
% h.FaceAlpha = 0.7;
% hold on

% h = histogram(reshape(REW_2796_dr7_PM0GP1(~isnan(REW_2796_dr7_PM0GP1) & REW_2796_dr7_PM0GP1>0),[],1), 'normalization', 'pdf');
% h.BinWidth=0.25;
% x2 = h.BinEdges;
% x2 = x2(2:end) - h.BinWidth/2;
% y2 = h.Values;
% h.FaceAlpha = 0.7;
% hold on 

% p = plot(x3,y3);
% p.LineWidth=3;
% p.Color = [0.2, 0.2, 0.2, 0.7]
% hold on 
% legend({'PM Only', 'GP Only', 'PM $\&$ GP'}, 'interpreter', 'latex')
% set(gca, 'fontsize', 15)
% ylabel('PDF')
% xlabel('${\rm W}_{r,2796}$', 'interpreter', 'latex')
% exportgraphics(fig, 'REW_PM0GP1_PM1GP0_PM1GP1.png', 'resolution', 800)


% % Scatter plot of W_r,GP,voigt vs. W_r,PM
% fig = figure()
% y = reshape(REW_2796_dr7_PM1GP1_GPVoigt(:,:), [], 1); % W_r,GP,voigt
% x = reshape(REW_2796_dr7_PM1GP1_PM(:,:), [], 1); % W_r,PM 
% e = reshape(REW_2796_dr7_PM1GP1_PMerr(:,:,3), [], 1); % err from PM
% s = scatter(x,y, 10, 'filled')
% s.MarkerFaceAlpha = 0.6;
% hold on 
% xLine = linspace(min(x), max(x), 1000);
% p = plot(xLine, xLine, 'LineWidth',1, 'LineStyle','-')
% p.Color = [1,0,0,0.6];
% xlabel('${{\rm W}_r}^{\rm PM}$~(\AA)', 'interpreter', 'latex')
% ylabel('${{\rm W}_r}^{\rm GP}$~(\AA)', 'interpreter', 'latex')
% xlim([min(x)-0.1, max(x)+.1])
% ylim([min(y)-0.1, max(y)+.1])

% exportgraphics(fig, 'REW_scatter.png', 'Resolution',800)




% % Scatter plot of W_r,GP vs. W_r,PM
% fig = figure()
% dW_r  = REW_2796_dr7_PM1GP1(:,:,1) - REW_2796_dr7_PM1GP1(:,:,2); % W_r,GP - W_r,PM
% y = reshape(dW_r, [], 1); % dW_r,PM,GP
% x = reshape(Dv_PM_GP3, [], 1);

% s = scatter(x,y, 10, 'filled')
% s.MarkerFaceAlpha =0.6;
% ylabel('${{\rm W}_r}^{\rm GP} - {{\rm W}_r}^{\rm PM}$~(\AA)', 'interpreter', 'latex')
% xlabel('$\delta v_{\rm GP,PM}~(km.s^{-1})$', 'interpreter', 'latex')
% % xlim([min(x)-0.1, max(x)+.1])
% % ylim([min(y)-0.1, max(y)+.1])

% exportgraphics(fig, 'v-dW_r_scatter.png', 'Resolution',800)


% % scatter of dW_r/err vs. W_r,PM
% fig =  figure();
% y = dW_r./REW_2796_dr7_PM1GP1(:,:,3);
% y = reshape(y,[],1);
% x = REW_2796_dr7_PM1GP1(:,:,2);
% x = reshape(x,[],1);

% s = scatter(x,y, 10, 'filled');
% s.MarkerFaceAlpha =0.7;
% hold on
% xLine = linspace(min(x), max(x), 1000);
% yLine = zeros(size(xLine));
% p=plot(xLine, yLine);
% p.Color = [1,0,0,0.7];
% hold on 
% yLine = zeros(size(xLine)) +2;

% p=plot(xLine, yLine, 'lineStyle','--');
% p.Color = [1,0,0,0.7];
% hold on 
% yLine = zeros(size(xLine)) -2;
% p=plot(xLine, yLine, 'lineStyle','--');
% p.Color = [1,0,0,0.7];
% hold on 

% xlim([min(x)-0.05, max(x)+.05])
% ylim([min(y)-0.05, max(y)+.05])
% xlabel('${{\rm W}_r}^{\rm PM}$~(\AA)', 'interpreter', 'latex')
% ylabel('${({\rm W}_r}^{\rm GP} - {{\rm W}_r}^{\rm PM})/err({\rm W}_r^{\rm PM})$', 'interpreter', 'latex')
% exportgraphics(fig, 'dW_r-err_scatter.png', 'Resolution',800)




% % Scatter plot of W_r,GP vs. W_r,PM
% fig = figure()
% dW_r_flux  = REW_2796_dr7_PM1GP1(:,:,4) - REW_2796_dr7_PM1GP1(:,:,2); % W_r,GP - W_r,PM
% y = reshape(dW_r, [], 1); % dW_r,PM,GP
% x = reshape(Dv_PM_GP3, [], 1);

% s = scatter(x,y, 10, 'filled')
% s.MarkerFaceAlpha =0.6;
% ylabel('${{\rm W}_r}^{\rm GP} - {{\rm W}_r}^{\rm PM}$~(\AA)', 'interpreter', 'latex')
% xlabel('$\delta v_{\rm GP,PM}~(km.s^{-1})$', 'interpreter', 'latex')
% % xlim([min(x)-0.1, max(x)+.1])
% % ylim([min(y)-0.1, max(y)+.1])

% exportgraphics(fig, 'v-dW_r_flux_scatter.png', 'Resolution',800)


% % scatter of dW_r/err vs. W_r,PM
% fig =  figure();
% y = dW_r_flux./REW_2796_dr7_PM1GP1(:,:,3);
% y = reshape(y,[],1);
% x = REW_2796_dr7_PM1GP1(:,:,2);
% x = reshape(x,[],1);

% s = scatter(x,y, 10, 'filled');
% s.MarkerFaceAlpha =0.7;
% hold on
% xLine = linspace(min(x), max(x), 1000);
% yLine = zeros(size(xLine));
% p=plot(xLine, yLine);
% p.Color = [1,0,0,0.7];
% hold on 
% yLine = zeros(size(xLine)) +2;

% p=plot(xLine, yLine, 'lineStyle','--');
% p.Color = [1,0,0,0.7];
% hold on 
% yLine = zeros(size(xLine)) -2;
% p=plot(xLine, yLine, 'lineStyle','--');
% p.Color = [1,0,0,0.7];
% hold on 

% xlim([min(x)-0.05, max(x)+.05])
% ylim([min(y)-0.05, max(y)+.05])
% xlabel('${{\rm W}_r}^{\m PM}$~(\AA)', 'interpreter', 'latex')
% ylabel('${({\rm W}_r}^{\m GP} - {{\rm W}_r}^{\rm PM})/err({\rm W}_r^{\rm PM})$', 'interpreter', 'latex')
% exportgraphics(fig, 'dW_r-err_flux_scatter.png', 'Resolution',800)


% %  Boxcar GP vs PM
% fig = figure()
% x = reshape(REW_2796_dr7_PM1GP1_GPflux(:,:), [], 1); % W_r,GP, flux 
% xe = reshape(REW_2796_dr7_PM1GP1_GPfluxErr(:,:), [], 1); % W_r,GP, flux error
% y = reshape(REW_2796_dr7_PM1GP1_PM(:,:), [], 1); % W_r,PM
% ye = reshape(REW_2796_dr7_PM1GP1_PMerr(:,:),[],1); % error(W_r,PM)
% s = scatter(x,y, 10, 'filled')
% s.MarkerFaceAlpha =0.6;
% % x_median = nanmedian(x)
% % x=x-x_median;
% % y = y - nanmedian(y);
% pp = errorbar(x,y,ye, 'LineStyle','none','Marker','.', 'MarkerSize',6 )

% hold on 
% xLine = linspace(min(x), max(x), 1000);
% p = plot(xLine, xLine, 'LineWidth',1, 'LineStyle','-')
% p.Color = [1,0,0,0.6];
% ylabel('${{\rm W}_r}^{\rm PM}$~(\AA)', 'interpreter', 'latex')
% xlabel('${{\rm W}_r}^{\rm GP,flux}$~(\AA)', 'interpreter', 'latex')
% xlim([min(x)-0.1, max(x)+.1])
% ylim([min(y)-0.1, max(y)+.1])
% hold on
% [a_coef, b_coef, R, a, b]= errFitter(x(~isnan(x)), y(~isnan(y)),...
%                                     zeros([length(~isnan(xe)),1]),...
%                                     ye(~isnan(ye)), 10000);
% ttl =  sprintf('%.3f_{%.3f}^{%.3f}x + (%.3f)_{%.3f}^{%.3f} <R>=%.3f',a, b, median(R));                                     
% title(ttl)
% hold on 
% plot(x, a(1)*x + b(1))

% exportgraphics(fig, 'EWGPfluxPM-sw4.png', 'Resolution',800)

% %  Boxcar GP vs PM both error
% fig = figure()
% x = reshape(REW_2796_dr7_PM1GP1_GPflux(:,:), [], 1); % W_r,GP, flux 
% xe = reshape(REW_2796_dr7_PM1GP1_GPfluxErr(:,:), [], 1); % W_r,GP, flux error
% y = reshape(REW_2796_dr7_PM1GP1_PM(:,:), [], 1); % W_r,PM
% ye = reshape(REW_2796_dr7_PM1GP1_PMerr(:,:),[],1); % error(W_r,PM)
% s = scatter(x,y, 10, 'filled')
% s.MarkerFaceAlpha =0.6;
% pp = errorbar(x,y,ye,ye,xe,xe, 'LineStyle','none','Marker','.', 'MarkerSize',6 )
% % x = x- nanmean;
% % y = y-nanmean(y);
% hold on 
% xLine = linspace(min(x), max(x), 1000);
% p = plot(xLine, xLine, 'LineWidth',1, 'LineStyle','-')
% p.Color = [1,0,0,0.6];
% ylabel('${{\rm W}_r}^{\rm PM}$~(\AA)', 'interpreter', 'latex')
% xlabel('${{\rm W}_r}^{\rm GP,flux}$~(\AA)', 'interpreter', 'latex')
% xlim([min(x)-0.1, max(x)+.1])
% ylim([min(y)-0.1, max(y)+.1])
% hold on
% [a_coef, b_coef, R, a, b]= errFitter(x(~isnan(x)), y(~isnan(y)),...
%                                     xe(~isnan(xe)),...
%                                     ye(~isnan(ye)), 10000);
% ttl =  sprintf('%.3f_{%.3f}^{%.3f}x + (%.3f)_{%.3f}^{%.3f}, <R>=%.3f',a, b, median(R));                                     
% title(ttl)
% hold on 
% plot(x, a(1)*x + b(1))

% exportgraphics(fig, 'EWGPfluxPM-sw4-both-e.png', 'Resolution',800)


% % Voigt vs Boxcar GP
% fig = figure()
% x = reshape(REW_2796_dr7_PM1GP1_GPVoigt(:,:), [], 1); % W_r,GP, voigt
% y = reshape(REW_2796_dr7_PM1GP1_GPflux(:,:), [], 1); % W_r,GP, flux
% ye = reshape(REW_2796_dr7_PM1GP1_GPfluxErr(:,:),[],1); % error(W_r,GP,flux)
% % s = scatter(x,y, 10, 'filled')
% % s.MarkerFaceAlpha =0.6;

% errorbar(x,y,ye, 'LineStyle','none','Marker','.', 'MarkerSize',10)
% hold on 
% xLine = linspace(min(x), max(x), 1000);
% p = plot(xLine, xLine, 'LineWidth',1, 'LineStyle','-')
% p.Color = [1,0,0,0.6];
% ylabel('${{\rm W}_r}^{\rm GP,flux}$~(\AA)', 'interpreter', 'latex')
% xlabel('${{\rm W}_r}^{\rm GP,Voigt}$~(\AA)', 'interpreter', 'latex')
% xlim([min(x)-0.1, max(x)+.1])
% ylim([min(y)-0.1, max(y)+.1])

% exportgraphics(fig, 'EWGPVoigtFlux-sw4.png', 'Resolution',800)




% % EW Voigt GP vs PM
% fig = figure()
% x = reshape(REW_2796_dr7_PM1GP1_GPVoigt, [], 1); % W_r,GP, voigt
% y = reshape(REW_2796_dr7_PM1GP1_PM, [], 1); % W_r,PM
% ye = reshape(REW_2796_dr7_PM1GP1_PMerr,[],1); % error(W_r,PM)
% errorbar(x,y, ye, 'LineStyle','none','Marker','.', 'MarkerSize',10)
% hold on 
% xLine = linspace(min(x), max(x), 1000);
% p = plot(xLine, xLine, 'LineWidth',1, 'LineStyle','-')
% p.Color = [1,0,0,0.6];
% xlabel('${{\rm W}_r}^{\rm GP, Voigt}$~(\AA)', 'interpreter', 'latex')
% ylabel('${{\rm W}_r}^{\rm PM}$~(\AA)', 'interpreter', 'latex')
% xlim([min(x)-0.1, max(x)+.1])
% ylim([min(y)-0.1, max(y)+.1])

% [a_coef, b_coef, R, a, b]= errFitter(x(~isnan(x)), y(~isnan(y)),...
%                                     zeros([length(x(isnan(x))),1]),...
%                                     ye(~isnan(ye)), 100);
% title(sprintf('%.3f_{%.3f}^{%.3f}x + (%.3f)_{%.3f}^{%.3f}, <R>:%.3f',...
%               a, b, median(R)))
% hold on
% plot(x, a(2)*x + b(2))
% exportgraphics(fig, 'EWGPVoigtPM-sw4.png', 'Resolution',800)



% % scatter (EW(GP,Voigt) - EW(PM))/total(err) vs EW(PM)
% fig = figure()
% y = reshape(REW_2796_dr7_PM1GP1_GPVoigt(:,:), [], 1); % W_r,GP, voigt
% x = reshape(REW_2796_dr7_PM1GP1_PM(:,:), [], 1); % W_r,PM
% difEW2796 = y-x;
% ye = reshape(REW_2796_dr7_PM1GP1_PMerr(:,:),[],1); % error(W_r,PM)
% xe = reshape(REW_2796_dr7_PM1GP1_GPfluxErr(:,:), [], 1); %error() W_r,GP, flux)
% errTotal = sqrt(xe.^2 + ye.^2);
% yPlot = difEW2796./errTotal; % W_r,PM - W_r,GP,flux / errTotal
% c = reshape(map_sigma_c4L2_PM1GP1,[],1);
% s = scatter(x,yPlot, 10, c, 'filled')
% colorbar
% s.MarkerFaceAlpha =0.6;
% % errorbar(x,y, 'LineStyle','none','Marker','.', 'MarkerSize',10)
% hold on 
% xLine = linspace(min(x), max(x), 1000);
% yLine = zeros(size(xLine));
% p=plot(xLine, yLine);
% p.Color = [1,0,0,0.7];
% hold on 
% yLine = zeros(size(xLine)) +1;

% p=plot(xLine, yLine, 'lineStyle','--');
% p.Color = [1,0,0,0.7];
% hold on 
% yLine = zeros(size(xLine)) -1;
% p=plot(xLine, yLine, 'lineStyle','--');
% p.Color = [1,0,0,0.7];
% hold on 

% xlim([min(x)-0.05, max(x)+.05])
% ylim([min(yPlot)-0.05, max(yPlot)+.05])
% xlabel('${{\rm W}_r}^{\rm PM}$~(\AA)', 'interpreter', 'latex')
% ylabel('${({\rm W}_r}^{\rm GP,Voigt} - {{\rm W}_r}^{\rm PM})/err_{\rm Total}$', 'interpreter', 'latex')
% exportgraphics(fig, 'dEWGPVoigt-sw4.png', 'Resolution',800)
% 
% % scatter  (EW(GP,Flux) - EW(PM))/total(err) vs EW(PM)
% fig = figure()
% x = reshape(REW_2796_dr7_PM1GP1_PM(:,:), [], 1); % W_r,PM
% yPlot = reshape(diffEW2796,[],1); % W_r,PM - W_r,GP,flux / errMax
% c = reshape(map_sigma_c4L2_PM1GP1/1e5,[],1);
% s = scatter(x,yPlot,  10, c, 'filled');
% s.MarkerFaceAlpha =0.6;
% cb = colorbar;
% cb.Label.String = '\sigma_{CIV} (kms^{-1})'
% % errorbar(x,y, 'LineStyle','none','Marker','.', 'MarkerSize',10)
% hold on 
% xLine = linspace(min(x), max(x), 1000);
% yLine = zeros(size(xLine));
% p=plot(xLine, yLine);
% p.Color = [1,0,0,0.7];
% hold on 
% yLine = zeros(size(xLine)) +2;
% 
% p=plot(xLine, yLine, 'lineStyle','--');
% p.Color = [1,0,0,0.7];
% hold on 
% yLine = zeros(size(xLine)) -2;
% p=plot(xLine, yLine, 'lineStyle','--');
% p.Color = [1,0,0,0.7];
% hold on 
% 
% xlim([min(x)-0.05, max(x)+.05])
% ylim([min(yPlot)-0.05, max(yPlot)+.05])
% xlabel('${\rm W}_{r,2796}^{\rm PM}$~(\AA)', 'interpreter', 'latex')
% ylabel('$({\rm W}_{r,2796}^{\rm GP,flux} - {\rm W}_{r,2796}^{\rm PM})/err_{\rm Max}$', 'interpreter', 'latex')
% exportgraphics(fig, 'dEW-GPflux-sw4.png', 'Resolution',800)
% nnz(abs(yPlot)<2 & ~isnan(yPlot))/nnz(~isnan(yPlot))

% % scatter  (EW(GP,Flux) - EW(PM))/total(err) vs EW(PM) for 2803A
% fig = figure()
% x = reshape(REW_2803_dr7_PM1GP1_PM(:,:), [], 1); % W_r,PM
% yPlot = reshape(diffEW2803,[],1); % W_r,PM - W_r,GP,flux / errTotal
% c = reshape(map_sigma_c4L2_PM1GP1,[],1);
% s = scatter(x,yPlot,  10, c, 'filled');
% s.MarkerFaceAlpha =0.6;
% cb = colorbar;
% cb.Label.String = '\sigma_{CIV} (kms^{-1})'
% % errorbar(x,y, 'LineStyle','none','Marker','.', 'MarkerSize',10)
% hold on 
% xLine = linspace(min(x), max(x), 1000);
% yLine = zeros(size(xLine));
% p=plot(xLine, yLine);
% p.Color = [1,0,0,0.7];
% hold on 
% yLine = zeros(size(xLine)) +2;

% p=plot(xLine, yLine, 'lineStyle','--');
% p.Color = [1,0,0,0.7];
% hold on 
% yLine = zeros(size(xLine)) -2;
% p=plot(xLine, yLine, 'lineStyle','--');
% p.Color = [1,0,0,0.7];
% hold on 

% xlim([min(x)-0.05, max(x)+.05])
% ylim([min(yPlot)-0.05, max(yPlot)+.05])
% xlabel('${\rm W}_{r,2803}^{\rm PM}$~(\AA)', 'interpreter', 'latex')
% ylabel('$({\rm W}_{r,1500}^{\rm GP,flux} - {\rm W}_{r,1500}^{\rm PM})/err_{\rm Total}$', 'interpreter', 'latex')
% exportgraphics(fig, 'dEW1500-GPflux-sw3.png', 'Resolution',800)
% nnz(abs(yPlot)<2 & ~isnan(yPlot))/nnz(~isnan(yPlot))


% % identifying outsider error ratio spectra 

vs = {'ind_low', 'ind_low_numMgII', 'ind_high', 'ind_high_numMgII', 'ind_0', 'ind_0_numMgII', 'diffEW2796', 'diffEW2803', ...
        'ind_dv_low', 'ind_dv_low_numMgII', 'ind_dv_high', 'ind_dv_high_numMgII', 'ind_dv_mid', 'ind_dv_mid_numMgII', 'Dv_PM1GP1',...
        'ind_EW_large_PM1GP0',  'DV_EW_large_PM1GP0', 'ind_PM1GP0'};
save('ind_err.mat', vs{:})

% vs = {'ind_low', 'ind_low_numMgII', 'ind_high', 'ind_high_numMgII', 'ind_0', 'ind_0_numMgII', 'diffEW2796', 'diffEW2803', ...
%         'ind_dv_low', 'ind_dv_low_numMgII', 'ind_dv_high', 'ind_dv_high_numMgII', 'ind_dv_mid', 'ind_dv_mid_numMgII', 'Dv_PM1GP1',...
%         'ind_EW_large_PM1GP0', 'ind_EW_large_PM1GP0_numMgII', 'DV_EW_large_PM1GP0', 'ind_PM1GP0'};
% save('ind_err.mat', vs{:})

% % SCatter of dv and rew(PM)
% fig = figure();
% box on
% x = reshape(REW_2796_dr7_PM1GP1_PM(:,:), [], 1); % W_r,PM
% y = reshape(Dv_PM1GP1, [], 1);

% s = scatter(x,y, 10, 'filled')
% s.MarkerFaceAlpha = 0.6;
% xlabel('${\rm W}_r^{\rm PM}$', 'interpreter', 'latex')
% ylabel('$\delta v_{\rm PM,GP}~({\rm km.s^{-1}})$', 'interpreter', 'latex')
% exportgraphics(fig, 'scatter-dv-EWPM.png', 'resolution', 800)