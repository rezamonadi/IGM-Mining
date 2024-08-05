set_parameters_dr7; 
build_catalog_dr7;
filename ='data/dr7/processed/processed_qsos_tst_N-1250-1610-S-35-115-civWVL.mat';
% filename ='data/dr7/processed/processed_qsos_tst_mask-1-prior-1-OccamRazor-1-nC4-30000-plt-0-MaskinP-0-fixedPr.mat';
load(filename);
load('EW/REW_DR7_sigma_width_4.mat')
EW_C13_1548             =       all_EW1(test_ind,:);
errEW_C13_1548         =       all_errEW1(test_ind,:);
EW_C13_1550             =       all_EW2(test_ind,:);
errEW_C13_1550         =       all_errEW2(test_ind,:);
num_quasars = nnz(test_ind);
% Saving Kathy's data      
variables_to_load= {'all_plate_dr7', 'all_mjd_dr7', 'all_fiber_dr7', ...
'all_QSO_ID', 'all_RA', 'all_DEC', 'all_zqso',...
'all_N_civ','all_z_civ1', 'all_z_civ2', 'all_z_civ3', 'all_RATING', 'c4_QSO_ID', 'all_EW1', 'all_errEW1', 'all_errEW2' };
load(sprintf('%s/catalog', processed_directory(release)), ...
    variables_to_load{:});          

z_qsos             =      all_zqso(test_ind);

Z_C13_1 = all_z_civ1(test_ind,:);
Z_C13_2 = all_z_civ2(test_ind,:);
Z_C13_3 = all_z_civ3(test_ind,:);
Rate_test = all_RATING(test_ind, :);
testID = all_QSO_ID(test_ind);
dv=350;
% num_quasars=1245;
indPM1GP1dv1_inGP = zeros([num_quasars, 7]);
indPM1GP1dv1_inPM = zeros([num_quasars, 17]);
indGP1PM1dv1_inGP = zeros([num_quasars, 7]);
indGP1PM1dv1_inPM = zeros([num_quasars, 17]);
indGP1PM1dv0_inGP = zeros([num_quasars, 7]);
indGP1PM1dv0_inPM = zeros([num_quasars, 17]);
indPM1GPUncertain_inGP = zeros([num_quasars, 7]);
indPM1GPUncertain_inPM = zeros([num_quasars, 17]);
indPM1GPmissed_inPM = zeros([num_quasars,17]);
Dv_PM1GP1          =nan([num_quasars, 17]);
zciv_PM1GP1_PM     = nan([num_quasars,17]);
REW_PM1GP1dv1_PM = nan([num_quasars,17]);
REW_PM1GPmissed = nan([num_quasars,17]);
REW_PM1GPUncertain = nan([num_quasars,17]);
REW_PM0GP1 = nan([num_quasars,7]);
map_sigma_c4L2_PM1GP1 = nan([num_quasars,17]);
diffEW1548 = nan([num_quasars,17]);
quasar_indPM1GPmissed = zeros([num_quasars,1]);
REW_1548_dr7_PM1GP1_PM = nan([num_quasars,17]);
ii=0;
for tr=0.95
    ii=ii+1;
    PM=0;
    GP=0;
    PM1GP1dv1=0; 
    PM1GP1dv0=0; 
    PM0GP1 = 0;
    PM1GP0 = 0;
    PM1GPUncertain=0;
    PM1GPFar =0;
    PM1GP1dv1V2=0;
    PM1GPmissed =0;
    PM1GP50 =0;
    GP1PM1dv1=0;
    GP1PM1dv0=0;
    DZ = nan([num_quasars, 17,7]);
  
    for quasar_ind=1:num_quasars
        for i=1:17
            if(Rate_test(quasar_ind,i)>=2)
                PM = PM+1; % --> 829
            end
        end
        for j=1:max_civ
            if (p_c4(quasar_ind, j)>=tr)
                GP = GP+1; % --> 822
            end
        end
    end
    for quasar_ind=1:num_quasars
        for i=1:17 % building DV matrix
            for j=1:7
                DZ(quasar_ind, i,j) = abs(Z_C13_1(quasar_ind, i) - map_z_c4L2(quasar_ind,j));
                
                if(Rate_test(quasar_ind,i)>=2) && ...
                    (p_c4(quasar_ind, j)>=tr) && ...
                    DZ(quasar_ind,i,j) <=kms_to_z(dv)*(1+Z_C13_1(quasar_ind, i))
                    
                
                    PM1GP1dv1 = PM1GP1dv1 +1;
                    indPM1GP1dv1_inPM(quasar_ind,i) = 1;
                    indPM1GP1dv1_inGP(quasar_ind,j) = 1;
                    REW_PM1GP1dv1_PM(quasar_ind,i) = EW_C13_1548(quasar_ind,i);
                    diffEW1548(quasar_ind,i) = ...
                    (REW_1548_DR7_flux(quasar_ind,j) - EW_C13_1548(quasar_ind,i))/ ...
                    max(ErrREW_1548_flux(quasar_ind,j), errEW_C13_1548(quasar_ind,i));
                    map_sigma_c4L2_PM1GP1(quasar_ind,i) = map_sigma_c4L2(quasar_ind,j);
                    Dv_PM1GP1(quasar_ind, i) = (Z_C13_1(quasar_ind, i) - map_z_c4L2(quasar_ind,j))/(1+Z_C13_1(quasar_ind,i))*speed_of_light/1000;
                    zciv_PM1GP1_PM(quasar_ind,i)= map_z_c4L2(quasar_ind,j);
                    REW_1548_dr7_PM1GP1_PM(quasar_ind, i) = EW_C13_1548(quasar_ind, i);

                end
            end
        end
    
    
        for i=1:17 % building DV matrix
            if(Rate_test(quasar_ind,i)>=2)
                DZ(quasar_ind, i,:) = abs(Z_C13_1(quasar_ind, i) - map_z_c4L2(quasar_ind,:));
                j = find( DZ(quasar_ind,i,:) == min(DZ(quasar_ind,i,:)));
                if size(j>1)
                    j=j(1);
                end
            
                
                if (p_c4(quasar_ind, j)<tr) && ...
                    DZ(quasar_ind, i,j)<= kms_to_z(dv)*(1+Z_C13_1(quasar_ind,i)) && ...
                    indPM1GP1dv1_inPM(quasar_ind,i)==0 && ...
                    indPM1GP1dv1_inGP(quasar_ind,j)==0 
            
                    

                        PM1GPUncertain = PM1GPUncertain + 1;
                        indPM1GPUncertain_inPM(quasar_ind,i) = 1;
                        indPM1GPUncertain_inGP(quasar_ind,j) = 1;
                        REW_PM1GPUncertain(quasar_ind,i) = REW_1548_DR7_flux(quasar_ind,j);
                        
                end

                if (p_c4(quasar_ind, j)<tr) && ...
                    p_c4(quasar_ind,j)>0.5 && ...
                    DZ(quasar_ind, i,j)<= kms_to_z(dv)*(1+Z_C13_1(quasar_ind,i)) && ...
                    indPM1GP1dv1_inPM(quasar_ind,i)==0 && ...
                    indPM1GP1dv1_inGP(quasar_ind,j)==0 
                    

                        PM1GP50 = PM1GP50 + 1;

                        
                        
                end
                
            end
        end
    

        for i=1:17 % building DV matrix
            if(Rate_test(quasar_ind,i)>=2) && ...
                indPM1GP1dv1_inPM(quasar_ind,i)==0 && ...
                indPM1GPUncertain_inPM(quasar_ind,i)== 0  %&& ...
                % all(DZ(quasar_ind, i,:)>kms_to_z(dv)*(1+Z_C13_1(quasar_ind,i)))

                    PM1GPmissed = PM1GPmissed + 1;
                    indPM1GPmissed_inPM(quasar_ind,i) = 1;
                    REW_PM1GPmissed(quasar_ind,i) = REW_1548_DR7_flux(quasar_ind,j);

                
            end
        end
   
        % for i=1:7 % building DV matrix
        %     if(p_c4(quasar_ind,i)>=tr)
            
        %         if all(Rate_test(quasar_ind,:)==-1)

        %             PM0GP1 = PM0GP1 +1;
        %             REW_PM0GP1(quasar_ind,i) = REW_1548_DR7_flux(quasar_ind,i);

        %         end
                
        %     end
        % end
        for j=1:7 % building DV matrix
            if (p_c4(quasar_ind, j)>=tr)
                for i=1:17
                    DZ(quasar_ind, i,j) = abs(Z_C13_1(quasar_ind, i) - map_z_c4L2(quasar_ind,j));
               
                     
                    if DZ(quasar_ind,i,j) <=kms_to_z(dv)*(1+Z_C13_1(quasar_ind, i)) && ...
                       Rate_test(quasar_ind,i)>=2 
                    
                
                        GP1PM1dv1 = GP1PM1dv1 +1;
                        indGP1PM1dv1_inPM(quasar_ind,i) = 1;
                        indGP1PM1dv1_inGP(quasar_ind,j) = 1;
                        % break;
                    end

                    %      GP1PM1dv0 = GP1PM1dv0 +1;
                    %      indGP1PM1dv0_inPM(quasar_ind,i) = 1;
                    %      indGP1PM1dv0_inGP(quasar_ind,j) = 1;
                    %      % REW_PM1GP1dv1(quasar_ind,j) = REW_1548_DR7_flux(quasar_ind,i);
                     
                    %  end


                    
                end
            end
        end
    



                  
    end
% 
end

% fprintf('PM:%d,PM1GP1dv1:%d, PM1GPUncertain:%d, PM1GPmissed:%d, D:%d\n',...
%          PM, PM1GP1dv1, PM1GPUncertain, PM1GPmissed, ...
%          PM-PM1GP1dv1- PM1GPUncertain- PM1GPmissed)
nnz(indGP1PM1dv1_inGP(p_c4>=tr)==0)
nnz(indGP1PM1dv1_inGP)

GP1PM1dv1


% REW of different categories 

h = histogram(reshape(REW_PM1GP1dv1_PM(~isnan(REW_PM1GP1dv1_PM) & REW_PM1GP1dv1_PM>=0),[],1), 'normalization', 'pdf');
h.BinWidth =0.25;
x3 = h.BinEdges;
x3 = x3(2:end) - h.BinWidth/2;
y3 = h.Values;


fig = figure();
clf();

h = histogram(reshape(REW_PM1GPmissed(~isnan(REW_PM1GPmissed) & REW_PM1GPmissed>=0),[],1), 'normalization', 'pdf');
h.BinWidth =0.25;
x1 = h.BinEdges;
x1 = x1(2:end) - h.BinWidth/2;
y1 = h.Values;
h.FaceAlpha = 0.5;
hold on
h = histogram(reshape(REW_PM1GPUncertain(~isnan(REW_PM1GPUncertain) & REW_PM1GPUncertain>=0),[],1), 'normalization', 'pdf');
h.BinWidth =0.25;
x4 = h.BinEdges;
x4 = x1(2:end) - h.BinWidth/2;
y4 = h.Values;
h.FaceAlpha = 0.5;
hold on

h = histogram(reshape(REW_1548_DR7_flux(p_c4>=tr & indGP1PM1dv1_inGP==0 & REW_1548_DR7_flux>=0),[],1), 'normalization', 'pdf');
h.BinWidth=0.25;
x2 = h.BinEdges;
x2 = x2(2:end) - h.BinWidth/2;
y2 = h.Values;
h.FaceAlpha = 0.5;
h.FaceColor = [0.3, 0.6, 0.0];
hold on 

p = plot(x3,y3);
p.LineWidth=3;
p.Color = [0.2, 0.2, 0.2, 0.7]
hold on 
legend({'PM Only (40)', 'GP uncertain (142)', 'GP only (175)', 'PM $\&$ GP (647)'}, 'interpreter', 'latex')
set(gca, 'fontsize', 15)

ylabel('PDF')
xlabel('${\rm W}_{r,1548}$', 'interpreter', 'latex')
exportgraphics(fig, 'REW_category.png', 'resolution', 800)



% % scatter  (EW(GP,Flux) - EW(PM))/total(err) vs EW(PM)
fig = figure()
x = reshape(REW_PM1GP1dv1_PM(:,:), [], 1); % W_r,PM
yPlot = reshape(diffEW1548,[],1); % W_r,PM - W_r,GP,flux / errMax
c = reshape(map_sigma_c4L2_PM1GP1/1e5,[],1);
s = scatter(x,yPlot,  10, c, 'filled');
s.MarkerFaceAlpha =0.6;
cb = colorbar;
cb.Label.String = '\sigma_{CIV} (km.s^{-1})'
hold on 
xLine = linspace(min(x), max(x), 1000);
yLine = zeros(size(xLine));
p=plot(xLine, yLine);
p.Color = [1,0,0,0.7];
hold on 
yLine = zeros(size(xLine)) +2;

p=plot(xLine, yLine, 'lineStyle','--');
p.Color = [1,0,0,0.7];
hold on 
yLine = zeros(size(xLine)) -2;
p=plot(xLine, yLine, 'lineStyle','--');
p.Color = [1,0,0,0.7];
hold on 

xlim([min(x)-0.05, max(x)+.05])
ylim([min(yPlot)-0.05, max(yPlot)+.05])
xlabel('${\rm W}_{r,1548}^{\rm PM}$~(\AA)', 'interpreter', 'latex')
ylabel('$({\rm W}_{r,1548}^{\rm GP,flux} - {\rm W}_{r,1548}^{\rm PM})/err_{\rm Max}$', 'interpreter', 'latex')
exportgraphics(fig, 'dEW-GPflux-sw4.png', 'Resolution',800)
nnz(abs(yPlot)<2 & ~isnan(yPlot))/nnz(~isnan(yPlot))

save('inMissed.mat', 'indPM1GPmissed_inPM');


 %   scatter dv_{PM,GP} vs.. Z_PM
fig = figure();
y = reshape(Dv_PM1GP1, [], 1);
y = y(~isnan(y));
x = reshape(zciv_PM1GP1_PM, [], 1);
x = x(~isnan(x));
p=scatter(x,y, 10, 'filled')
p.MarkerFaceAlpha = 0.7;
hold on 
xLine = linspace(min(x), max(x), 1000);
p = plot(xLine, zeros([1000,1])+150, 'LineWidth',1, 'LineStyle','--')
p.Color = [1,0,0,0.6];
hold on 
xLine = linspace(min(x), max(x), 1000);
p = plot(xLine, zeros([1000,1])-150, 'LineWidth',1, 'LineStyle','--')
p.Color = [1,0,0,0.6];
hold on 
xLine = linspace(min(x), max(x), 1000);
p = plot(xLine, zeros([1000,1]), 'LineWidth',1.5, 'LineStyle','-')
p.Color = [1,0,0,0.6];
ylabel('$\delta {\rm v}_{\rm PM,GP}~(km/s)$', 'interpreter', 'latex')
xlabel('$z_{CIV}^{\rm PM}$', 'interpreter', 'latex')
xlim([min(x)-0.05, max(x)+0.01])
ylim([min(y)-5, max(y)+5])
set(gca, 'fontsize', 12)
exportgraphics(fig, 'dcatterDV_PM1GP1_ZPM.png', 'Resolution', 800)


% SCatter of dv and rew(PM)
fig = figure();
box on
x = reshape(REW_1548_dr7_PM1GP1_PM(:,:), [], 1); % W_r,PM
y = reshape(Dv_PM1GP1, [], 1);
% x = x(~isnan(x));
% y= y(~isnan(x));

s = scatter(x,y, 10, 'filled')
s.MarkerFaceAlpha = 0.6;
xlabel('${\rm W}_{r,1548}^{\rm PM}$', 'interpreter', 'latex')
ylabel('$\delta {\rm v}_{\rm PM,GP}~(km/s)$', 'interpreter', 'latex')
xlim([min(x)-0.01, max(x)+0.01])
ylim([min(y)-5, max(y)+5])
exportgraphics(fig, 'scatter-dv-EWPM.png', 'resolution', 800)

