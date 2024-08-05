clear; 
set_parameters_dr12;
release = 'dr12';
% build_catalog_dr12;
catDR7 = load(sprintf('%s/catalog', processed_directory(releasePrior)));
filter_flagsDR7 = load(sprintf('%s/filter_flags', processed_directory(releasePrior)), ...
'filter_flags');
prior_ind = (filter_flagsDR7.filter_flags==0);
all_z_civ_C13 = catDR7.all_z_civ1;
PM_catalog = ...                                   
    load(sprintf('%s/catalog', processed_directory(releasePrior)));
% train_ind -> those LOSs without any problem (filter_flag==0) that does not have
% Civ and useful for training null model (a model without Civ absorption line)
% prior_ind -> those LOSs with Civ absorption and in the half part of test
% test_ind -> second half without any filter flag for testing null and absorption model
% on the LOSs that we know have Civ or not. So, we asses our algorithm in this way.

if (ischar(prior_ind))
    prior_ind = eval(prior_ind);
end

% My prior_ind here is already those OK sight of lines that have CIV
prior.z_qsos  = PM_catalog.all_zqso(prior_ind);
prior.c4_ind = prior_ind;
prior.z_c4 = PM_catalog.all_z_civ1(prior_ind);
 
z_PM_prior = all_z_civ_C13(prior_ind,:);
i=0;

for z_qso=1.7:0.2:5.3
    i=i+1;
    less_ind = (prior.z_qsos < (z_qso + prior_z_qso_increase));
    less_systems = z_PM_prior(less_ind,:);

    this_num_quasars = nnz(less_ind);
    this_p_c4(i,1) = nnz(less_systems(:,1)>0)/this_num_quasars; % at least 1 CIV
    for j=2:max_civ-1
        this_p_c4(i,j) = nnz(less_systems(:,j)>0)/nnz(less_systems(:,j-1)>0);
        
        
    end
    this_p_c4(i,max_civ) = nnz(less_systems(:,max_civ)>0)/nnz(less_systems(:,max_civ-1)>0);

end
z  = 1.7:0.2:5.3;
close all
fig = figure();
clf();
for j=1:7
    y = this_p_c4(:,j);
    ll = plot(z, y , "LineWidth",2);
    if j==1
        ll.Color =[1, 0.1, 0.1, 0.7]
    end
    if j==2
        ll.Color =[0.8500, 0.3250, 0.0980, 0.7]
    end
    if j==3
        ll.Color =[0.4940, 0.1840, 0.5560, 0.7]
    end
    if j==4
        ll.Color =[0.4660, 0.6740, 0.1880, 0.7]
    end
    if j==5
        ll.Color = [0.3010, 0.7450, 0.9330, 0.7]
    end
    if j==7
        ll.Color =[0, 0.4470, 0.7410, 0.7]
    end
    if j==6
        ll.Color =[0.6350, 0.0780, 0.1840, 0.7]
    end


    hold on
    
end
% legend('Pr(1 CIV)', 'Pr(2 CIV)', 'Pr(3 CIV)', 'Pr(4 CIV)',...
%     'Pr(5 CIV)', 'Pr(6 CIV)', 'Pr(7 CIV)', 'Location','southeast',...
%     'color', 'none')
xlabel('$z_{\rm QSO}$', 'interpreter', 'latex')
ylabel('Pr$({\rm M}_{D}$) or Pr$({\rm M}_{S}$)', 'Interpreter','latex')
xlim([1.7 5.3])
text(4.7, this_p_c4(end-3,1),sprintf('$k=%d$', 1), "FontSize",10, 'Color','k', 'Interpreter', 'latex')
text(4.7, this_p_c4(end-3,2),sprintf('$k=%d$', 2), "FontSize",10, 'Color','k', 'Interpreter', 'latex')
text(4.7, this_p_c4(end-3,3),sprintf('$k=%d$', 3), "FontSize",10, 'Color','k', 'Interpreter', 'latex')
text(3.5, this_p_c4(end-9,4),sprintf('$k=%d$', 4), "FontSize",10, 'Color','k', 'Interpreter', 'latex')
text(4.7, this_p_c4(end-3,5),sprintf('$k=%d$', 5), "FontSize",10, 'Color','k', 'Interpreter', 'latex')
text(2.3, this_p_c4(end-15,6),sprintf('$k=%d$', 6), "FontSize",10, 'Color','k', 'Interpreter', 'latex')
text(4.7, this_p_c4(end-3,7),sprintf('$k=%d$', 7), "FontSize",10, 'Color','k', 'Interpreter', 'latex')
set(gca, 'FontSize', 15)
exportgraphics(fig, 'Priors.png', 'Resolution',800)
