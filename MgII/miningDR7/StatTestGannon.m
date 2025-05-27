%% StatTesting - Gannon Hughes
% Loading Data
        load('data/dr7/processed/processed_sigma_125_binednorm_1%-Masking-All-samp-20k.mat');


% Indicies of Randomly selected quasars
       inds = find((test_ind == 1));
       used_QSOS = all_QSO_ID(test_ind);


% Initializing Comparison Variables
       PM1GP1=0;
       PM1GP0 = 0;
       GP1PM0 = 0;
       PMn1GP1 = 0;
       tr = 0.85;


% Initializing Arrays
        Z_Seyffert = all_z_MgII2(test_ind,:);
       Z_MgII_2_sort = map_z_MgIIL2 .* (p_MgII > tr);
       Rate_test = all_RATING(test_ind, :);


% Finding Number of absorbers PM found
   PM_ind = find((Rate_test>=2)==1);
   PM = nnz(Rate_test>=2);


% Finding the 1% indicie (1-694) of QSO where absorber was found
   [row,col] = find((Rate_test>=2)==1);
   PM_QSO_ind = sort(row)';
   unique_PM_QSO_ind = unique(PM_QSO_ind);


% Finding the QSO catalog indice (1-79595) of QSO where absorber was found
   PM_QSO_number = inds(PM_QSO_ind);
   PM_QSO_number = unique(PM_QSO_number);


% Number of Absorbers PM found in each Quasar Spectra
   x = [row col];
   [unique_vals, ~, id_x] = unique(x(:,1));
   max_vals = accumarray(id_x, x(:,2), [], @max);
   new_x = [unique_vals, max_vals];


% Displays PM_Absorbers Table (# of Quasars Found, QSO_ID (1-79595), 1%_Indicie (1-694), Number of Absorbers Found)
   num = [linspace(1,length(new_x),length(new_x))];
   PM_Absorbers = [num',PM_QSO_number, new_x];
   PM_Absorbers;


% Counting GP Detected Absorbers
    GP_ind = find((p_MgII >= tr));
    GP = nnz(GP_ind);


% Finding the 1% indicie (1-694) of QSO where absorber was found
    [row,col] = find((p_MgII >= tr)==1);
    GP_QSO_ind = sort(row)';
    unique_GP_QSO_ind = unique(GP_QSO_ind);


% Finding the (1-79595) indicie for GP
   GP_QSO_number = inds(unique_GP_QSO_ind);


% Number of Absorbers GP found in each Quasar
    x = [row col];
    new_x = zeros(length(unique_GP_QSO_ind),2);
    for ii = 1:length(unique_GP_QSO_ind)
        y = nnz(p_MgII(unique_GP_QSO_ind(ii),:) >= tr);
        new_x(ii,1) = unique_GP_QSO_ind(ii);
        new_x(ii,2) = y;
    end


 % Displays GP_Absorbers Table (#, QSO_ID (1-79595), 1%_Indicie (1-694), Absorbers Found)
       num = [linspace(1,length(new_x),length(new_x))];
       GP_Absorbers = [num',GP_QSO_number, new_x];
       GP_Absorbers;


% PM1GP1
zc = 0.005; %redshift condition
for ii = 1:length(PM_Absorbers)
    for jj = 1:length(GP_Absorbers)
        if PM_Absorbers(ii,3) == GP_Absorbers(jj,3) % Finds when the indicie (0-694) are equal for GP and PM
            for aa = 1:PM_Absorbers(ii,4)
                val = Z_Seyffert(PM_Absorbers(ii,3),aa);
                if  any(abs(Z_MgII_2_sort(GP_Absorbers(jj,3),:) - val) <= zc) % if the redshift for a given absorber in PM matches one in GP spectra
                                PM1GP1 = PM1GP1 + 1;
                end
            end
        end
    end
end

PM1GP1;
P = PM1GP1/GP
C = PM1GP1/PM


% PM1GP0 PM1GP0 PM1GP0 PM1GP0 PM1GP0 PM1GP0 PM1GP0 PM1GP0 PM1GP0 PM1GP0

% Array of (1-694) indicies that appear in GP but not PM, and vice versa
    GPnotinPM = setdiff(GP_Absorbers(:,3),PM_Absorbers(:,3));
    PMnotinGP = setdiff(PM_Absorbers(:,3),GP_Absorbers(:,3));
    diff_inds = union(GPnotinPM,PMnotinGP);

PM1GP0_GP_inds = [];
GP_missed_abs_z = [];
% Takes the spectra indicies (1-694) that PM found absorbers in BUT GP
% didnt find absorbers in, and adds all the absorbers found in those specs
    for ii = 1:length(PM_Absorbers(:,3))
        for jj = 1:length(PMnotinGP)
            if PM_Absorbers(ii,3) == PMnotinGP(jj)
                PM1GP0 = PM1GP0 + PM_Absorbers(ii,4);
                PM1GP0_GP_inds(end+1) = PM_Absorbers(ii,3);
                GP_missed_abs_z(end+1,1:PM_Absorbers(ii,4)) = Z_Seyffert(PM_Absorbers(ii,3),1:PM_Absorbers(ii,4));
            end
        end
    end

% Not taking into account when PM found DIFFERENT absorbers than GP in a given 
% spectra (when I decrease tr, PM1GP0 for this section goes up)
    for ii = 1:length(PM_Absorbers)
        for jj = 1:length(GP_Absorbers)
            if PM_Absorbers(ii,3) == GP_Absorbers(jj, 3)
            % if the ID (0-694) is equal for both PM and GM need to check if they found any different absorbers
                 for ff = 1:PM_Absorbers(ii,4)
                        val = Z_Seyffert(PM_Absorbers(ii,3), ff);  % val is a specific redshift of a PM absorber
                        if ~any(abs(Z_MgII_2_sort(GP_Absorbers(jj,3),:) - val) <= zc) % if no redshift values in GP for the same quasar match PM absorber
                              PM1GP0 = PM1GP0 + 1;

                              if jj == 1
                                  PM1GP0_GP_inds(end+1) = GP_Absorbers(jj,3);
                                  GP_missed_abs_z(end+1,1) = val;

                              elseif PM1GP0_GP_inds(end) ~= GP_Absorbers(jj,3)
                                  PM1GP0_GP_inds(end+1) = GP_Absorbers(jj,3);
                                  GP_missed_abs_z(end+1,1) = val;
                                  a = 1;
                              elseif PM1GP0_GP_inds(end) == GP_Absorbers(jj,3)
                                  a = a+1;
                                  GP_missed_abs_z(end,a) = val;

                              end
                       end

                 end
             end
         end
    end
PM1GP0

missed_abs = [PM1GP0_GP_inds' GP_missed_abs_z];
j = 1:1:length(missed_abs);
missed_abs = [j' PM1GP0_GP_inds' GP_missed_abs_z];

% GP1PM0 GP1PM0 GP1PM0 GP1PM0 GP1PM0 GP1PM0 GP1PM0 GP1PM0 GP1PM0 GP1PM0 GP1PM0 GP1PM0 GP1PM0 GP1PM0
% Takes the spectra indicies (1-694) that GP found absorbers in BUT PM
% didnt find absorbers in, and adds all the absorbers found in those specs
    for ii = 1:length(GP_Absorbers(:,3))
        for jj = 1:length(GPnotinPM)
            if GP_Absorbers(ii,3) == GPnotinPM(jj)
                GP1PM0 = GP1PM0 + GP_Absorbers(ii,4);
            end
        end
    end
%Not taking into account when PM found DIFFERENT absorbers than GP in a given
%spectra
GP1PM0_PM_inds_diff_abs = [];
    for ii = 1:length(PM_Absorbers)
        for jj = 1:length(GP_Absorbers)
            if PM_Absorbers(ii,3) == GP_Absorbers(jj, 3)
            % if the ID (0-694) is equal for both PM and GM need to check if they found any different absorbers
                 for ff = 1:length(Z_MgII_2_sort(GP_Absorbers(jj,3), :))
                        val = Z_MgII_2_sort(GP_Absorbers(jj,3), ff);
                        if val~= 0 && ~isnan(val) && ~any(abs(Z_Seyffert(PM_Absorbers(ii,3),:) - val) <= zc) % if val isnt nan and if no PM absorbers match the redshift absorber
                              GP1PM0 = GP1PM0 + 1;
                              GP1PM0_PM_inds_diff_abs(end+1) = PM_Absorbers(ii,3);
                        end
                 end
             end
         end
    end
GP1PM0;
%% SDSS QSO IDS For PM1GP0_GP_inds
used_QSOS(PM1GP0_GP_inds);

missed_abs
used_QSOS(missed_abs(:,2))



