
% test_flux = all_flux(test_ind);
% test_noise2 = all_noise_variance(test_ind);
% test_pixel_mask = all_pixel_mask(test_ind);

% num_quasars = nnz(test_ind);
% Saving Kathy's data                
zQSO_test = all_zqso(test_ind);

% comparing Z and N
% ID = all_QSO_ID(test_ind);
Z_C13 = all_z_civ(test_ind,:);
% N_C13 = all_N_civ(test_ind,:);
% EW1_c13 = all_EW1(test_ind, :);
% EW2_c13 = all_EW2(test_ind, :);
Rate_test = all_RATING(test_ind, :);
% all_wavelengths = all_wavelengths(test_ind);
% fig = figure(); 
dv= 350;
tr=0.95;
% for lw=[0.01, 0.1, 0.5, 1, 2, 3, 4, 5]
for thisMaxCiv = 1:max_civ
    PM=0;
    GP=0;
    PM1GP1=0; 
    PM0GP1 = 0;
    PM1GP0 = 0;
    PMn1GP1 = 0; 
    lw =5;
    indPM1GP1_in_PM = zeros([num_quasars, 17]);
    indPM1GP1_in_GP = zeros([num_quasars, thisMaxCiv]);

    for quasar_ind=1:num_quasars
        for i=1:17
            if(Rate_test(quasar_ind,i)>=2)
                PM = PM+1;
            end
        end
        for j=1:thisMaxCiv
            if (p_c4(quasar_ind, j)>=tr)
                GP = GP+1;
            end
        end

        for i=1:17
            if(Rate_test(quasar_ind,i)>=2)
                for j=1:thisMaxCiv
                    DZ = abs(Z_C13(quasar_ind, i) - map_z_c4L2(quasar_ind, j, 1));
                    if p_c4(quasar_ind, j)>=tr & DZ<kms_to_z(dv)*(1+zQSO_test(quasar_ind))
                        PM1GP1 = PM1GP1 + 1;
                        indPM1GP1(PM1GP1,1) = i;
                        indPM1GP1(PM1GP1,2) = j;
                        break;
                    end
                end
            end
        end

    %     if(all(Rate_test(quasar_ind, :))==-1)
    %         for j=1:max_civ
    %             if (p_c4_all(quasar_ind, j)>=tr)
    %                 PM0GP1 = PM0GP1 + 1;
    %             end
    %         end
    %     end


    end
    fprintf('mx:%d\nPM = %d\nGP = %d\nPM1GP1 = %d\nPurity = %.2f, Completness = %.2f\n',...
            thisMaxCiv, PM, GP, PM1GP1, 100*PM1GP1/GP, 100*PM1GP1/PM);
end
