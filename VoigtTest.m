% Testing Voigt Profile



fN = 1.2;



sigma_samples = offset_sigma_samples*(max_sigma-min_sigma)+min_sigma;

EW = nan(num_C4_samples,1);
iEW=0;
parfor iSample=1:num_C4_samples
     L1 = voigt_iP(padded_wavelengths,...
                    map_z_c4L2(quasar_ind,1), fN*nciv_samples(iSample), 1,...
                    sigma_samples(iSample), padded_sigma_pixels);
      EW(iSample) =trapz(this_unmasked_wavelengths, 1-L1)./z_qso;

    
end

fig = figure();
p=histogram(EW, 'Normalization','pdf', 'NumBins', 20);
bins = p.BinEdges;
% title(sprintf('Center:%.2e, range:%.2e, fN:%.2f, max(EW)=%.3f',...
    % Center, range, fN, max(EW)))
hold on
EW_PM = reshape(all_EW1, 17*26030,1);
EW_PM = EW_PM(EW_PM>0);
p=histogram(EW_PM, 'Normalization','pdf', 'NumBins',20);
p.BinEdges= bins;
hold on 
legend('REW_{GP}', 'REW_{PM}')
xlabel('REW')
ylabel('pdf')
exportgraphics(fig, 'REW_Prior.png', 'Resolution', 800)
% hold on
% % hold on
