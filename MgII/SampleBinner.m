function y = SampleBinner(Weights,...
                          ColumnDensitySamples,...
                          minEdgeColumnDensityBins,...
                          widthColumnDensityBins,...
                          numColumnDensityBins,...
                          NumWeightedSamples)
                          

% This function gets the weights and the samples of 
% column densities and outputs the counts in each bin of 
% column density

randomSampledColumnDensitySamples = randsample(ColumnDensitySamples, NumWeightedSamples, true, Weights);

y = zeros([numColumnDensityBins, 1]);
for i=1:length(randomSampledColumnDensitySamples)
    This_ind = floor((randomSampledColumnDensitySamples(i) - minEdgeColumnDensityBins)/widthColumnDensityBins) + 1;
    y(This_ind) = y(This_ind) + 1;
end

