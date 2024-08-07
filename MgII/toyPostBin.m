clear
clc
nMAP=250;
nSample =1000;
nBin=20;
MAP_values = rand([nMAP,1])*10 -5; % MAP values are randomly selected between 5 and 10
Posterior_Samples = zeros([nMAP,nSample]);
random_dist = randn(nMAP, nSample)*0.1;
 for i=1:nMAP
    Posterior_Samples(i,:) = MAP_values(i,1) + random_dist(i,:); 
    % Creating a Gaussian dist around the MAP values
    % with a variance of 0.4 
end


h = (max(MAP_values)-min(MAP_values))/nBin;

Sample_hist_all = zeros([nBin-1,nSample]);
for nB=1:nBin
    xBins(nB) = min(MAP_values) + (nB-1)*h;
end

MAP_hist = zeros([nBin-1,1]);
for nB=1:nBin-1
    for nM =1:nMAP
        if  (MAP_values(nM)>=xBins(nB) & MAP_values(nM)<=xBins(nB+1))
                MAP_hist(nB) = MAP_hist(nB) + 1;
        end
    end
end

for nS=1:nSample
    % bLow = min(MAP_values);
    for nB=1:nBin-1
        for nM =1:nMAP
            if  (Posterior_Samples(nM,nS)>=xBins(nB) & Posterior_Samples(nM,nS)<=xBins(nB+1))
                Sample_hist_all(nB,nS) = Sample_hist_all(nB,nS) + 1;
            end
        end
    end
    % bHigh
end

x_mid_Bins = xBins(2:end) - h/2;

fig = figure();
for i=1:nSample
    plot(x_mid_Bins, Sample_hist_all(:,i), '.','Color','r')
    hold on 
end
plot(x_mid_Bins, MAP_hist)
