function [a_coef, b_coef, R, a, b]= errFitter(x,y,xe,ye, NumSample);
    a_coef = nan([NumSample,1]);
    b_coef = nan([NumSample,1]);
    DimData = length(x);
    x_sample = nan([DimData, NumSample]);
    y_sample = nan([DimData, NumSample]);
    for i=1:DimData
        % generate a Gaussian sample from each input data point 
        % given their error as sigma and the value of data point as mu
        x_sample(i,:) = randn([1,NumSample])*xe(i) + x(i); 
        y_sample(i,:) = randn([1,NumSample])*ye(i) + y(i);
    end
    for i=1:NumSample
        [fit,S] = polyfit(x_sample(:,i), y_sample(:,i),1);
        a_coef(i) = fit(1);
        b_coef(i) = fit(2);
        R(i) = S.normr;
    end

    a(2) = quantile(a_coef, 0.05);
    a(1) = median(a_coef);
    a(3) = quantile(a_coef,0.95);

    b(2) = quantile(b_coef, 0.05);
    b(1) = median(b_coef);
    b(3) = quantile(b_coef,0.95);
    
end

        
    
