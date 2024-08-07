function y = bin_averager(data_x,data_y, nBin)
    
    h = (max(data_x) - min(data_x))/nBin;

    y = zeros([nBin,2]);
    for i=1:nBin;
        lo_x = min(data_x) + (i-1)*h;
        hi_x = lo_x + h; 
        y(i,2) = mean(data_y(data_x>=lo_x & data_x<=hi_x));
        y(i,1) = lo_x + h/2;
    end
end
