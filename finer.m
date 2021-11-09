function y_fine = finer(y, n)
    L = length(y);
    L_fine = (L-1)*n +L;
    y_fine = zeros(1,L_fine);
    y_fine(1) = y(1);
    y_fine(L_fine)= y(L);
    
    for i=1:L-1
        h = (y(i+1)-y(i))/(n+1);
        for j=0:n
            y_fine((n+1)*i-j) = y(i) + (n-j)*h;
        end
    end

end