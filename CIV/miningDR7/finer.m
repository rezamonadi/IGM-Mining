function y_fine = finer(y, n)
    L = length(y);
    if L==0
        y_fine = y;
    else
        L_fine = (L-1)*n +L; 
        % the length of the fined array is 1 +(n+1)*(L-1) 
        y_fine = zeros(1,L_fine); % initializing y_fine with zeros
        y_fine(1) = y(1); % y_fine in the beginig is the same as original y
        y_fine(L_fine)= y(L);% y_fine at the end is the same as oroginal y

        for i=1:L-1
            h = (y(i+1)-y(i))/(n+1);
            % the distance between each segments in the fined array
            % is distance between two sussecive y elements devided by n+1
            for j=0:n
                y_fine((n+1)*i-j) = y(i) + (n-j)*h; 
                % I checked this with several examples and this works
                % j=0 gives the y_fine falling on the original y
                % j=1 gives y_fine just one step behind y(i+1) 
                % ...
                % j=n gives y_fine where it is equal to y(i)
            end
        end
    end

end