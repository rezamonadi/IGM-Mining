 function y_averaged = Averager(y_fine, nAVG, L_org);
    
    % edit: In the avearging prpcess points do not get blended
    % nAVG --> the number of interpolating points needed between 
    % two succesive points in the origianl vector

    % L_org --> the length of the original vector wich was fined and 
    % it is equal to the length of the final averaged vector

    
    % size of the fined vector

    y_averaged = zeros(L_org, 1);
    % initializing the averaged vector
    y_averaged(1,1)   =y_fine(1);
    y_averaged(end,1) = y_fine(end);
    % boundary condition --> the first (last) element of 
    % the averaged vector is the same as the first (last) element of 
    % the fine vector
    for i=2:L_org-1
        sum=0;
        for j=0:nAVG
            sum=sum+y_fine(nAVG*(i-1/2)+i-j);
            % We are averaging all of the elements of y_fine corresponding to 
            % y_org (i)  so we sum nAVG arrays around y_fine((nAVG+1)*i) --> total of 
            % 2*nAVG +1 arrays 
        end
            y_averaged(i,1) = sum/(nAVG + 1);
            % total number of the element within the 
            % averaging windwo is 2*nAVG +1 (the total number 
            % of iterations in the j loop)
    end
    

 end