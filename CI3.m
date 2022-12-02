function y = CI3(theta, dTheta, L, level, MAPtheta)
      assert(length(dTheta)==3, 'Data is not 3D')
      iMax = 10000;
      jMax = 10000;
      kMax = 10000;
      SL = sum(L); % total sume of Likelihood
      [~,indV] = max(L); % initialize index for the credible region is just the single 
      % pixel that the max(L) is residing 
      i=0;
      while sum(L(indV))<level*SL
            i=i+1;

            % in each step a larger rectangle 
            % surrounding the maximum pixel is found 
            indV = (abs(theta(:,1)-MAPtheta(1))<dTheta(1)*i) & ...
                   (abs(theta(:,2)-MAPtheta(2))<dTheta(2)*i) & ...
                   (abs(theta(:,3)-MAPtheta(3))<dTheta(3)*i);
            


      end
      fprintf('L-%f, L/SL-%f\n', sum(L(indV)), sum(L(indV))/SL);
      % flag = 0;
      % for ii=1:iMax
            
      %       for jj=1:jMax
                  
      %             for kk=1:kMax
      %                   indV = (abs(theta(:,1)-MAPtheta(1))<dTheta(1)*ii) & ...
      %                                (abs(theta(:,2)-MAPtheta(2))<dTheta(2)*jj) & ...
      %                                (abs(theta(:,3)-MAPtheta(3))<dTheta(3)*kk);

      %                   if sum(L(indV))>=level*SL
      %                         for i=1:3
      %                               y(i,1) = min(theta(indV,i));
      %                               y(i,2) = max(theta(indV,i));
      %                         end
      %                         flag = 1;
      %                         break;
      %                   end
      %             end
      %             if flag==1
      %                   break;
      %             end
      %       end
      %       if flag==1
      %             break;
      %       end
      % end
      for i=1:3
            y(i,1) = min(theta(indV,i));
            y(i,2) = max(theta(indV,i));
      end
      y(4,1) = nnz(indV);
      
end
