function y = credibleIntervalFunction(theta, L, level)
       
%       L = L - max(L);
%       L = exp(L);
      [L_sorted,ind_L_sorted] = sort(L, 'descend'); % sort L so the max(L) is the 1st one
      theta_L_sorted = theta(ind_L_sorted); % sort theta based on weights of L so MAP is the 1st element
      
      L_sorted = L_sorted/sum(L); % normalize the L
      
      this_sum_L = zeros([length(L),1]);
      for i=1:length(L)
          this_sum_L(i) = sum(L_sorted(1:i));
%           if (this_sum_L>=level)
            % finding the index of sorted L that at least covers level% of credibility
%             break;  
%           end
      end

      [~,i] = min(abs(this_sum_L-level))
      y(1) = min(theta_L_sorted(1:i)); % lower limit of Credible Interval
      y(2) = max(theta_L_sorted(1:i)); % upper limit 
      
     
end
