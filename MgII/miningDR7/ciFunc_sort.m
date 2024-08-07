function y = ciFunc_sort(Weights,...
                         ConfidenceLevel,...
                         sampleColumnDensityCIV, ...
                         sampleSigmaParameterCIV,...
                         sampleRedshiftCIV,...
                         sampleW_r_1548,...
                         sampleW_r_1550)
    % Weights (float num_c4_samples)
    % ConfidenceLevel can be 0.68 for ~ 1-sigma  (float [1])
    % Theta --> parameter space (float num_c4_samples X 3)

    % y--> output is a lower and upper limits for each parameter given the ConfidenceLevel
    y = nan([5,2]);

    % Sorting the weights so that the wight corresponding to 
    % the Maximum a posteriori value is the first element in the sorted vector 
    Weights = Weights./sum(Weights);
    [Weights_sorted, ind_sorted_weights] = sort(Weights, 'descend');


    % Assuring that Weights are summed to 1
    % if sum(Weights) ~=1
        % warning('Weights are not normalized to 1 ...')
    % end

    thisSum = 0; % initializing the current sum of sorted weights
    i_sort = 1; % counter for summing the weights
    while thisSum <= ConfidenceLevel
        thisSum = thisSum + Weights_sorted(i_sort);
        i_sort = i_sort + 1;
        i_max = i_sort;
        if thisSum > ConfidenceLevel & i_sort > 2
            i_max = i_sort - 1;
        end
             

    end
    

    sampleColumnDensityCIV_sorted  = sampleColumnDensityCIV(ind_sorted_weights);
    sampleSigmaParameterCIV_sorted = sampleSigmaParameterCIV(ind_sorted_weights);
    sampleRedshiftCIV_sorted       = sampleRedshiftCIV(ind_sorted_weights);
    sampleW_r_1548_sorted          = sampleW_r_1548(ind_sorted_weights);
    sampleW_r_1550_sorted          = sampleW_r_1550(ind_sorted_weights);
    % Sorting the parameter space based on the weights

    y(1,1) = min(sampleColumnDensityCIV_sorted(1:i_max));
    y(1,2) = max(sampleColumnDensityCIV_sorted(1:i_max));

    y(2,1) = min(sampleSigmaParameterCIV_sorted(1:i_max));
    y(2,2) = max(sampleSigmaParameterCIV_sorted(1:i_max));

    y(3,1) = min(sampleRedshiftCIV_sorted(1:i_max));
    y(3,2) = max(sampleRedshiftCIV_sorted(1:i_max));

    y(4,1) = min(sampleW_r_1548_sorted(1:i_max));
    y(4,2) = max(sampleW_r_1548_sorted(1:i_max));

    y(5,1) = min(sampleW_r_1550_sorted(1:i_max));
    y(5,2) = max(sampleW_r_1550_sorted(1:i_max));


    % fprintf('Sum(Weight)=%.3f, i_sort:%d, L(Wight):%d, CI_W:%.5f\n',thisSum,...
    %          i_max, length(Weights), y(4,2)-y(4,1))


end




