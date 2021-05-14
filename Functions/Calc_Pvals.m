function p_vals = Calc_Pvals(HPAProtResults_a01_score200, N)
    %% Stats for each of 500 top scorers
    P = HPAProtResults_a01_score200{2, 1}.hiv_info.Index;
    % creating null distribution of scores to compare against
    K = HPAProtResults_a01_score200{2, 1}.diff_mat;
    GOI_ind = HPAProtResults_a01_score200{2, 1}.prot_info.Index(1:N);
    ScoreDist = zeros(N,1000);
    for i = 1:1000
        P_prime = setdiff(1:length(K(:,1)), P); % exclude actual Positive set
        y = zeros(length(K(:,1)), 1); 
        % set random diffusion sources
        y(P_prime(randperm(length(P_prime), length(P)))) = 1;
        scores = K*y;
        % calculate scores of 500 GOI
        ScoreDist(:,i) = scores(GOI_ind);
    end

    %% T-test of score against random dist
    p_vals = (sum(ScoreDist > repmat(HPAProtResults_a01_score200{2, 1}.prot_info.Score(1:500), 1, 1000),2))/1000;
end