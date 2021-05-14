function [Prot_info, krog_info, K] = DiffScore(PPI, K_ens, alpha)
    % Creates Regularized Laplacian matrix to model diffusion in the
    % network and muliplies K by y to show diffusion from select nodes
    % (denoted as 1s in vector y) throughout the remainder of the network
    
    % PPI -- 2 column string array of ENSEMBL gene IDs indicating protein
    % interactions
    % K_ens -- diffusion source nodes in ENSEMBL gene format
    % alpha -- random walk parameter (between 0 and 1)
    G_list = unique(PPI); % create gene list from unique PPI entries
    % preallocate indices and adjacency matrix
    PPI_ind = zeros(size(PPI)); 
    A = zeros(length(G_list));
    tic % track time for large datasets
    % find indices in list
    for i = 1:length(G_list)
        PPI_ind(PPI == G_list(i)) = i;
    end
    % assign interactions 1 in A
    for ii = 1:length(PPI_ind)
        A(PPI_ind(ii,1),PPI_ind(ii,2)) = 1;
        A(PPI_ind(ii,2),PPI_ind(ii,1)) = 1;
    end
    toc
    any(A(find(eye(length(G_list))==1)==1)) %Verify that no self-interaction included
    % calculate Regularized Laplacian
    D = diag(sum(A,2));
    L = D - A;
    I = eye(length(G_list));
    K = inv(I + alpha*L);
    y = ismember(G_list, K_ens);
    y_score = (I + alpha*L)\y; % equivalent to K*y
    gene_score = y_score(y == 0);
    
    % assign results to output variables
    Prot_info = table(G_list(y == 0), find(y == 0), gene_score, 'VariableNames',{'ID','Index','Score'});
    krog_info = table(G_list(y==1), find(y==1), 'VariableNames',{'ID','Index'});
end