% --- Main script calls all other functions --- %
% must add DataSets and Functions folders to path before running
%% ALDK matrices for various tissue types
load WholeBlood_PPI.mat
PPI = {WholeBlood}; % additional PPIs can be added in cell array
% preallocate cellarray to store results
HPAProtResults_a01_score200 = cell(2,length(PPI));
HPAProtResults_a01_score200(1,:) = {'WholeBlood'};

%% Select HIV-human interactors based on edge threshold cutoff of 200
load AnnotatedProteinLinks.mat
load ENSGtoENSPConvert.mat
HS_int_ENSG = HIVInteractors(HIV_stringproteinannotations, ...
    STRING_ProteinLinkes_HomoSapiens, string9606ENSGENSP10allT, 200);

%% Runs diffusion analysis and saves results to HPAProtResults
for i = 1:length(PPI)   
    [prot_info, hiv_info, K] = DiffScore(PPI{i}, HS_int_ENSG, 0.1); % alpha = 0.1
    HPAProtResults_a01_score200{2,i} = struct('prot_info', prot_info, 'hiv_info', hiv_info, 'diff_mat', K);
end

% sort results by score in descending order
for k = 1:length(PPI)
    HPAProtResults_a01_score200{2,k}.prot_info = sortrows(HPAProtResults_a01_score200{2,k}.prot_info, 3, 'descend');
end

%% Write excel file of diffusion from sources to high scoring nodes
HeatMapSubset(HPAProtResults_a01_score200, 'HeatMapSubset_05052021.xlsx')

%% Calculate p-value for top N scorers 
p_vals = Calc_Pvals(HPAProtResults_a01_score200, 500); % p val calculation for top 500 scorers