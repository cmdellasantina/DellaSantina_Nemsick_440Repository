## HIV_TissueNetResults
Files used to run diffusion analysis on tissueNet GTex PPI for whole blood and store all results in a cell array of structs. 

## Installation
Code written in MATLab_R2021a and no add-ons required. Must add all folders in directory to path when running main script (DiffusionAnalysis_MainScript.m)

## DiffusionAnalysis_MainScript.m 
Loads all variables and calls all functions necessary to generate all results found in HPAProtResults_a01_score200.mat and generate excel file with data found in DiffusionMatrixSubset_ENStoGeneNameConv.xlsx.

## Functions
HIVInteractors.m - Finds all human interactions with 24 HIV proteins listed in HIV_stringproteinannotations table in AnnotatedProteinLinks.mat.Takes 4 inputs: (1) HIV_stringproteinannotations - a 24x789 table from STRING.Viruses 'Human immunodeficiency virus type 1 [HIV-1] group M subtype B (isolate HXB2)'; (2) STRING_ProteinLinkes_HomoSapiens - 11437065x3 table where, for each row, columns 1 and 2 indicate two proteins that interact and column 3 their combined confidence score; (3) string9606ENSGENSP10allT - 23448x2 table of the conversion between ensemble gene and ensemble protein ID for all nodes in string; (4) threshold_value - minimum numerical value imposed on combined score for an edge to be included in analysis. Output is HS_int_ENSG a string array of HIV-Human interactions in ensemble gene ID.

DiffScore.m - Performs network diffusion analysis using regularized laplacian matrix. Takes PPI, Seed labels and alpha value as input. Alpha is a parameter that roughly corresponds to the length of the random walk from a diffusion source throughout the network. For this analysis it is set to 0.1. More detailed explanation of algorithm is provided in: Lancour D, Naj A, Mayeux R, Haines JL, Pericak-Vance MA, et al. (2018) One for all and all for One: Improving replication of genetic studies through network diffusion. PLOS Genetics 14(4): e1007306.

HeatMapSubset.m - Generates excel file of diffusion matrix subset used to generate Figure 2 in accompanying report. 2 inputs: (1) Struct output from DiffScore.m; (2) excel filename for to save under.

Calc_Pvals.m - Generates and tests scores calculated in DiffScore.m agains null distribution. Inputs are result from DiffScore and number of top scorers for which p-values should be calculated (N). Returns vector p_vals with statistical significance of proteins in DiffScore output.

## DataSets
AnnotatedProteinLinks.mat - HIV_stringproteinannotations and STRING_ProteinLinkes_HomoSapiens (described in HIVInteractors.m).

ENSGtoENSPConvert.mat - string9606ENSP10allT (described in HIVInteractors.m).

WholeBlood_PPI.mat - Protein-Protein interaction network from TissueNet 2.0 for whold blood. Each row represents and interaction between the two proteins in the row. Edges are unweighted and all IDs are in ENSEMBL ID format.

## Results
DiffusionMatrixSubset_ENStoGeneNameConv.xlsx - 3 sheet excel file with diffusion matrix subset output by HeatMapSubset.m labeled with the gene name conversions on the subsequent sheets. First column of first sheet contains p-values calculated for protein in each row with Calc_Pvals.m
