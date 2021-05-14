function HeatMapSubset(HPAProtResults_a01_score200, filename)
%% Heatmap of diffusion from sources to high scoring nodes. Biclusters?
    col = HPAProtResults_a01_score200{2,1}.hiv_info.Index;
    row = HPAProtResults_a01_score200{2,1}.prot_info.Index(1:500);
    matSubset = HPAProtResults_a01_score200{2, 1}.diff_mat(row,col);
    %%
    matSubset_RowNorm = (matSubset)./repmat(sum(matSubset,2),1,length(col));
    col_label = HPAProtResults_a01_score200{2,1}.hiv_info.ID;
    row_label = HPAProtResults_a01_score200{2,1}.prot_info.ID(1:500);
    T = array2table(matSubset_RowNorm);
    T.Properties.VariableNames = col_label;
    T.Properties.RowNames = row_label;
    writetable(T,filename,'WriteRowNames',true)
end