function HS_int_ENSG = HIVInteractors(HIV_stringproteinannotations, STRING_ProteinLinkes_HomoSapiens,string9606ENSGENSP10allT, threshold_value)
    %% Find all human interactors with 24 HIV proteins listed in 
    %  HIV_stringproteinannotations table
    HIVProt = categorical(HIV_stringproteinannotations.identifier);
    network = STRING_ProteinLinkes_HomoSapiens;
    LiA1 = ismember(network.Protein1, HIVProt);
    LiA2 = ismember(categorical(network.Protein2), HIVProt);

    % check that only Human-Virus interactions
    sum(LiA1) > 0 && sum(LiA2) > 0; % no virus-virus based on indices

    %% Human interactors vector
    com_score = double(string(network.Combined_Score));
    L_ind = (((com_score>=threshold_value) + LiA2) == 2);
    HS_int = unique(network.Protein1(L_ind));

    %%
    % load(ENSGtoENSPConvert.mat)
    HS_int_edit = cellstr(HS_int);
    % remove 9606. from HS_int
    for i = 1:length(HS_int)
        HS_int_edit{i} = HS_int_edit{i}(6:end);
    end
    %% Find in String network
    HS_int_ENSG = HS_int_edit;
    for j = 1:length(HS_int)
        HS_int_ENSG{j} = string9606ENSGENSP10allT.ENSG(ismember(string9606ENSGENSP10allT.ENSP, HS_int_edit{j}));
    end
    HS_int_ENSG = string(HS_int_ENSG);
end