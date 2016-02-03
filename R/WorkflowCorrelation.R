x<-"hello"
require(IdMappingAnalysis)
#### suggested to call: options(stringsAsFactors = FALSE) 

################RUN ONLY ONCE UNTIL CORRECTED
Formatted_primary_data<-add_row_names_col(RPPADATA)
Formatted_secondary_data<-add_row_names_col(RNASEQDATA,exp_colname = "RNAseqWorkflow")
#################################

primaryIDs <- IdMapBase$primaryIDs(Formatted_primary_data);
secondaryIDs <- IdMapBase$primaryIDs(Formatted_secondary_data);

###########FIX THE ORDER OF ARGS ITS IS ANNOYING
idmap_rnaseqworkflow_df<-MakeWorkflowMap(Formatted_secondary_data,Formatted_primary_data)

#uniquePairs.rnaseq<-UniquePairs(idmap_rnaseqworkflow_df,primaryKey="Symbol",secondaryKey="RNAseqWorkflow")

IdMap_rna_rppa<-IdMap(DF=idmap_rnaseqworkflow_df,name="IdMap_rna_rppa", primaryKey="Symbol",secondaryKey="RNAseqWorkflow")


######################### unique pairs from the jointIdMap
#uniquePairsOvarian <- as.UniquePairs.IdMap(IdMap_NA_RPPA,secondaryIDsOvarian)

uniquePairs_workflow <- as.UniquePairs.IdMap(IdMap_rna_rppa,secondaryIDs)
corrData.rnaseq <- CorrData(uniquePairs_workflow,Formatted_primary_data,Formatted_secondary_data)
corr.rnarppa.try1 <- Corr(corrData.rnaseq,method="spearman",verbose=TRUE);
