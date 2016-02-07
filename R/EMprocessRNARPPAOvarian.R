x<-"hello"
require(IdMappingAnalysis)
#### suggested to call: options(stringsAsFactors = FALSE) 

CancerType<-"READ"
Workflow<-"V1"

RNAseqfile_name<-function(x,y){
  z<-paste(as.character(x),"_TCGADATA",as.character(y),sep="")
  return(z)
}

#### suggested to call: options(stringsAsFactors = FALSE) 
assign(RNAseqfile_name(CancerType,Workflow),read.table(file="./QuickStartGuide_Results/AdvancedProcessingResult/CombinedMultiPlatformDataV1.txt",header=TRUE))


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



bootstrap.rnarppa<-Bootstrap(corrData.rnaseq,Fisher=TRUE,verbose=TRUE);

bootModel.rna.rppa <-as.data.frame(bootstrap.rnarppa)
bootModel.rna.rppa<-bootModel.rna.rppa[c(1:48,50:115,117:134),]
postProbs.rnarppa<- NULL
postProbVar.rnarppa<- NULL

EMtest.rnarppa<- fit2clusters.v2(bootModel.rna.rppa$corr, bootModel.rna.rppa$sd^2,psi0Constraint=0, sameV=T,estimatesOnly=F,seed=Random.seed.save)
postProbsOvarian<-as.vector(EMtestOvarian[[1]]/(1+EMtestOvarian[[1]]))  #not needed at the output is a datafrmae from fit2clusters
postProbVarOvarian <-as.vector(EMtestOvarian[[2]])