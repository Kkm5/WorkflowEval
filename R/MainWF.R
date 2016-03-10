
require(IdMappingAnalysis)
options(stringsAsFactors = FALSE)

#### suggested to call: options(stringsAsFactors = FALSE) 

################ExperimentDesign Class S3 object methods must be able to:
###This is list of executed workflow steps for a given object
####example  


#################ExpressionFormat Class S3 object methods must be able to 
#1 recognize identifiers and create primary ids on the input
#2 parse the data to include only subjects and identifiers that are present in both objects
#3 create an expanded secondary list where identifiers are retagged with the workflow ID

################ method to ExpressionFormat Class



#Add row to the front of the data to specify symbol

modify_data <-function(dataset, exp_colname="Symbol"){
  modified_dataset<-cbind(row.names(dataset),dataset)
  names(modified_dataset)[1]<-paste(exp_colname)
  return(modified_dataset)
}




WorkflowData <- function(...) {
  input_list <- list(...)
  output_list <- lapply(X=input_list, function(x) {modify_data(x)})
  class(output_list)<-"WorkflowData"
  return(output_list)
}

cname <-function(x,...) {
  args<-list(...)
  object<-list(attribute.name= "something",...)
  class(object)<-"cname"
  return(object)
}

myobject<-cname(x,2,3,4)

Dataset.Collection<-WorkflowData(RPPADATA,RNASEQDATAV1,RNASEQDATAV2)

print.WorkflowData<-function(WFdata) {
  length(WFdata)
}



#Formatted_primary_data<-add_row_names_col(RPPADATA)
#Formatted_secondary_data<-add_row_names_col(RNASEQDATA,exp_colname = "RNAseqWorkflow")
#################################

extract.identifiers<- function(...){
  input_list <- list(...)
  output_list <-lapply(X=input_list, function(x) {row.names(x)})
  return(output_list)
}

Identifier.Collection<-extract.identifiers(RPPADATA,RNASEQDATAV1,RNASEQDATAV2)

#primaryIDs <- IdMapBase$primaryIDs(Formatted_primary_data);
#secondaryIDs <- IdMapBase$primaryIDs(Formatted_secondary_data);



###plan split in RNASEQDATA in half for testing, remove tags, re-write MakeWorkFlowMap for multiple datasets
####Make a funtion that splits the similar data into chunks
#RNASEQDATAV1<-RNASEQDATA[1:67,]
#RNASEQDATANAMESV1<-strsplit(row.names(RNASEQDATAV1),"_")
#row.names(RNASEQDATAV1)<-sapply(RNASEQDATANAMESV1, "[", 1)
#RNASEQDATAV2<-RNASEQDATA[68:134,]
#RNASEQDATANAMESV2<-strsplit(row.names(RNASEQDATAV2),"_")
#row.names(RNASEQDATAV2)<-sapply(RNASEQDATANAMESV2, "[", 1)
####################################################################################
assign.status<-function(x,status,data.type,version){
  if (status == "workflow_option")
    status_tag<-paste(row.names(x),"_WFO","_",as.character(data.type),"_",as.character(version),sep="")
  if (status == "driver")
    status_tag<-paste(row.names(x),"_DRIVER",sep="")
  return(status_tag)
}


######create a function to carry out this effect
merge.options<-function(...){
    input_list <- list(...)
    output_list <-lapply(X=input_list, function(x) {})
    return(output_list)
  }

P<-RPPADATA
RS1<-RNASEQDATAV1
RS2<-RNASEQDATAV2

row.names(P)<-assign.status(RPPADATA,status="driver",data.type="P",1)
row.names(RS1)<-assign.status(RNASEQDATAV1,status="workflow_option",data.type="RS",1)
row.names(RS2)<-assign.status(RNASEQDATAV2,status="workflow_option",data.type="RS",2)

Merged.options<-rbind(P,RS1,RS2)


make.workflow.map <- function(Merged.options){
  drivers<-row.names(Merged.options[sapply(strsplit(row.names(Merged.options),"_"), "[", 2) == "DRIVER",])
  workflow_options_data<-Merged.options[sapply(strsplit(row.names(Merged.options),"_"), "[", 2) == "WFO",]
  workflow_options<-row.names(workflow_options_data[sapply(strsplit(row.names(workflow_options_data),"_"), "[", 2) == "WFO",])
  imax<-max(unique(sapply(strsplit(row.names(workflow_options_data),"_"), "[", 4)),na.rm = TRUE)
  imin<-min(unique(sapply(strsplit(row.names(workflow_options_data),"_"), "[", 4)),na.rm = TRUE)
  workflow_options_merged<-paste(row.names(workflow_options_data[sapply(strsplit(row.names(workflow_options_data),"_"), "[", 4) == imin,]),row.names(workflow_options_data[sapply(strsplit(row.names(workflow_options_data),"_"), "[", 4) == imax,]),sep=",") 
  WorkflowMap<-data.frame(drivers,workflow_options_merged)
  ###when your ready class(Merged.options)<- "WorkflowMap"
}


  ####use unique...
  #  workflow_options<-paste(row.names(gene_expression_data[1:(rows_gene_data/2),]),row.names(gene_expression_data[((rows_gene_data/2)+1):rows_gene_data,]),sep=",")
  #  IdMap_df<-data.frame(primary_symbols,workflow_options)
  #  row.names(IdMap_df)<-primary_symbols
  #  return(IdMap_df)
  #}
#MakeWorkflowMap<-function(gene_expression_data,protein_expression_data){
#  rows_gene_data<-dim(gene_expression_data)[[1]]
#  primary_symbols<-row.names(protein_expression_data)
#  workflow_options<-paste(row.names(gene_expression_data[1:(rows_gene_data/2),]),row.names(gene_expression_data[((rows_gene_data/2)+1):rows_gene_data,]),sep=",")
#  IdMap_df<-data.frame(primary_symbols,workflow_options)
#  row.names(IdMap_df)<-primary_symbols
#  return(IdMap_df)
#}



###########FIX THE ORDER OF ARGS ITS IS ANNOYING
idmap_rnaseqworkflow_df<-MakeWorkflowMap(Formatted_secondary_data,Formatted_primary_data)

#uniquePairs.rnaseq<-UniquePairs(idmap_rnaseqworkflow_df,primaryKey="Symbol",secondaryKey="RNAseqWorkflow")

IdMap_rna_rppa<-IdMap(DF=idmap_rnaseqworkflow_df,name="IdMap_rna_rppa", primaryKey="Symbol",secondaryKey="RNAseqWorkflow")




######################### unique pairs from the jointIdMap
#uniquePairsOvarian <- as.UniquePairs.IdMap(IdMap_NA_RPPA,secondaryIDsOvarian)

uniquePairs_workflow <- as.UniquePairs.IdMap(IdMap_rna_rppa,secondaryIDs)
corrData.rnaseq <- CorrData(uniquePairs_workflow,Formatted_primary_data,Formatted_secondary_data)
corr.rnarppa.try1 <- Corr(corrData.rnaseq,method="spearman",verbose=TRUE)

bootstrap.rnarppa<-Bootstrap(corrData.rnaseq,Fisher=TRUE,verbose=TRUE)

bootModel.rna.rppa <-as.data.frame(bootstrap.rnarppa)
bootModel.rna.rppa<-bootModel.rna.rppa[c(1:96,99:134),]
postProbs.rnarppa<- NULL
postProbVar.rnarppa<- NULL

EMtest.rnarppa<- fit2clusters.RNAseq(bootModel.rna.rppa$corr, bootModel.rna.rppa$sd^2,psi0Constraint=0, sameV=T,estimatesOnly=F,seed=Random.seed.save)
postProbsOvarian<-as.vector(EMtestOvarian[[1]]/(1+EMtestOvarian[[1]]))  #not needed at the output is a datafrmae from fit2clusters
postProbVarOvarian <-as.vector(EMtestOvarian[[2]])


jointIdMap_from_dataset <- JointIdMap(examples$identDfList[c(6)], primaryIDs_from_dataset, verbose=FALSE);


######################### unique pairs from the jointIdMap
uniquePairs <- as.UniquePairs(getUnionIdMap(jointIdMap_from_dataset,verbose=FALSE),secondaryIDs);


######################### data for correlations are placed in a CorrData object
corrData <- CorrData(uniquePairs,msmsExperimentSet,mrnaExperimentSet,verbose=FALSE);


######################### Joint unique pairs
jointUniquePairs <- JointUniquePairs(uniquePairs,getIdMapList(jointIdMap_from_dataset,verbose=FALSE),verbose=FALSE);


######################## correlation values
corr <- Corr(corrData,method="pearson",verbose=FALSE);
corr.spearman <- Corr(corrData,method="spearman",verbose=FALSE);

corrSet <- getCorr(jointUniquePairs,corr,
                   groups=c("union","EnVision_Q","NetAffx_Q","DAVID_Q","DAVID_F"),
                   full.group=TRUE,verbose=FALSE);


pairs <- jointUniquePairs$as.data.frame()

##########################
pairMatchTable<-pairs
#bootModel.fulldata<-as.data.frame(bootstrap)
bootModel <-as.data.frame(bootstrap)
#load("dotRandom.seed.rda")
Random.seed.save <- sync.Random.seed <- .Random.seed

bootMergedWithPairs = merge(data.frame(postProbs=postProbs, postProbVar=postProbVar,bootModel), pairMatchTable)


