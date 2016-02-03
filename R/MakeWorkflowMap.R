############function MakeWorkflowMap

############Description: From user Protein Expression Data (or any other primary key) create a IdMap object
### args 1: Protein Expression data 3: first gene expression data 4-?: ... following gene expression data
###
###


#paste(row.names(RNASEQDATA[1:dim(RNASEQDATA)[[1]]/2,]),row.names(RNASEQDATA[(dim(RNASEQDATA)[[1]]/2)+1:dim(RNASEQDATA)[[1]],]),sep=",")
#paste(row.names(RNASEQDATA[1:67,]),row.names(RNASEQDATA[68:134,]),sep=",")

MakeWorkflowMap<-function(gene_expression_data,protein_expression_data){
  rows_gene_data<-dim(gene_expression_data)[[1]]
  primary_symbols<-row.names(protein_expression_data)
  workflow_options<-paste(row.names(gene_expression_data[1:(rows_gene_data/2),]),row.names(gene_expression_data[((rows_gene_data/2)+1):rows_gene_data,]),sep=",")
  IdMap_df<-data.frame(primary_symbols,workflow_options)
  row.names(IdMap_df)<-primary_symbols
  return(IdMap_df)
  }

  