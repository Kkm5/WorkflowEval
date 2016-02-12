#############Function to create a file name for each RNAseq download

RNAseqfile_name<-function(x,y){
  z<-paste(as.character(x),"_TCGADATA",as.character(y),sep="")
  return(z)
}

RNA_File<-RNAseqfile_name("READ","V1")

RNA_File


eval(RNAseqfile_name("READ","V1"))
