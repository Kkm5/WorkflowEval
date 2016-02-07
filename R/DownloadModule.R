##############TCGA-Assembler download process

RPPARawData<-DownloadRPPAData(traverseResultFile = "./DirectoryTraverseResult_Jul-08-2014.rda", saveFolderName = 
                                 "C:/Users/Mcdade/Desktop/WorkflowEval", cancerType = "READ", assayPlatform = "mda_rppa_core", 
                               inputPatientIDs = NULL);
#########write a function to grep TCGA
experiment_names<-c((RPPARawData[1,])[53:182])

READ_RSEM_NORM_RNASeqRawData<-DownloadRNASeqData(traverseResultFile = "./DirectoryTraverseResult_Jul-08-2014.rda", saveFolderName = 
                                     "C:/Users/Mcdade/Desktop/WorkflowEval", cancerType = "READ", assayPlatform = "RNASeqV2", 
                                   dataType = "rsem.genes.normalized_results", inputPatientIDs = NULL);

GeneExpData = ProcessRNASeqData(inputFilePath = 
                                  "./QuickStartGuide_Results/RawData/READ__unc.edu__illuminahiseq_rnaseqv2__rsem.genes.normalized_results__Jul-08-2014.txt", 
                                outputFileName = "READ__illuminahiseq_rnaseqv2__GeneExp", outputFileFolder = 
                                  "./QuickStartGuide_Results/BasicProcessingResult", dataType = "GeneExp", verType = "RNASeqV2");
