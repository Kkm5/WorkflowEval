############## This code checks the names of the TCGA cancer types for use in the RTCGA package, must install magrittr package to run
(cohorts <- infoTCGA() %>% 
   rownames() %>% 
   sub("-counts", "", x=.))

############# code to check for available datasets

checkTCGA("DataSets", "OV")

############# This code creates a directory for the RTCGA package

# dir.create( "data2" )
releaseDate <- "2015-08-21"
sapply( cohorts, function(element){
  tryCatch({
    downloadTCGA( cancerTypes = element, 
                  dataSet = "rnaseqv2__illuminahiseq_rnaseqv2__unc_edu__Level_3__RSEM_genes_normalized__data.Level",
                  destDir = "data2/", 
                  date = releaseDate )},
    error = function(cond){
      cat("Error: Maybe there weren't rnaseq data for ", element, " cancer.\n")
    }
  )
})



################## my version to test the code
# dir.create( "data2" )
releaseDate <- "2015-08-21"
sapply( cohorts, function(element){
  tryCatch({
    downloadTCGA( cancerTypes = "OV", 
                  dataSet = "OV.Merge_rnaseqv2__illuminahiseq_rnaseqv2__unc_edu__Level_3__junction_quantification__data.Level_3.2015110100.0.0.tar.gz",
                  destDir = "C:/Users/Mcdade/Desktop/WorkflowEval", 
                  date = releaseDate )},
    error = function(cond){
      cat("Error: Maybe there weren't rnaseq data for ", element, " cancer.\n")
    }
  )
})



dir.create('data2')

# downloading rnaseq data
downloadTCGA( cancerTypes = 'BRCA',
              dataSet = 'Merge_rnaseqv2__illuminahiseq_rnaseqv2__unc_edu__Level_3__RSEM_genes_normalized__data.Level',
              destDir = 'data2' )

# shortening paths and directories
list.files( 'data2/') %>%
  file.path( 'data2', .) %>%
  file.rename( to = substr(.,start=1,stop=50))

# reading data
list.files( 'data2/') %>%
  file.path( 'data2', .) -> folder

folder %>%
  list.files %>%
  file.path( folder, .) %>%
  grep( pattern = 'illuminahiseq', x = ., value = TRUE) -> pathRNA
readTCGA( path = pathRNA, dataType = 'rnaseq' ) -> my_data


sapply( cohorts, function( element ){
  folder <- grep( paste0( "(_",element,"\\.", "|","_",element,"-FFPE)", ".*rnaseqv2"), 
                  list.files("data2"),value = TRUE)
  file <- grep( paste0(element, ".*rnaseqv2"), list.files( file.path( "data2",folder ) ),
                value = TRUE)
  path <- file.path( "data2", folder, file )
  assign( value = path, x = paste0(element, ".rnaseq.path"), envir = .GlobalEnv)
}) 


dir.create( 'hre')

downloadTCGA( cancerTypes = 'ACC', dataSet = 'miR_gene_expression',
              destDir = 'hre', date =  tail( checkTCGA('Dates'), 2 )[1],untarFile=TRUE, removeTar = TRUE )





