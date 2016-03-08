modify_data <-function(dataset, exp_colname="Symbol"){
  modified_dataset<-cbind(row.names(dataset),dataset)
  names(modified_dataset)[1]<-paste(exp_colname)
  return(modified_dataset)
}




assemble.data <- function(...) {
  input_list <- list(...)
  output_list <- lapply(X=input_list, function(x) {head(modify_data(x),n=2)})
  return(output_list)
}
