###############FormatExpressionData

####Prepare expression data for analysis

add_row_names_col<-function(ExpressionObject){
  formated_ExpressionObject<-cbind(row.names(ExpressionObject),ExpressionObject)
  names(formated_ExpressionObject)[1]<-paste("Symbol")
  return(formated_ExpressionObject)
}
