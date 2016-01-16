#' workflowpath.R
#' 
#'  Executes expected utility calculation for each path provided by workflowdesign
#'  @param 
#'   
#' 
workflowpath <-function(PairsDB,Lfp=1,Utp=2,deltaPlus=1,guarantee=1e-5,ruleColStart=14){
  x<-expectedUtility(label="Complete Set", dataset=PairsDB)
  y<-rbind(x,sapply(PairsDB[,14:23],expectedUtility,label=as.character(colnames(PairsDB[ruleColStart])),dataset = PairsDB %where% PairsDB[,ruleColStart]))
  return(y)
}