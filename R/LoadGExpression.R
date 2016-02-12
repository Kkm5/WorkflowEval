####first data loading function

#library("magrittr", lib.loc="~/R/win-library/3.2")
#library("RTCGA", lib.loc="~/R/win-library/3.2")
(cohorts <- infoTCGA() %>% 
  rownames() %>% 
  sub("-counts", "", x=.))
