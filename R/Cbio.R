rnaseqv2data.final<-t(rnaseqv2data.edit.t[row.names(rnaseqv2data.edit.t) %in% row.names(t(rppa.data.for.rnaseq)),])
rppa.data.for.rnaseq<-t(rppa.data.v2.t[row.names(rppa.data.v2.t) %in% row.names(t(rnaseqv2data.edit)),])
rnaseqv2data.edit.t<-t(rnaseqv2data.edit)
rppa.data.v2.t<-t(rppa.data.v2)
mycgds = CGDS("http://www.cbioportal.org/public-portal/")
library("cgdsr", lib.loc="C:/Users/kkm5/Documents/R/win-library/2.15")
getCancerStudies(mycgds)
mycancerstudy = getCancerStudies(mycgds)[49,1]
mycaselist = getCaseLists(mycgds,mycancerstudy)[9,1]
getGeneticProfiles(mycgds,mycancerstudy)[5,1]
