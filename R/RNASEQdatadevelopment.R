ovarian_RNAseqv1<-read.delim(file="ov_rnaseqv1_final.txt",header=FALSE)
ovarian_RNAseqv1$V1<-as.character(ovarian_RNAseqv1$V1)
ov_matches<-regexpr(pattern="^TCGA-[0-9]{2}-[0-9]{4}",ovarian_RNAseqv1$V1)
ov_text<-regmatches(ovarian_RNAseqv1$V1,ov_matches)
ovarian_RNAseqv1$TCGAcode<-regmatches(ovarian_RNAseqv1$V1,ov_matches)
ov_genename<-strsplit(ovarian_RNAseqv1$V1,split=":")
ov_genename.final<-unlist(ov_genename)[2*(1:length(ovarian_RNAseqv1$V1))]
ovarian_RNAseqv1$GeneName<-unlist(ov_genename_list)[2*(1:length(ovarian_RNAseqv1$V1))]
ov_name_matches<-regexpr(pattern="^[A-Z]{1,5}[0-9]{0,2}",ov_genename.final)
ovarian_RNAseqv1$GeneName<-regmatches(ov_genename.final,ov_name_matches)
ovarian_RNAseqv1_small<-ovarian_RNAseqv1[,4:6]
ov_name_matches<-regexpr(pattern="^[A-Z]{1,5}[0-9]{0,2}[A-Z]{0,3}[0-9]{0,1}",ov_genename.final)
ovarian_RNAseqv1$GeneName<-regmatches(ov_genename.final,ov_name_matches)
ovarian_RNAseqv1_small<-ovarian_RNAseqv1[,4:6]
ov_name_matches<-regexpr(pattern="^[A-Z]{1,5}[0-9]{0,2}[A-Z]{0,3}[0-9]{0,2}",ov_genename.final)
ovarian_RNAseqv1$GeneName<-regmatches(ov_genename.final,ov_name_matches)
ovarian_RNAseqv1_small<-ovarian_RNAseqv1[,4:6]
gsub(".","-",unique(ovarian_RNAseqv1_small$TCGAcode),fixed=TRUE)
ovarian_RNAseqv1_removed<-ovarian_RNAseqv1_fit[-as.numeric(rownames(ovarian_RNAseqv1_fit[ovarian_RNAseqv1_fit$TCGAcode=="TCGA.29.1770" | ovarian_RNAseqv1_fit$TCGAcode=="TCGA.29.1705" | ovarian_RNAseqv1_fit$TCGAcode=="TCGA.13.0913",])),]
ovarian_RNAseqv1_removed<-ovarian_RNAseqv1_fit[setdiff(rownames(ovarian_RNAseqv1_fit),index_missingTCGA),]
melted_ovarian_RNAseqv1_removed<-melt(ovarian_RNAseqv1_removed)
cast_ovarian<-cast(melted_ovarian_RNAseqv1_removed,TCGAcode~GeneName)
cast.ovarian.t<-t(cast_ovarian)
ovarian_rnaseq_v1_recasted<-as.data.frame(as.matrix(cast.ovarian.t))
rnaseqv2data.final.matchgenename<-rnaseqv2data.final[row.names(rnaseqv2data.final) %in% row.names(ovarian_rnaseq_v1_recasted),]
rnaseqv2data.final.matchgenename.t<-t(rnaseqv2data.final.matchgenename)
ovarian_rnaseq_v1_recasted.t<-t(ovarian_rnaseq_v1_recasted)
rnaseqv2data.removed.edited<-t(rnaseqv2data.final.matchgenename.t[row.names(rnaseqv2data.final.matchgenename.t) %in% row.names(ovarian_rnaseq_v1_recasted.t),])
idmap_rnaseqworkflow<-c(row.names(ovarian_rnaseq_v1_recasted),row.names(rnaseqv2data.removed.edited))
paste(row.names(ovarian_rnaseq_v1_recasted),"_v1",sep="")
row.names(rnaseqv2data.removed.edited)<-paste(row.names(rnaseqv2data.removed.edited),"_v2",sep="")
row.names(ovarian_rnaseq_v1_recasted)<-paste(row.names(ovarian_rnaseq_v1_recasted),"_v1",sep="")

rnaseqv2data.removed.edited.log2<-log2(rnaseqv2data.removed.edited)
ovarian_rnaseq_v1_recasted.log2<-log2(ovarian_rnaseq_v1_recasted)
RNASEQDATA<-rbind(ovarian_rnaseq_v1_recasted.log2,rnaseqv2data.removed.edited.log2)
idmap_rnaseqworkflow<-c(names(as.data.frame(ovarian_rnaseq_v1_recasted.t)),row.names(rnaseqv2data.removed.edited))

names(idmap_rnaseqworkflow)<-c("Symbol","RNAseqWorkflow")

rppa.edit<-rppa.data.for.rnaseq[row.names(rppa.data.for.rnaseq) %in% unique(idmap_rnaseqworkflow),]
rppa.edit<-as.data.frame(t(rppa.data.for.rnaseq[row.names(rppa.data.for.rnaseq) %in% unique(idmap_rnaseqworkflow),]))
                        
RPPADATA<-as.data.frame(t(rppa.edit[row.names(rppa.edit) %in% names(RNASEQDATA),]))
                        

                        