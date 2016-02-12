

############################setup the msms proteomic data set from the 98 Endometrial samples
############################ primary ids are determined from the IdMapBase function
primaryIDs <- IdMapBase$primaryIDs(msmsExperimentSet)

########################### secondary ids are from the affymetrix chip unlogged MAS 5.0 normalized
secondaryIDs <- IdMapBase$primaryIDs(mrnaExperimentSet);


########################## A data filter is applied to the msms data to remove all data points that
########################## average less than .52 spectral counts per Uniprot identifier,
########################## NA series removed as well
msmsExperimentSet <- DataFilter$do.apply(msmsExperimentSet,byRows=TRUE,filterFun=DataFilter$minAvgCountConstraint,filtParams=0.52,verbose=FALSE);
msmsExperimentSet <- DataFilter$removeNASeries(msmsExperimentSet,byRows=TRUE,verbose=FALSE);


######################### primary ids limited to ids that pass the filter
primaryIDs_from_dataset <- IdMapBase$primaryIDs(msmsExperimentSet);
secondaryIDs <- IdMapBase$primaryIDs(mrnaExperimentSet);

######################### joint ids determined from multiple annotation platforms
jointIdMap_from_dataset <- JointIdMap(examples$identDfList[c(6)], primaryIDs_from_dataset, verbose=FALSE);


######################### unique pairs from the jointIdMap
uniquePairs <- as.UniquePairs(getUnionIdMap(jointIdMap_from_dataset,verbose=FALSE),secondaryIDs);


######################### data for correlations are placed in a CorrData object
corrData <- CorrData(uniquePairs,msmsExperimentSet,mrnaExperimentSet,verbose=FALSE);


######################### Joint unique pairs
jointUniquePairs <- JointUniquePairs(uniquePairs,getIdMapList(jointIdMap_from_dataset,verbose=FALSE),verbose=FALSE);


######################## correlation values
corr <- Corr(corrData,method="pearson",verbose=FALSE);
corr.spearman <- Corr(corrData,method="spearman",verbose=FALSE);

corrSet <- getCorr(jointUniquePairs,corr,
                   groups=c("union","EnVision_Q","NetAffx_Q","DAVID_Q","DAVID_F"),
                   full.group=TRUE,verbose=FALSE);

########################## Not used here old method using mClust
mixture <- Mixture(corr,G=2,verbose=FALSE);
mixture.spearman <- Mixture(corr.spearman,G=2,verbose=FALSE);



bootstrap<-Bootstrap(corrData,Fisher=TRUE,verbose=TRUE);
#bootstrap$plot(new.plot=TRUE,file.copy=TRUE,copy.zoom=2,bg="white");

#fitLinear <- with(as.data.frame(jointUniquePairs),
#                  glm(qualityMeasure ~ DAVID_Q + EnVision_Q + NetAffx_Q,
#                      weights=(as.data.frame(bootstrap)$sd) ^ (-2))
#                  )
#coefficients(summary(fitLinear));

#head(jointUniquePairs$as.data.frame() )

pairs <- jointUniquePairs$as.data.frame()





##########################
pairMatchTable<-pairs
#bootModel.fulldata<-as.data.frame(bootstrap)
bootModel <-as.data.frame(bootstrap)
load("dotRandom.seed.rda")
Random.seed.save <- sync.Random.seed <- .Random.seed

#EMtest <- fit2clusters.v2(bootModel.fulldata$corr, bootModel.fulldata$sd^2,psi0Constraint=0, sameV=T,estimatesOnly=F, seed = sync.Random.seed)
#plot(EMtest[[1]]/(1+EMtest[[1]]), EMtest[[2]], main = "PosteriorProb vs PosterProbVar", xlab = "PosteriorProb", ylab = "PosteriorProbVar")
# plot(bootModel$corr, bootModel$sd^2)


#setcompare(bootModel.fulldata$Affy, PROBEFILTERDB.v3$Affy)

#postProbs<-as.vector(EMtest[[1]]/(1+EMtest[[1]]))  #not needed at the output is a datafrmae from fit2clusters
#postProbVar <-as.vector(EMtest[[2]])
#bootMergedWithPairs = merge(data.frame(postProbs=postProbs, postProbVar=postProbVar,bootModel.fulldata), pairMatchTable)



postProbs<- NULL
postProbVar<- NULL


EMtest<- fit2clusters.cancerinformatics(bootModel$corr, bootModel$sd^2,psi0Constraint=0, sameV=T,estimatesOnly=F, seed = Random.seed.save)
postProbs<-as.vector(EMtest[[1]]/(1+EMtest[[1]]))  #not needed at the output is a datafrmae from fit2clusters
postProbVar <-as.vector(EMtest[[2]])


bootMergedWithPairs = merge(data.frame(postProbs=postProbs, postProbVar=postProbVar,bootModel), pairMatchTable)

subsetForProbeFiltering = !is.na(match(bootModel$Affy, PROBEFILTERDB.v3$Affy))

#PDBA_accepts = PROBEFILTERDB.v3$Pass25PercentThreshold[match(bootModel$Affy, PROBEFILTERDB.v3$Affy)] == TRUE
PDBA_accepts = PROBEFILTERDB.v3$Match.pldb.proportion[match(bootModel$Affy, PROBEFILTERDB.v3$Affy)] > .30

Masker_accepts = PROBEFILTERDB.v3$MaskedByNCI[match(bootModel$Affy, PROBEFILTERDB.v3$Affy)] == FALSE
Complete_Masker_accepts = PROBEFILTERDB.v3$NCI.Acceptance.Measure[match(bootModel$Affy, PROBEFILTERDB.v3$Affy)] == TRUE
ENCODE_accepts = PROBEFILTERDB.v3$Gene.element.filter[match(bootModel$Affy, PROBEFILTERDB.v3$Affy)] == TRUE
Affy_tag_accepts = PROBEFILTERDB.v3$Affy.tag.filter[match(bootModel$Affy, PROBEFILTERDB.v3$Affy)] == TRUE
Affy_Grade_accepts = PROBEFILTERDB.v3$Affy.grade.filter[match(bootModel$Affy, PROBEFILTERDB.v3$Affy)] == TRUE
Geneannot_quality_accepts = PROBEFILTERDB.v3$Geneannot.quality.filter[match(bootModel$Affy, PROBEFILTERDB.v3$Affy)] == TRUE

#Geneannot_sensitivity_accepts = PROBEFILTERDB.v3$Geneannot.sensitivity.filter[match(bootModel$Affy, PROBEFILTERDB.v3$Affy)] == TRUE
Geneannot_sensitivity_accepts = PROBEFILTERDB.v3$Sensitivity[match(bootModel$Affy, PROBEFILTERDB.v3$Affy)] > .75

#Geneannot_specificity_accepts = PROBEFILTERDB.v3$Geneannot.specificity.filter[match(bootModel$Affy, PROBEFILTERDB.v3$Affy)] == TRUE
Geneannot_specificity_accepts = PROBEFILTERDB.v3$Specificity[match(bootModel$Affy, PROBEFILTERDB.v3$Affy)] > .90

Jetset_accepts = PROBEFILTERDB.v3$jetset.choice[match(bootModel$Affy, PROBEFILTERDB.v3$Affy)] == TRUE



#bootMergedWithPairs <- merge(bootModel,pairs, by="Affy")

#bootMergedWithPairs$PDBA_accepts = bootMergedWithPairs$Affy %in% (PROBEFILTERDB.v3$Affy[PROBEFILTERDB.v3$Pass25PercentThreshold==TRUE])
bootMergedWithPairs$PDBA_accepts = bootMergedWithPairs$Affy %in% (PROBEFILTERDB.v3$Affy[PROBEFILTERDB.v3$Match.pldb.proportion > 0.30])

bootMergedWithPairs$Masker_accepts = bootMergedWithPairs$Affy %in% (PROBEFILTERDB.v3$Affy[PROBEFILTERDB.v3$MaskedByNCI==FALSE])
bootMergedWithPairs$Complete_Masker_accepts = bootMergedWithPairs$Affy %in% (PROBEFILTERDB.v3$Affy[PROBEFILTERDB.v3$NCI.Acceptance.Measure==TRUE])
bootMergedWithPairs$ENCODE_accepts = bootMergedWithPairs$Affy %in% (PROBEFILTERDB.v3$Affy[PROBEFILTERDB.v3$Gene.element.filter==TRUE])
bootMergedWithPairs$Affy_tag_accepts = bootMergedWithPairs$Affy %in% (PROBEFILTERDB.v3$Affy[PROBEFILTERDB.v3$Affy.tag.filter==TRUE])
bootMergedWithPairs$Affy_Grade_accepts = bootMergedWithPairs$Affy %in% (PROBEFILTERDB.v3$Affy[PROBEFILTERDB.v3$Affy.grade.filter==TRUE])
bootMergedWithPairs$Geneannot_quality_accepts = bootMergedWithPairs$Affy %in% (PROBEFILTERDB.v3$Affy[PROBEFILTERDB.v3$Geneannot.quality.filter==TRUE])

#bootMergedWithPairs$Geneannot_sensitivity_accepts = bootMergedWithPairs$Affy %in% (PROBEFILTERDB.v3$Affy[PROBEFILTERDB.v3$Geneannot.sensitivity.filter==TRUE])
bootMergedWithPairs$Geneannot_sensitivity_accepts = bootMergedWithPairs$Affy %in% (PROBEFILTERDB.v3$Affy[PROBEFILTERDB.v3$Sensitivity > 0.75])

#bootMergedWithPairs$Geneannot_specificity_accepts = bootMergedWithPairs$Affy %in% (PROBEFILTERDB.v3$Affy[PROBEFILTERDB.v3$Geneannot.specificity.filter==TRUE])
bootMergedWithPairs$Geneannot_specificity_accepts = bootMergedWithPairs$Affy %in% (PROBEFILTERDB.v3$Affy[PROBEFILTERDB.v3$Specificity > 0.90])
bootMergedWithPairs$Jetset_accepts = bootMergedWithPairs$Affy %in% (PROBEFILTERDB.v3$Affy[PROBEFILTERDB.v3$jetset.choice==TRUE])




############JETSET GREEDY SEARCH

Lfp=1
Utp=2
deltaPlus=1
guarantee = 1e-5

resultAll = rbind(
  expectedUtility(label="Use All", dataset=bootMergedWithPairs),  # All (no filtering)
  ###level1
  expectedUtility(label="PD", dataset=bootMergedWithPairs %where% PDBA_accepts),   # PDBA accepts
  expectedUtility(label="M", dataset=bootMergedWithPairs %where% Masker_accepts),   # “Masker accepts”
  expectedUtility(label="AT", dataset=bootMergedWithPairs %where% Affy_tag_accepts),
  expectedUtility(label="AG", dataset=bootMergedWithPairs %where% Affy_Grade_accepts),
  expectedUtility(label="E", dataset=bootMergedWithPairs %where% ENCODE_accepts),
  expectedUtility(label="GQ", dataset=bootMergedWithPairs %where% Geneannot_quality_accepts),
  expectedUtility(label="GSEN", dataset=bootMergedWithPairs %where% Geneannot_sensitivity_accepts),
  expectedUtility(label="GSPE", dataset=bootMergedWithPairs %where% Geneannot_specificity_accepts),
  expectedUtility(label="J", dataset=bootMergedWithPairs %where% Jetset_accepts),
  ###Level2
  expectedUtility(label="+int(PD)", dataset=bootMergedWithPairs %where% (PDBA_accepts & Jetset_accepts)),
  expectedUtility(label="+int(M)", dataset=bootMergedWithPairs %where% (Masker_accepts & Jetset_accepts)),
  expectedUtility(label="+int(AT)", dataset=bootMergedWithPairs %where% (Affy_tag_accepts & Jetset_accepts)),
  expectedUtility(label="+int(AG)", dataset=bootMergedWithPairs %where% (Affy_Grade_accepts & Jetset_accepts)),
  expectedUtility(label="+int(E)", dataset=bootMergedWithPairs %where% (ENCODE_accepts & Jetset_accepts)),
  expectedUtility(label="+int(GQ)", dataset=bootMergedWithPairs %where% (Geneannot_quality_accepts & Jetset_accepts)),
  expectedUtility(label="+int(GSEN)", dataset=bootMergedWithPairs %where% (Geneannot_sensitivity_accepts & Jetset_accepts)),
  expectedUtility(label="+int(GSPE)", dataset=bootMergedWithPairs %where% (Geneannot_specificity_accepts & Jetset_accepts)),
  expectedUtility(label="+union(PD)", dataset=bootMergedWithPairs %where% (PDBA_accepts | Jetset_accepts)),
  expectedUtility(label="+union(M)", dataset=bootMergedWithPairs %where% (Masker_accepts | Jetset_accepts)),
  expectedUtility(label="+union(AT)", dataset=bootMergedWithPairs %where% (Affy_tag_accepts | Jetset_accepts)),
  expectedUtility(label="+union(AG)", dataset=bootMergedWithPairs %where% (Affy_Grade_accepts | Jetset_accepts)),
  expectedUtility(label="+union(E)", dataset=bootMergedWithPairs %where% (ENCODE_accepts | Jetset_accepts)),
  expectedUtility(label="+union(GQ)", dataset=bootMergedWithPairs %where% (Geneannot_quality_accepts | Jetset_accepts)),
  expectedUtility(label="+union(GSEN)", dataset=bootMergedWithPairs %where% (Geneannot_sensitivity_accepts | Jetset_accepts)),
  expectedUtility(label="+union(GSPE)", dataset=bootMergedWithPairs %where% (Geneannot_specificity_accepts | Jetset_accepts)),

  #######Level3
  expectedUtility(label="+int(M)", dataset=bootMergedWithPairs %where% (Masker_accepts & Jetset_accepts & PDBA_accepts)),
  expectedUtility(label="+int(AT)", dataset=bootMergedWithPairs %where% (Affy_tag_accepts & Jetset_accepts & PDBA_accepts)),
  expectedUtility(label="+int(AG)", dataset=bootMergedWithPairs %where% (Affy_Grade_accepts & Jetset_accepts & PDBA_accepts)),
  expectedUtility(label="+int(E)", dataset=bootMergedWithPairs %where% (ENCODE_accepts & Jetset_accepts & PDBA_accepts)),
  expectedUtility(label="+int(GQ)", dataset=bootMergedWithPairs %where% (Geneannot_quality_accepts & Jetset_accepts & PDBA_accepts)),
  expectedUtility(label="+int(GSEN)", dataset=bootMergedWithPairs %where% (Geneannot_sensitivity_accepts & Jetset_accepts & PDBA_accepts)),
  expectedUtility(label="+int(GSPE)", dataset=bootMergedWithPairs %where% (Geneannot_specificity_accepts & Jetset_accepts & PDBA_accepts)),
  expectedUtility(label="+union(PD)", dataset=bootMergedWithPairs %where% (PDBA_accepts | Jetset_accepts & PDBA_accepts)),
  expectedUtility(label="+union(M)", dataset=bootMergedWithPairs %where% (Masker_accepts | Jetset_accepts & PDBA_accepts)),
  expectedUtility(label="+union(AT)", dataset=bootMergedWithPairs %where% (Affy_tag_accepts | Jetset_accepts & PDBA_accepts)),
  expectedUtility(label="+union(AG)", dataset=bootMergedWithPairs %where% (Affy_Grade_accepts | Jetset_accepts & PDBA_accepts)),
  expectedUtility(label="+union(E)", dataset=bootMergedWithPairs %where% (ENCODE_accepts | Jetset_accepts & PDBA_accepts)),
  expectedUtility(label="+union(GQ)", dataset=bootMergedWithPairs %where% (Geneannot_quality_accepts | Jetset_accepts & PDBA_accepts)),
  expectedUtility(label="+union(GSEN)", dataset=bootMergedWithPairs %where% (Geneannot_sensitivity_accepts | Jetset_accepts & PDBA_accepts)),
  expectedUtility(label="+union(GSPE)", dataset=bootMergedWithPairs %where% (Geneannot_specificity_accepts | Jetset_accepts & PDBA_accepts)),


  #######Level4
  expectedUtility(label="+int(M)", dataset=bootMergedWithPairs %where% (Masker_accepts & Jetset_accepts & PDBA_accepts & Affy_Grade_accepts)),
  expectedUtility(label="+int(AT)", dataset=bootMergedWithPairs %where% (Affy_tag_accepts & Jetset_accepts & PDBA_accepts & Affy_Grade_accepts)),
  expectedUtility(label="+int(AG)", dataset=bootMergedWithPairs %where% (Affy_Grade_accepts & Jetset_accepts & PDBA_accepts & Affy_Grade_accepts)),
  expectedUtility(label="+int(E)", dataset=bootMergedWithPairs %where% (ENCODE_accepts & Jetset_accepts & PDBA_accepts & Affy_Grade_accepts)),
  expectedUtility(label="+int(GQ)", dataset=bootMergedWithPairs %where% (Geneannot_quality_accepts & Jetset_accepts & PDBA_accepts & Affy_Grade_accepts)),
  expectedUtility(label="+int(GSEN)", dataset=bootMergedWithPairs %where% (Geneannot_sensitivity_accepts & Jetset_accepts & PDBA_accepts & Affy_Grade_accepts)),
  expectedUtility(label="+int(GSPE)", dataset=bootMergedWithPairs %where% (Geneannot_specificity_accepts & Jetset_accepts & PDBA_accepts & Affy_Grade_accepts)),
  expectedUtility(label="+union(PD)", dataset=bootMergedWithPairs %where% (PDBA_accepts | Jetset_accepts & PDBA_accepts & Affy_Grade_accepts)),
  expectedUtility(label="+union(M)", dataset=bootMergedWithPairs %where% (Masker_accepts | Jetset_accepts & PDBA_accepts & Affy_Grade_accepts)),
  expectedUtility(label="+union(AT)", dataset=bootMergedWithPairs %where% (Affy_tag_accepts | Jetset_accepts & PDBA_accepts & Affy_Grade_accepts)),
  expectedUtility(label="+union(AG)", dataset=bootMergedWithPairs %where% (Affy_Grade_accepts | Jetset_accepts & PDBA_accepts & Affy_Grade_accepts)),
  expectedUtility(label="+union(E)", dataset=bootMergedWithPairs %where% (ENCODE_accepts | Jetset_accepts & PDBA_accepts & Affy_Grade_accepts)),
  expectedUtility(label="+union(GQ)", dataset=bootMergedWithPairs %where% (Geneannot_quality_accepts | Jetset_accepts & PDBA_accepts & Affy_Grade_accepts)),
  expectedUtility(label="+union(GSEN)", dataset=bootMergedWithPairs %where% (Geneannot_sensitivity_accepts | Jetset_accepts & PDBA_accepts & Affy_Grade_accepts)),
  expectedUtility(label="+union(GSPE)", dataset=bootMergedWithPairs %where% (Geneannot_specificity_accepts | Jetset_accepts & PDBA_accepts & Affy_Grade_accepts)),

  ###############Level5

  expectedUtility(label="+int(M)", dataset=bootMergedWithPairs %where% (Masker_accepts & Jetset_accepts & PDBA_accepts & Affy_Grade_accepts & Masker_accepts)),
  expectedUtility(label="+int(AT)", dataset=bootMergedWithPairs %where% (Affy_tag_accepts & Jetset_accepts & PDBA_accepts & Affy_Grade_accepts & Masker_accepts)),
  expectedUtility(label="+int(AG)", dataset=bootMergedWithPairs %where% (Affy_Grade_accepts & Jetset_accepts & PDBA_accepts & Affy_Grade_accepts & Masker_accepts)),
  expectedUtility(label="+int(E)", dataset=bootMergedWithPairs %where% (ENCODE_accepts & Jetset_accepts & PDBA_accepts & Affy_Grade_accepts & Masker_accepts)),
  expectedUtility(label="+int(GQ)", dataset=bootMergedWithPairs %where% (Geneannot_quality_accepts & Jetset_accepts & PDBA_accepts & Affy_Grade_accepts & Masker_accepts)),
  expectedUtility(label="+int(GSEN)", dataset=bootMergedWithPairs %where% (Geneannot_sensitivity_accepts & Jetset_accepts & PDBA_accepts & Affy_Grade_accepts & Masker_accepts)),
  expectedUtility(label="+int(GSPE)", dataset=bootMergedWithPairs %where% (Geneannot_specificity_accepts & Jetset_accepts & PDBA_accepts & Affy_Grade_accepts & Masker_accepts)),
  expectedUtility(label="+union(PD)", dataset=bootMergedWithPairs %where% (PDBA_accepts | Jetset_accepts & PDBA_accepts & Affy_Grade_accepts & Masker_accepts)),
  expectedUtility(label="+union(M)", dataset=bootMergedWithPairs %where% (Masker_accepts | Jetset_accepts & PDBA_accepts & Affy_Grade_accepts & Masker_accepts)),
  expectedUtility(label="+union(AT)", dataset=bootMergedWithPairs %where% (Affy_tag_accepts | Jetset_accepts & PDBA_accepts & Affy_Grade_accepts & Masker_accepts)),
  expectedUtility(label="+union(AG)", dataset=bootMergedWithPairs %where% (Affy_Grade_accepts | Jetset_accepts & PDBA_accepts & Affy_Grade_accepts & Masker_accepts)),
  expectedUtility(label="+union(E)", dataset=bootMergedWithPairs %where% (ENCODE_accepts | Jetset_accepts & PDBA_accepts & Affy_Grade_accepts & Masker_accepts)),
  expectedUtility(label="+union(GQ)", dataset=bootMergedWithPairs %where% (Geneannot_quality_accepts | Jetset_accepts & PDBA_accepts & Affy_Grade_accepts & Masker_accepts)),
  expectedUtility(label="+union(GSEN)", dataset=bootMergedWithPairs %where% (Geneannot_sensitivity_accepts | Jetset_accepts & PDBA_accepts & Affy_Grade_accepts & Masker_accepts)),
  expectedUtility(label="+union(GSPE)", dataset=bootMergedWithPairs %where% (Geneannot_specificity_accepts | Jetset_accepts & PDBA_accepts & Affy_Grade_accepts & Masker_accepts)),


  ################Level6
  expectedUtility(label="+int(M)", dataset=bootMergedWithPairs %where% (Masker_accepts & Jetset_accepts & PDBA_accepts & Affy_Grade_accepts & Masker_accepts & Geneannot_specificity_accepts)),
  expectedUtility(label="+int(AT)", dataset=bootMergedWithPairs %where% (Affy_tag_accepts & Jetset_accepts & PDBA_accepts & Affy_Grade_accepts & Masker_accepts & Geneannot_specificity_accepts)),
  expectedUtility(label="+int(AG)", dataset=bootMergedWithPairs %where% (Affy_Grade_accepts & Jetset_accepts & PDBA_accepts & Affy_Grade_accepts & Masker_accepts & Geneannot_specificity_accepts)),
  expectedUtility(label="+int(E)", dataset=bootMergedWithPairs %where% (ENCODE_accepts & Jetset_accepts & PDBA_accepts & Affy_Grade_accepts & Masker_accepts & Geneannot_specificity_accepts)),
  expectedUtility(label="+int(GQ)", dataset=bootMergedWithPairs %where% (Geneannot_quality_accepts & Jetset_accepts & PDBA_accepts & Affy_Grade_accepts & Masker_accepts & Geneannot_specificity_accepts)),
  expectedUtility(label="+int(GSEN)", dataset=bootMergedWithPairs %where% (Geneannot_sensitivity_accepts & Jetset_accepts & PDBA_accepts & Affy_Grade_accepts & Masker_accepts & Geneannot_specificity_accepts)),
  expectedUtility(label="+int(GSPE)", dataset=bootMergedWithPairs %where% (Geneannot_specificity_accepts & Jetset_accepts & PDBA_accepts & Affy_Grade_accepts & Masker_accepts & Geneannot_specificity_accepts)),
  expectedUtility(label="+union(PD)", dataset=bootMergedWithPairs %where% (PDBA_accepts | Jetset_accepts & PDBA_accepts & Affy_Grade_accepts & Masker_accepts & Geneannot_specificity_accepts)),
  expectedUtility(label="+union(M)", dataset=bootMergedWithPairs %where% (Masker_accepts | Jetset_accepts & PDBA_accepts & Affy_Grade_accepts & Masker_accepts & Geneannot_specificity_accepts)),
  expectedUtility(label="+union(AT)", dataset=bootMergedWithPairs %where% (Affy_tag_accepts | Jetset_accepts & PDBA_accepts & Affy_Grade_accepts & Masker_accepts & Geneannot_specificity_accepts)),
  expectedUtility(label="+union(AG)", dataset=bootMergedWithPairs %where% (Affy_Grade_accepts | Jetset_accepts & PDBA_accepts & Affy_Grade_accepts & Masker_accepts & Geneannot_specificity_accepts)),
  expectedUtility(label="+union(E)", dataset=bootMergedWithPairs %where% (ENCODE_accepts | Jetset_accepts & PDBA_accepts & Affy_Grade_accepts & Masker_accepts & Geneannot_specificity_accepts)),
  expectedUtility(label="+union(GQ)", dataset=bootMergedWithPairs %where% (Geneannot_quality_accepts | Jetset_accepts & PDBA_accepts & Affy_Grade_accepts & Masker_accepts & Geneannot_specificity_accepts)),
  expectedUtility(label="+union(GSEN)", dataset=bootMergedWithPairs %where% (Geneannot_sensitivity_accepts | Jetset_accepts & PDBA_accepts & Affy_Grade_accepts & Masker_accepts & Geneannot_specificity_accepts)),
  expectedUtility(label="+union(GSPE)", dataset=bootMergedWithPairs %where% (Geneannot_specificity_accepts | Jetset_accepts & PDBA_accepts & Affy_Grade_accepts & Masker_accepts & Geneannot_specificity_accepts))



  )


########################################################


resultJGreedy = rbind(
  expectedUtility(label="Use All", dataset=bootMergedWithPairs),  # All (no filtering)
  expectedUtility(label="J", dataset=bootMergedWithPairs %where% Jetset_accepts),
  expectedUtility(label="+int(PD)", dataset=bootMergedWithPairs %where% (PDBA_accepts & Jetset_accepts)),
  expectedUtility(label="+int(GSEN)", dataset=bootMergedWithPairs %where% (Jetset_accepts & PDBA_accepts & Geneannot_sensitivity_accepts)),
  expectedUtility(label="+int(GSPE)", dataset=bootMergedWithPairs %where% (Jetset_accepts & PDBA_accepts & Geneannot_sensitivity_accepts & Geneannot_specificity_accepts)),
  expectedUtility(label="+int(M)", dataset=bootMergedWithPairs %where% (Jetset_accepts & PDBA_accepts & Geneannot_sensitivity_accepts & Geneannot_specificity_accepts & Masker_accepts))
  )

plot(resultJGreedy$label,resultJGreedy$Eutility1, type="p",cex=.05, las=1,main="Greedy Search Summary", xlab="Probeset Resource Set", ylab="Mean Expected Utility")



##########################################################

resultLevel1 = rbind(
  expectedUtility(label="PD", dataset=bootMergedWithPairs %where% PDBA_accepts),   # PDBA accepts
  expectedUtility(label="M", dataset=bootMergedWithPairs %where% Masker_accepts),   # “Masker accepts”
  expectedUtility(label="AT", dataset=bootMergedWithPairs %where% Affy_tag_accepts),
  expectedUtility(label="AG", dataset=bootMergedWithPairs %where% Affy_Grade_accepts),
  expectedUtility(label="E", dataset=bootMergedWithPairs %where% ENCODE_accepts),
  expectedUtility(label="GQ", dataset=bootMergedWithPairs %where% Geneannot_quality_accepts),
  expectedUtility(label="GSEN", dataset=bootMergedWithPairs %where% Geneannot_sensitivity_accepts),
  expectedUtility(label="GSPE", dataset=bootMergedWithPairs %where% Geneannot_specificity_accepts),
  expectedUtility(label="J", dataset=bootMergedWithPairs %where% Jetset_accepts)
  )

###########################################################

resultLevel2 = rbind(
  expectedUtility(label="+int(PD)", dataset=bootMergedWithPairs %where% (PDBA_accepts & Jetset_accepts)),
  expectedUtility(label="+int(M)", dataset=bootMergedWithPairs %where% (Masker_accepts & Jetset_accepts)),
  expectedUtility(label="+int(AT)", dataset=bootMergedWithPairs %where% (Affy_tag_accepts & Jetset_accepts)),
  expectedUtility(label="+int(AG)", dataset=bootMergedWithPairs %where% (Affy_Grade_accepts & Jetset_accepts)),
  expectedUtility(label="+int(E)", dataset=bootMergedWithPairs %where% (ENCODE_accepts & Jetset_accepts)),
  expectedUtility(label="+int(GQ)", dataset=bootMergedWithPairs %where% (Geneannot_quality_accepts & Jetset_accepts)),
  expectedUtility(label="+int(GSEN)", dataset=bootMergedWithPairs %where% (Geneannot_sensitivity_accepts & Jetset_accepts)),
  expectedUtility(label="+int(GSPE)", dataset=bootMergedWithPairs %where% (Geneannot_specificity_accepts & Jetset_accepts)),
  expectedUtility(label="+union(PD)", dataset=bootMergedWithPairs %where% (PDBA_accepts | Jetset_accepts)),
  expectedUtility(label="+union(M)", dataset=bootMergedWithPairs %where% (Masker_accepts | Jetset_accepts)),
  expectedUtility(label="+union(AT)", dataset=bootMergedWithPairs %where% (Affy_tag_accepts | Jetset_accepts)),
  expectedUtility(label="+union(AG)", dataset=bootMergedWithPairs %where% (Affy_Grade_accepts | Jetset_accepts)),
  expectedUtility(label="+union(E)", dataset=bootMergedWithPairs %where% (ENCODE_accepts | Jetset_accepts)),
  expectedUtility(label="+union(GQ)", dataset=bootMergedWithPairs %where% (Geneannot_quality_accepts | Jetset_accepts)),
  expectedUtility(label="+union(GSEN)", dataset=bootMergedWithPairs %where% (Geneannot_sensitivity_accepts | Jetset_accepts)),
  expectedUtility(label="+union(GSPE)", dataset=bootMergedWithPairs %where% (Geneannot_specificity_accepts | Jetset_accepts))
  )


#################################
resultLevel3 =rbind(
  expectedUtility(label="+int(M)", dataset=bootMergedWithPairs %where% (Masker_accepts & Jetset_accepts & PDBA_accepts)),
  expectedUtility(label="+int(AT)", dataset=bootMergedWithPairs %where% (Affy_tag_accepts & Jetset_accepts & PDBA_accepts)),
  expectedUtility(label="+int(AG)", dataset=bootMergedWithPairs %where% (Affy_Grade_accepts & Jetset_accepts & PDBA_accepts)),
  expectedUtility(label="+int(E)", dataset=bootMergedWithPairs %where% (ENCODE_accepts & Jetset_accepts & PDBA_accepts)),
  expectedUtility(label="+int(GQ)", dataset=bootMergedWithPairs %where% (Geneannot_quality_accepts & Jetset_accepts & PDBA_accepts)),
  expectedUtility(label="+int(GSEN)", dataset=bootMergedWithPairs %where% (Geneannot_sensitivity_accepts & Jetset_accepts & PDBA_accepts)),
  expectedUtility(label="+int(GSPE)", dataset=bootMergedWithPairs %where% (Geneannot_specificity_accepts & Jetset_accepts & PDBA_accepts)),
  expectedUtility(label="+union(PD)", dataset=bootMergedWithPairs %where% (PDBA_accepts | Jetset_accepts & PDBA_accepts)),
  expectedUtility(label="+union(M)", dataset=bootMergedWithPairs %where% (Masker_accepts | Jetset_accepts & PDBA_accepts)),
  expectedUtility(label="+union(AT)", dataset=bootMergedWithPairs %where% (Affy_tag_accepts | Jetset_accepts & PDBA_accepts)),
  expectedUtility(label="+union(AG)", dataset=bootMergedWithPairs %where% (Affy_Grade_accepts | Jetset_accepts & PDBA_accepts)),
  expectedUtility(label="+union(E)", dataset=bootMergedWithPairs %where% (ENCODE_accepts | Jetset_accepts & PDBA_accepts)),
  expectedUtility(label="+union(GQ)", dataset=bootMergedWithPairs %where% (Geneannot_quality_accepts | Jetset_accepts & PDBA_accepts)),
  expectedUtility(label="+union(GSEN)", dataset=bootMergedWithPairs %where% (Geneannot_sensitivity_accepts | Jetset_accepts & PDBA_accepts)),
  expectedUtility(label="+union(GSPE)", dataset=bootMergedWithPairs %where% (Geneannot_specificity_accepts | Jetset_accepts & PDBA_accepts))

  )


##############################
resultLevel4 = rbind(
  expectedUtility(label="+int(M)", dataset=bootMergedWithPairs %where% (Masker_accepts & Jetset_accepts & PDBA_accepts & Geneannot_sensitivity_accepts)),
  expectedUtility(label="+int(AT)", dataset=bootMergedWithPairs %where% (Affy_tag_accepts & Jetset_accepts & PDBA_accepts & Geneannot_sensitivity_accepts)),
  expectedUtility(label="+int(AG)", dataset=bootMergedWithPairs %where% (Affy_Grade_accepts & Jetset_accepts & PDBA_accepts & Geneannot_sensitivity_accepts)),
  expectedUtility(label="+int(E)", dataset=bootMergedWithPairs %where% (ENCODE_accepts & Jetset_accepts & PDBA_accepts & Geneannot_sensitivity_accepts)),
  expectedUtility(label="+int(GQ)", dataset=bootMergedWithPairs %where% (Geneannot_quality_accepts & Jetset_accepts & PDBA_accepts & Geneannot_sensitivity_accepts)),
  expectedUtility(label="+int(GSEN)", dataset=bootMergedWithPairs %where% (Geneannot_sensitivity_accepts & Jetset_accepts & PDBA_accepts & Geneannot_sensitivity_accepts)),
  expectedUtility(label="+int(GSPE)", dataset=bootMergedWithPairs %where% (Geneannot_specificity_accepts & Jetset_accepts & PDBA_accepts & Geneannot_sensitivity_accepts)),
  expectedUtility(label="+union(PD)", dataset=bootMergedWithPairs %where% (PDBA_accepts | Jetset_accepts & PDBA_accepts & Geneannot_sensitivity_accepts)),
  expectedUtility(label="+union(M)", dataset=bootMergedWithPairs %where% (Masker_accepts | Jetset_accepts & PDBA_accepts & Geneannot_sensitivity_accepts)),
  expectedUtility(label="+union(AT)", dataset=bootMergedWithPairs %where% (Affy_tag_accepts | Jetset_accepts & PDBA_accepts & Geneannot_sensitivity_accepts)),
  expectedUtility(label="+union(AG)", dataset=bootMergedWithPairs %where% (Affy_Grade_accepts | Jetset_accepts & PDBA_accepts & Geneannot_sensitivity_accepts)),
  expectedUtility(label="+union(E)", dataset=bootMergedWithPairs %where% (ENCODE_accepts | Jetset_accepts & PDBA_accepts & Geneannot_sensitivity_accepts)),
  expectedUtility(label="+union(GQ)", dataset=bootMergedWithPairs %where% (Geneannot_quality_accepts | Jetset_accepts & PDBA_accepts & Geneannot_sensitivity_accepts)),
  expectedUtility(label="+union(GSEN)", dataset=bootMergedWithPairs %where% (Geneannot_sensitivity_accepts | Jetset_accepts & PDBA_accepts & Geneannot_sensitivity_accepts)),
  expectedUtility(label="+union(GSPE)", dataset=bootMergedWithPairs %where% (Geneannot_specificity_accepts | Jetset_accepts & PDBA_accepts & Geneannot_sensitivity_accepts))
  )

###########################


resultLevel5 = rbind(
  expectedUtility(label="+int(M)", dataset=bootMergedWithPairs %where% (Masker_accepts & Jetset_accepts & PDBA_accepts & Geneannot_sensitivity_accepts & Geneannot_specificity_accepts)),
  expectedUtility(label="+int(AT)", dataset=bootMergedWithPairs %where% (Affy_tag_accepts & Jetset_accepts & PDBA_accepts & Geneannot_sensitivity_accepts & Geneannot_specificity_accepts)),
  expectedUtility(label="+int(AG)", dataset=bootMergedWithPairs %where% (Affy_Grade_accepts & Jetset_accepts & PDBA_accepts & Geneannot_sensitivity_accepts & Geneannot_specificity_accepts)),
  expectedUtility(label="+int(E)", dataset=bootMergedWithPairs %where% (ENCODE_accepts & Jetset_accepts & PDBA_accepts & Geneannot_sensitivity_accepts & Geneannot_specificity_accepts)),
  expectedUtility(label="+int(GQ)", dataset=bootMergedWithPairs %where% (Geneannot_quality_accepts & Jetset_accepts & PDBA_accepts & Geneannot_sensitivity_accepts & Geneannot_specificity_accepts)),
  expectedUtility(label="+int(GSEN)", dataset=bootMergedWithPairs %where% (Geneannot_sensitivity_accepts & Jetset_accepts & PDBA_accepts & Geneannot_sensitivity_accepts & Geneannot_specificity_accepts)),
  expectedUtility(label="+int(GSPE)", dataset=bootMergedWithPairs %where% (Geneannot_specificity_accepts & Jetset_accepts & PDBA_accepts & Geneannot_sensitivity_accepts & Geneannot_specificity_accepts)),
  expectedUtility(label="+union(PD)", dataset=bootMergedWithPairs %where% (PDBA_accepts | Jetset_accepts & PDBA_accepts & Geneannot_sensitivity_accepts & Geneannot_specificity_accepts)),
  expectedUtility(label="+union(M)", dataset=bootMergedWithPairs %where% (Masker_accepts | Jetset_accepts & PDBA_accepts & Geneannot_sensitivity_accepts & Geneannot_specificity_accepts)),
  expectedUtility(label="+union(AT)", dataset=bootMergedWithPairs %where% (Affy_tag_accepts | Jetset_accepts & PDBA_accepts & Geneannot_sensitivity_accepts & Geneannot_specificity_accepts)),
  expectedUtility(label="+union(AG)", dataset=bootMergedWithPairs %where% (Affy_Grade_accepts | Jetset_accepts & PDBA_accepts & Geneannot_sensitivity_accepts & Geneannot_specificity_accepts)),
  expectedUtility(label="+union(E)", dataset=bootMergedWithPairs %where% (ENCODE_accepts | Jetset_accepts & PDBA_accepts & Geneannot_sensitivity_accepts & Geneannot_specificity_accepts)),
  expectedUtility(label="+union(GQ)", dataset=bootMergedWithPairs %where% (Geneannot_quality_accepts | Jetset_accepts & PDBA_accepts & Geneannot_sensitivity_accepts & Geneannot_specificity_accepts)),
  expectedUtility(label="+union(GSEN)", dataset=bootMergedWithPairs %where% (Geneannot_sensitivity_accepts | Jetset_accepts & PDBA_accepts & Geneannot_sensitivity_accepts & Geneannot_specificity_accepts)),
  expectedUtility(label="+union(GSPE)", dataset=bootMergedWithPairs %where% (Geneannot_specificity_accepts | Jetset_accepts & PDBA_accepts & Geneannot_sensitivity_accepts & Geneannot_specificity_accepts))

  )
Geneannot_sensitivity_accepts & Geneannot_specificity_accepts & Masker_accepts

resultLevel6 = rbind(
  expectedUtility(label="+int(M)", dataset=bootMergedWithPairs %where% (Masker_accepts & Jetset_accepts & PDBA_accepts & Geneannot_sensitivity_accepts & Geneannot_specificity_accepts & Masker_accepts)),
  expectedUtility(label="+int(AT)", dataset=bootMergedWithPairs %where% (Affy_tag_accepts & Jetset_accepts & PDBA_accepts & Geneannot_sensitivity_accepts & Geneannot_specificity_accepts & Masker_accepts)),
  expectedUtility(label="+int(AG)", dataset=bootMergedWithPairs %where% (Affy_Grade_accepts & Jetset_accepts & PDBA_accepts & Geneannot_sensitivity_accepts & Geneannot_specificity_accepts & Masker_accepts)),
  expectedUtility(label="+int(E)", dataset=bootMergedWithPairs %where% (ENCODE_accepts & Jetset_accepts & PDBA_accepts & Geneannot_sensitivity_accepts & Geneannot_specificity_accepts & Masker_accepts)),
  expectedUtility(label="+int(GQ)", dataset=bootMergedWithPairs %where% (Geneannot_quality_accepts & Jetset_accepts & PDBA_accepts & Geneannot_sensitivity_accepts & Geneannot_specificity_accepts & Masker_accepts)),
  expectedUtility(label="+int(GSEN)", dataset=bootMergedWithPairs %where% (Geneannot_sensitivity_accepts & Jetset_accepts & PDBA_accepts & Geneannot_sensitivity_accepts & Geneannot_specificity_accepts & Masker_accepts)),
  expectedUtility(label="+int(GSPE)", dataset=bootMergedWithPairs %where% (Geneannot_specificity_accepts & Jetset_accepts & PDBA_accepts & Geneannot_sensitivity_accepts & Geneannot_specificity_accepts & Masker_accepts)),
  expectedUtility(label="+union(PD)", dataset=bootMergedWithPairs %where% (PDBA_accepts | Jetset_accepts & PDBA_accepts & Geneannot_sensitivity_accepts & Geneannot_specificity_accepts & Masker_accepts)),
  expectedUtility(label="+union(M)", dataset=bootMergedWithPairs %where% (Masker_accepts | Jetset_accepts & PDBA_accepts & Geneannot_sensitivity_accepts & Geneannot_specificity_accepts & Masker_accepts)),
  expectedUtility(label="+union(AT)", dataset=bootMergedWithPairs %where% (Affy_tag_accepts | Jetset_accepts & PDBA_accepts & Geneannot_sensitivity_accepts & Geneannot_specificity_accepts & Masker_accepts)),
  expectedUtility(label="+union(AG)", dataset=bootMergedWithPairs %where% (Affy_Grade_accepts | Jetset_accepts & PDBA_accepts & Geneannot_sensitivity_accepts & Geneannot_specificity_accepts & Masker_accepts)),
  expectedUtility(label="+union(E)", dataset=bootMergedWithPairs %where% (ENCODE_accepts | Jetset_accepts & PDBA_accepts & Geneannot_sensitivity_accepts & Geneannot_specificity_accepts & Masker_accepts)),
  expectedUtility(label="+union(GQ)", dataset=bootMergedWithPairs %where% (Geneannot_quality_accepts | Jetset_accepts & PDBA_accepts & Geneannot_sensitivity_accepts & Geneannot_specificity_accepts & Masker_accepts)),
  expectedUtility(label="+union(GSEN)", dataset=bootMergedWithPairs %where% (Geneannot_sensitivity_accepts | Jetset_accepts & PDBA_accepts & Geneannot_sensitivity_accepts & Geneannot_specificity_accepts & Masker_accepts)),
  expectedUtility(label="+union(GSPE)", dataset=bootMergedWithPairs %where% (Geneannot_specificity_accepts | Jetset_accepts & PDBA_accepts & Geneannot_sensitivity_accepts & Geneannot_specificity_accepts & Masker_accepts))




  )

plot(resultAll$PrPlus, resultAll$nPairs, pch = "", main = "Greedy search across IDF resources by mean expected utility", xlab = "Estimated Proportion of Concordance", ylab = "Number of Pairs")
lines(resultJGreedy$PrPlus, resultJGreedy$nPairs, col ="black")
points(resultLevel1$PrPlus, resultLevel1$nPairs, col="red",pch = "1")
points(resultLevel2$PrPlus, resultLevel2$nPairs, col="blue",pch = "2")
points(resultLevel3$PrPlus, resultLevel3$nPairs, col="green", pch = "3")
points(resultLevel4$PrPlus, resultLevel4$nPairs, col="orange", pch = "4")
points(resultLevel5$PrPlus, resultLevel5$nPairs, col="purple", pch = "5")
points(resultLevel6$PrPlus, resultLevel6$nPairs, col="pink", pch = "6")
text(resultLevel1$PrPlus, resultLevel1$nPairs, labels = resultLevel1$label, pos = 1)

###use pch to legend
### chaange color on 1,2,3

rownames(resultJGreedy) = NULL
options(digits=3)
utilitytableJGreedy<- resultJGreedy[order(resultJGreedy$Eutility),]
utilitytable.EutilitymeanJGreedy<- resultJGreedy[order(resultJGreedy$Eutility1),]
utilitytable.EutilitymeanJGreedy








rownames(result) = NULL
options(digits=3)
utilitytable<- result[order(result$Eutility),]
utilitytable.Eutilitymean<- result[order(result$Eutility1),]
utilitytable.Eutilitymean
names(utilitytable)

with(PROBEFILTERDB.v3, table(jetset.choice, Pass25PercentThreshold))
with(PROBEFILTERDB.v3, table(jetset.choice, Affy.tag.filter))
with(PROBEFILTERDB.v3, table(Pass25PercentThreshold, Affy.tag.filter))

fisher.test(with(PROBEFILTERDB.v3, table(jetset.choice, Pass25PercentThreshold)))
fisher.test(with(PROBEFILTERDB.v3, table(jetset.choice, Affy.tag.filter)))
fisher.test(with(PROBEFILTERDB.v3, table(Pass25PercentThreshold, Affy.tag.filter)))


PROBEFILTER.envisononly<-PROBEFILTERDB.v3[PROBEFILTERDB.v3$Affy %in% as.data.frame(uniquePairs)$Affy,]

with(PROBEFILTER.envisononly, table(jetset.choice, Pass25PercentThreshold))
with(PROBEFILTER.envisononly, table(jetset.choice, Affy.tag.filter))
with(PROBEFILTER.envisononly, table(Pass25PercentThreshold, Affy.tag.filter))

fisher.test(with(PROBEFILTERDB.v3, table(jetset.choice, Gene.element.filter)))$estimate

fisher.test(with(PROBEFILTERDB.v3, table(Gene.element.filter, Pass25PercentThreshold)))$estimate

fisher.test(with(PROBEFILTERDB.v3, table(Pass25PercentThreshold, Geneannot.sensitivity.filter)))$estimate

fisher.test(with(PROBEFILTERDB.v3, table(Geneannot.sensitivity.filter, Geneannot.specificity.filter)))$estimate

fisher.test(with(PROBEFILTERDB.v3, table(Geneannot.specificity.filter, Geneannot.quality.filter)))$estimate

fisher.test(with(PROBEFILTERDB.v3, table(Geneannot.quality.filter, Affy.tag.filter)))$estimate

fisher.test(with(PROBEFILTERDB.v3, table(Affy.tag.filter, Affy.grade.filter)))$estimate

fisher.test(with(PROBEFILTERDB.v3, table(Affy.grade.filter, MaskedByNCI)))$estimate


table(PROBEFILTERDB.v3$Geneannot.sensitivity.filter,PROBEFILTERDB.v3$Geneannot.quality.filter)
#######################################################

fisher.test(table(Jetset_accepts,ENCODE_accepts))$estimate

fisher.test(table(ENCODE_accepts,PDBA_accepts))$estimate

fisher.test(table(PDBA_accepts,Geneannot_sensitivity_accepts))$estimate

fisher.test(table(Geneannot_sensitivity_accepts,Geneannot_specificity_accepts))$estimate

fisher.test(table(Geneannot_specificity_accepts,Geneannot_quality_accepts))$estimate

fisher.test(table(Geneannot_quality_accepts,Affy_tag_accepts))$estimate

fisher.test(table(Affy_tag_accepts,Affy_Grade_accepts))$estimate

fisher.test(table(Affy_Grade_accepts,Masker_accepts))$estimate



with(bootMergedWithPairs[bootMergedWithPairs$Affy %in% PROBEFILTER.envisononly$Affy,], points(postProbs ~ jitter(as.numeric(Jetset_accepts+1), amount=.1)))
with(bootMergedWithPairs[bootMergedWithPairs$Affy %in% PROBEFILTER.envisononly$Affy,], boxplot(postProbs ~ Jetset_accepts))






################Threshold selection

PDBA_accepts0 = PROBEFILTERDB.v3$Match.pldb.proportion[match(bootModel$Affy, PROBEFILTERDB.v3$Affy)] > 0
PDBA_accepts0.05 = PROBEFILTERDB.v3$Match.pldb.proportion[match(bootModel$Affy, PROBEFILTERDB.v3$Affy)] > .05
PDBA_accepts0.10 = PROBEFILTERDB.v3$Match.pldb.proportion[match(bootModel$Affy, PROBEFILTERDB.v3$Affy)] > .1
PDBA_accepts0.15 = PROBEFILTERDB.v3$Match.pldb.proportion[match(bootModel$Affy, PROBEFILTERDB.v3$Affy)] > .15
PDBA_accepts0.20 = PROBEFILTERDB.v3$Match.pldb.proportion[match(bootModel$Affy, PROBEFILTERDB.v3$Affy)] > .2
PDBA_accepts0.25 = PROBEFILTERDB.v3$Match.pldb.proportion[match(bootModel$Affy, PROBEFILTERDB.v3$Affy)] > .25
PDBA_accepts0.30 = PROBEFILTERDB.v3$Match.pldb.proportion[match(bootModel$Affy, PROBEFILTERDB.v3$Affy)] > .30
PDBA_accepts0.35 = PROBEFILTERDB.v3$Match.pldb.proportion[match(bootModel$Affy, PROBEFILTERDB.v3$Affy)] > .35
PDBA_accepts0.40 = PROBEFILTERDB.v3$Match.pldb.proportion[match(bootModel$Affy, PROBEFILTERDB.v3$Affy)] > .40
PDBA_accepts0.45 = PROBEFILTERDB.v3$Match.pldb.proportion[match(bootModel$Affy, PROBEFILTERDB.v3$Affy)] > .45
PDBA_accepts0.50 = PROBEFILTERDB.v3$Match.pldb.proportion[match(bootModel$Affy, PROBEFILTERDB.v3$Affy)] > .5
PDBA_accepts0.55 = PROBEFILTERDB.v3$Match.pldb.proportion[match(bootModel$Affy, PROBEFILTERDB.v3$Affy)] > .55
PDBA_accepts0.60 = PROBEFILTERDB.v3$Match.pldb.proportion[match(bootModel$Affy, PROBEFILTERDB.v3$Affy)] > .60
PDBA_accepts0.65 = PROBEFILTERDB.v3$Match.pldb.proportion[match(bootModel$Affy, PROBEFILTERDB.v3$Affy)] > .65
PDBA_accepts0.70 = PROBEFILTERDB.v3$Match.pldb.proportion[match(bootModel$Affy, PROBEFILTERDB.v3$Affy)] > .70
PDBA_accepts0.75 = PROBEFILTERDB.v3$Match.pldb.proportion[match(bootModel$Affy, PROBEFILTERDB.v3$Affy)] > .75
PDBA_accepts0.80 = PROBEFILTERDB.v3$Match.pldb.proportion[match(bootModel$Affy, PROBEFILTERDB.v3$Affy)] > .80
PDBA_accepts0.85 = PROBEFILTERDB.v3$Match.pldb.proportion[match(bootModel$Affy, PROBEFILTERDB.v3$Affy)] > .85
PDBA_accepts0.90 = PROBEFILTERDB.v3$Match.pldb.proportion[match(bootModel$Affy, PROBEFILTERDB.v3$Affy)] > .90
PDBA_accepts0.95 = PROBEFILTERDB.v3$Match.pldb.proportion[match(bootModel$Affy, PROBEFILTERDB.v3$Affy)] > .95
PDBA_accepts1 = PROBEFILTERDB.v3$Match.pldb.proportion[match(bootModel$Affy, PROBEFILTERDB.v3$Affy)] == 1




#bootMergedWithPairs <- merge(bootModel,pairs, by="Affy")

bootMergedWithPairs$PDBA_accepts0 = bootMergedWithPairs$Affy %in% (PROBEFILTERDB.v3$Affy[PROBEFILTERDB.v3$Match.pldb.proportion > 0])
bootMergedWithPairs$PDBA_accepts0.05 = bootMergedWithPairs$Affy %in% (PROBEFILTERDB.v3$Affy[PROBEFILTERDB.v3$Match.pldb.proportion > 0.05])
bootMergedWithPairs$PDBA_accepts0.10 = bootMergedWithPairs$Affy %in% (PROBEFILTERDB.v3$Affy[PROBEFILTERDB.v3$Match.pldb.proportion > 0.10])
bootMergedWithPairs$PDBA_accepts0.15 = bootMergedWithPairs$Affy %in% (PROBEFILTERDB.v3$Affy[PROBEFILTERDB.v3$Match.pldb.proportion > 0.15])
bootMergedWithPairs$PDBA_accepts0.20 = bootMergedWithPairs$Affy %in% (PROBEFILTERDB.v3$Affy[PROBEFILTERDB.v3$Match.pldb.proportion > 0.20])
bootMergedWithPairs$PDBA_accepts0.25 = bootMergedWithPairs$Affy %in% (PROBEFILTERDB.v3$Affy[PROBEFILTERDB.v3$Match.pldb.proportion > 0.25])
bootMergedWithPairs$PDBA_accepts0.30 = bootMergedWithPairs$Affy %in% (PROBEFILTERDB.v3$Affy[PROBEFILTERDB.v3$Match.pldb.proportion > 0.30])
bootMergedWithPairs$PDBA_accepts0.35 = bootMergedWithPairs$Affy %in% (PROBEFILTERDB.v3$Affy[PROBEFILTERDB.v3$Match.pldb.proportion > 0.35])
bootMergedWithPairs$PDBA_accepts0.40 = bootMergedWithPairs$Affy %in% (PROBEFILTERDB.v3$Affy[PROBEFILTERDB.v3$Match.pldb.proportion > 0.40])
bootMergedWithPairs$PDBA_accepts0.45 = bootMergedWithPairs$Affy %in% (PROBEFILTERDB.v3$Affy[PROBEFILTERDB.v3$Match.pldb.proportion > 0.45])
bootMergedWithPairs$PDBA_accepts0.50 = bootMergedWithPairs$Affy %in% (PROBEFILTERDB.v3$Affy[PROBEFILTERDB.v3$Match.pldb.proportion > 0.50])
bootMergedWithPairs$PDBA_accepts0.55 = bootMergedWithPairs$Affy %in% (PROBEFILTERDB.v3$Affy[PROBEFILTERDB.v3$Match.pldb.proportion > 0.55])
bootMergedWithPairs$PDBA_accepts0.60 = bootMergedWithPairs$Affy %in% (PROBEFILTERDB.v3$Affy[PROBEFILTERDB.v3$Match.pldb.proportion > 0.60])
bootMergedWithPairs$PDBA_accepts0.65 = bootMergedWithPairs$Affy %in% (PROBEFILTERDB.v3$Affy[PROBEFILTERDB.v3$Match.pldb.proportion > 0.65])
bootMergedWithPairs$PDBA_accepts0.70 = bootMergedWithPairs$Affy %in% (PROBEFILTERDB.v3$Affy[PROBEFILTERDB.v3$Match.pldb.proportion > 0.70])
bootMergedWithPairs$PDBA_accepts0.75 = bootMergedWithPairs$Affy %in% (PROBEFILTERDB.v3$Affy[PROBEFILTERDB.v3$Match.pldb.proportion > 0.75])
bootMergedWithPairs$PDBA_accepts0.80 = bootMergedWithPairs$Affy %in% (PROBEFILTERDB.v3$Affy[PROBEFILTERDB.v3$Match.pldb.proportion > 0.80])
bootMergedWithPairs$PDBA_accepts0.85 = bootMergedWithPairs$Affy %in% (PROBEFILTERDB.v3$Affy[PROBEFILTERDB.v3$Match.pldb.proportion > 0.85])
bootMergedWithPairs$PDBA_accepts0.90 = bootMergedWithPairs$Affy %in% (PROBEFILTERDB.v3$Affy[PROBEFILTERDB.v3$Match.pldb.proportion > 0.90])
bootMergedWithPairs$PDBA_accepts0.95 = bootMergedWithPairs$Affy %in% (PROBEFILTERDB.v3$Affy[PROBEFILTERDB.v3$Match.pldb.proportion > 0.95])
bootMergedWithPairs$PDBA_accepts1 = bootMergedWithPairs$Affy %in% (PROBEFILTERDB.v3$Affy[PROBEFILTERDB.v3$Match.pldb.proportion == 1])






table(PDBA_accepts, Masker_accepts)
table(bootMergedWithPairs$PDBA_accepts, bootMergedWithPairs$Masker_accepts)
table(Geneannot_specificity_accepts, Geneannot_sensitivity_accepts)
table(bootMergedWithPairs$Geneannot_specificity_accepts, bootMergedWithPairs$Geneannot_sensitivity_accepts)


Lfp=1
Utp=2
deltaPlus=1/2
guarantee = 1e-5

resultPDBAQC = rbind(
  expectedUtility(label="0", dataset=bootMergedWithPairs %where% PDBA_accepts0),   # PDBA accepts
  expectedUtility(label=".05", dataset=bootMergedWithPairs %where% PDBA_accepts0.05),
  expectedUtility(label=".10", dataset=bootMergedWithPairs %where% PDBA_accepts0.10),
  expectedUtility(label=".15", dataset=bootMergedWithPairs %where% PDBA_accepts0.15),
  expectedUtility(label=".20", dataset=bootMergedWithPairs %where% PDBA_accepts0.20),
  expectedUtility(label=".25", dataset=bootMergedWithPairs %where% PDBA_accepts0.25),
  expectedUtility(label=".30", dataset=bootMergedWithPairs %where% PDBA_accepts0.30),
  expectedUtility(label=".35", dataset=bootMergedWithPairs %where% PDBA_accepts0.35),
  expectedUtility(label=".40", dataset=bootMergedWithPairs %where% PDBA_accepts0.40),
  expectedUtility(label=".45", dataset=bootMergedWithPairs %where% PDBA_accepts0.45),
  expectedUtility(label=".50", dataset=bootMergedWithPairs %where% PDBA_accepts0.50),
  expectedUtility(label=".55", dataset=bootMergedWithPairs %where% PDBA_accepts0.55),
  expectedUtility(label=".60", dataset=bootMergedWithPairs %where% PDBA_accepts0.60),
  expectedUtility(label=".65", dataset=bootMergedWithPairs %where% PDBA_accepts0.65),
  expectedUtility(label=".70", dataset=bootMergedWithPairs %where% PDBA_accepts0.70),
  expectedUtility(label=".75", dataset=bootMergedWithPairs %where% PDBA_accepts0.75),
  expectedUtility(label=".80", dataset=bootMergedWithPairs %where% PDBA_accepts0.80),
  expectedUtility(label=".85", dataset=bootMergedWithPairs %where% PDBA_accepts0.85),
  expectedUtility(label=".90", dataset=bootMergedWithPairs %where% PDBA_accepts0.90),
  expectedUtility(label=".95", dataset=bootMergedWithPairs %where% PDBA_accepts0.95),
  expectedUtility(label="1", dataset=bootMergedWithPairs %where% PDBA_accepts1)
  )




rownames(resultPDBAQC) = NULL
options(digits=3)
utilitytablePDBAQC<- resultPDBAQC[order(resultPDBAQC$Eutility),]
utilitytablePDBAQC.Eutilitymean<- resultPDBAQC[order(resultPDBAQC$Eutility1),]
utilitytablePDBAQC


plot(resultPDBAQC$label,resultPDBAQC$Eutility,type="p",cex=.1, las=1,main="Optimize Plandbaffy", xlab="Proportion Threshold", ylab="Total Expected Utility")


GSEN_accepts0 = PROBEFILTERDB.v3$Sensitivity[match(bootModel$Affy, PROBEFILTERDB.v3$Affy)] > 0
GSEN_accepts0.05 = PROBEFILTERDB.v3$Sensitivity[match(bootModel$Affy, PROBEFILTERDB.v3$Affy)] > .05
GSEN_accepts0.10 = PROBEFILTERDB.v3$Sensitivity[match(bootModel$Affy, PROBEFILTERDB.v3$Affy)] > .1
GSEN_accepts0.15 = PROBEFILTERDB.v3$Sensitivity[match(bootModel$Affy, PROBEFILTERDB.v3$Affy)] > .15
GSEN_accepts0.20 = PROBEFILTERDB.v3$Sensitivity[match(bootModel$Affy, PROBEFILTERDB.v3$Affy)] > .2
GSEN_accepts0.25 = PROBEFILTERDB.v3$Sensitivity[match(bootModel$Affy, PROBEFILTERDB.v3$Affy)] > .25
GSEN_accepts0.30 = PROBEFILTERDB.v3$Sensitivity[match(bootModel$Affy, PROBEFILTERDB.v3$Affy)] > .30
GSEN_accepts0.35 = PROBEFILTERDB.v3$Sensitivity[match(bootModel$Affy, PROBEFILTERDB.v3$Affy)] > .35
GSEN_accepts0.40 = PROBEFILTERDB.v3$Sensitivity[match(bootModel$Affy, PROBEFILTERDB.v3$Affy)] > .40
GSEN_accepts0.45 = PROBEFILTERDB.v3$Sensitivity[match(bootModel$Affy, PROBEFILTERDB.v3$Affy)] > .45
GSEN_accepts0.50 = PROBEFILTERDB.v3$Sensitivity[match(bootModel$Affy, PROBEFILTERDB.v3$Affy)] > .5
GSEN_accepts0.55 = PROBEFILTERDB.v3$Sensitivity[match(bootModel$Affy, PROBEFILTERDB.v3$Affy)] > .55
GSEN_accepts0.60 = PROBEFILTERDB.v3$Sensitivity[match(bootModel$Affy, PROBEFILTERDB.v3$Affy)] > .60
GSEN_accepts0.65 = PROBEFILTERDB.v3$Sensitivity[match(bootModel$Affy, PROBEFILTERDB.v3$Affy)] > .65
GSEN_accepts0.70 = PROBEFILTERDB.v3$Sensitivity[match(bootModel$Affy, PROBEFILTERDB.v3$Affy)] > .70
GSEN_accepts0.75 = PROBEFILTERDB.v3$Sensitivity[match(bootModel$Affy, PROBEFILTERDB.v3$Affy)] > .75
GSEN_accepts0.80 = PROBEFILTERDB.v3$Sensitivity[match(bootModel$Affy, PROBEFILTERDB.v3$Affy)] > .80
GSEN_accepts0.85 = PROBEFILTERDB.v3$Sensitivity[match(bootModel$Affy, PROBEFILTERDB.v3$Affy)] > .85
GSEN_accepts0.90 = PROBEFILTERDB.v3$Sensitivity[match(bootModel$Affy, PROBEFILTERDB.v3$Affy)] > .90
GSEN_accepts0.95 = PROBEFILTERDB.v3$Sensitivity[match(bootModel$Affy, PROBEFILTERDB.v3$Affy)] > .95
GSEN_accepts1 = PROBEFILTERDB.v3$Sensitivity[match(bootModel$Affy, PROBEFILTERDB.v3$Affy)] == 1


#bootMergedWithPairs <- merge(bootModel,pairs, by="Affy")

bootMergedWithPairs$GSEN_accepts0 = bootMergedWithPairs$Affy %in% (PROBEFILTERDB.v3$Affy[PROBEFILTERDB.v3$Sensitivity > 0])
bootMergedWithPairs$GSEN_accepts0.05 = bootMergedWithPairs$Affy %in% (PROBEFILTERDB.v3$Affy[PROBEFILTERDB.v3$Sensitivity > 0.05])
bootMergedWithPairs$GSEN_accepts0.10 = bootMergedWithPairs$Affy %in% (PROBEFILTERDB.v3$Affy[PROBEFILTERDB.v3$Sensitivity > 0.10])
bootMergedWithPairs$GSEN_accepts0.15 = bootMergedWithPairs$Affy %in% (PROBEFILTERDB.v3$Affy[PROBEFILTERDB.v3$Sensitivity > 0.15])
bootMergedWithPairs$GSEN_accepts0.20 = bootMergedWithPairs$Affy %in% (PROBEFILTERDB.v3$Affy[PROBEFILTERDB.v3$Sensitivity > 0.20])
bootMergedWithPairs$GSEN_accepts0.25 = bootMergedWithPairs$Affy %in% (PROBEFILTERDB.v3$Affy[PROBEFILTERDB.v3$Sensitivity > 0.25])
bootMergedWithPairs$GSEN_accepts0.30 = bootMergedWithPairs$Affy %in% (PROBEFILTERDB.v3$Affy[PROBEFILTERDB.v3$Sensitivity > 0.30])
bootMergedWithPairs$GSEN_accepts0.35 = bootMergedWithPairs$Affy %in% (PROBEFILTERDB.v3$Affy[PROBEFILTERDB.v3$Sensitivity > 0.35])
bootMergedWithPairs$GSEN_accepts0.40 = bootMergedWithPairs$Affy %in% (PROBEFILTERDB.v3$Affy[PROBEFILTERDB.v3$Sensitivity > 0.40])
bootMergedWithPairs$GSEN_accepts0.45 = bootMergedWithPairs$Affy %in% (PROBEFILTERDB.v3$Affy[PROBEFILTERDB.v3$Sensitivity > 0.45])
bootMergedWithPairs$GSEN_accepts0.50 = bootMergedWithPairs$Affy %in% (PROBEFILTERDB.v3$Affy[PROBEFILTERDB.v3$Sensitivity > 0.50])
bootMergedWithPairs$GSEN_accepts0.55 = bootMergedWithPairs$Affy %in% (PROBEFILTERDB.v3$Affy[PROBEFILTERDB.v3$Sensitivity > 0.55])
bootMergedWithPairs$GSEN_accepts0.60 = bootMergedWithPairs$Affy %in% (PROBEFILTERDB.v3$Affy[PROBEFILTERDB.v3$Sensitivity > 0.60])
bootMergedWithPairs$GSEN_accepts0.65 = bootMergedWithPairs$Affy %in% (PROBEFILTERDB.v3$Affy[PROBEFILTERDB.v3$Sensitivity > 0.65])
bootMergedWithPairs$GSEN_accepts0.70 = bootMergedWithPairs$Affy %in% (PROBEFILTERDB.v3$Affy[PROBEFILTERDB.v3$Sensitivity > 0.70])
bootMergedWithPairs$GSEN_accepts0.75 = bootMergedWithPairs$Affy %in% (PROBEFILTERDB.v3$Affy[PROBEFILTERDB.v3$Sensitivity > 0.75])
bootMergedWithPairs$GSEN_accepts0.80 = bootMergedWithPairs$Affy %in% (PROBEFILTERDB.v3$Affy[PROBEFILTERDB.v3$Sensitivity > 0.80])
bootMergedWithPairs$GSEN_accepts0.85 = bootMergedWithPairs$Affy %in% (PROBEFILTERDB.v3$Affy[PROBEFILTERDB.v3$Sensitivity > 0.85])
bootMergedWithPairs$GSEN_accepts0.90 = bootMergedWithPairs$Affy %in% (PROBEFILTERDB.v3$Affy[PROBEFILTERDB.v3$Sensitivity > 0.90])
bootMergedWithPairs$GSEN_accepts0.95 = bootMergedWithPairs$Affy %in% (PROBEFILTERDB.v3$Affy[PROBEFILTERDB.v3$Sensitivity > 0.95])
bootMergedWithPairs$GSEN_accepts1 = bootMergedWithPairs$Affy %in% (PROBEFILTERDB.v3$Affy[PROBEFILTERDB.v3$Sensitivity == 1])







Lfp=1
Utp=2
deltaPlus=1
guarantee = 1e-5

resultGSEN = rbind(
  expectedUtility(label="0", dataset=bootMergedWithPairs %where% GSEN_accepts0),   # PDBA accepts
  expectedUtility(label=".05", dataset=bootMergedWithPairs %where% GSEN_accepts0.05),
  expectedUtility(label=".10", dataset=bootMergedWithPairs %where% GSEN_accepts0.10),
  expectedUtility(label=".15", dataset=bootMergedWithPairs %where% GSEN_accepts0.15),
  expectedUtility(label=".20", dataset=bootMergedWithPairs %where% GSEN_accepts0.20),
  expectedUtility(label=".25", dataset=bootMergedWithPairs %where% GSEN_accepts0.25),
  expectedUtility(label=".30", dataset=bootMergedWithPairs %where% GSEN_accepts0.30),
  expectedUtility(label=".35", dataset=bootMergedWithPairs %where% GSEN_accepts0.35),
  expectedUtility(label=".40", dataset=bootMergedWithPairs %where% GSEN_accepts0.40),
  expectedUtility(label=".45", dataset=bootMergedWithPairs %where% GSEN_accepts0.45),
  expectedUtility(label=".50", dataset=bootMergedWithPairs %where% GSEN_accepts0.50),
  expectedUtility(label=".55", dataset=bootMergedWithPairs %where% GSEN_accepts0.55),
  expectedUtility(label=".60", dataset=bootMergedWithPairs %where% GSEN_accepts0.60),
  expectedUtility(label=".65", dataset=bootMergedWithPairs %where% GSEN_accepts0.65),
  expectedUtility(label=".70", dataset=bootMergedWithPairs %where% GSEN_accepts0.70),
  expectedUtility(label=".75", dataset=bootMergedWithPairs %where% GSEN_accepts0.75),
  expectedUtility(label=".80", dataset=bootMergedWithPairs %where% GSEN_accepts0.80),
  expectedUtility(label=".85", dataset=bootMergedWithPairs %where% GSEN_accepts0.85),
  expectedUtility(label=".90", dataset=bootMergedWithPairs %where% GSEN_accepts0.90),
  expectedUtility(label=".95", dataset=bootMergedWithPairs %where% GSEN_accepts0.95),
  expectedUtility(label="1", dataset=bootMergedWithPairs %where% GSEN_accepts1)
  )


rownames(resultGSEN) = NULL
options(digits=3)
utilitytableGSEN<- resultGSEN[order(resultGSEN$Eutility),]
utilitytableGSEN


plot(resultGSEN$label,resultGSEN$Eutility,type="p",cex=.1, las=1,main="Optimize GeneAnnot Sensitivity", xlab="Sensitivity Threshold", ylab="Total Expected Utility")

plot(resultGSEN$PrPlus,resultGSEN$nPairs)

plot(resultGSEN$PrPlus, resultGSEN$nPairs, pch = "", main = "Greedy search across IDF resources by mean expected utility", xlab = "Estimated Proportion of Concordance", ylab = "Number of Pairs")
#lines(resultGSEN$PrPlus, resultGSEN$nPairs, col ="blue")
points(resultLevel1$PrPlus, resultLevel1$nPairs, col="red",pch = "1")
points(resultLevel2$PrPlus, resultLevel2$nPairs, col="orange",pch = "2")
points(resultLevel3$PrPlus, resultLevel3$nPairs, col="green", pch = "3")
points(resultLevel4$PrPlus, resultLevel4$nPairs, col="blue", pch = "4")
points(resultLevel5$PrPlus, resultLevel5$nPairs, col="purple", pch = "5")
points(resultLevel6$PrPlus, resultLevel6$nPairs, col="pink", pch = "6")
text(resultLevel1$PrPlus, resultLevel1$nPairs, labels = resultLevel1$label, pos = 1)




GSPE_accepts0 = PROBEFILTERDB.v3$Specificity[match(bootModel$Affy, PROBEFILTERDB.v3$Affy)] > 0
GSPE_accepts0.05 = PROBEFILTERDB.v3$Specificity[match(bootModel$Affy, PROBEFILTERDB.v3$Affy)] > .05
GSPE_accepts0.10 = PROBEFILTERDB.v3$Specificity[match(bootModel$Affy, PROBEFILTERDB.v3$Affy)] > .1
GSPE_accepts0.15 = PROBEFILTERDB.v3$Specificity[match(bootModel$Affy, PROBEFILTERDB.v3$Affy)] > .15
GSPE_accepts0.20 = PROBEFILTERDB.v3$Specificity[match(bootModel$Affy, PROBEFILTERDB.v3$Affy)] > .2
GSPE_accepts0.25 = PROBEFILTERDB.v3$Specificity[match(bootModel$Affy, PROBEFILTERDB.v3$Affy)] > .25
GSPE_accepts0.30 = PROBEFILTERDB.v3$Specificity[match(bootModel$Affy, PROBEFILTERDB.v3$Affy)] > .30
GSPE_accepts0.35 = PROBEFILTERDB.v3$Specificity[match(bootModel$Affy, PROBEFILTERDB.v3$Affy)] > .35
GSPE_accepts0.40 = PROBEFILTERDB.v3$Specificity[match(bootModel$Affy, PROBEFILTERDB.v3$Affy)] > .40
GSPE_accepts0.45 = PROBEFILTERDB.v3$Specificity[match(bootModel$Affy, PROBEFILTERDB.v3$Affy)] > .45
GSPE_accepts0.50 = PROBEFILTERDB.v3$Specificity[match(bootModel$Affy, PROBEFILTERDB.v3$Affy)] > .5
GSPE_accepts0.55 = PROBEFILTERDB.v3$Specificity[match(bootModel$Affy, PROBEFILTERDB.v3$Affy)] > .55
GSPE_accepts0.60 = PROBEFILTERDB.v3$Specificity[match(bootModel$Affy, PROBEFILTERDB.v3$Affy)] > .60
GSPE_accepts0.65 = PROBEFILTERDB.v3$Specificity[match(bootModel$Affy, PROBEFILTERDB.v3$Affy)] > .65
GSPE_accepts0.70 = PROBEFILTERDB.v3$Specificity[match(bootModel$Affy, PROBEFILTERDB.v3$Affy)] > .70
GSPE_accepts0.75 = PROBEFILTERDB.v3$Specificity[match(bootModel$Affy, PROBEFILTERDB.v3$Affy)] > .75
GSPE_accepts0.80 = PROBEFILTERDB.v3$Specificity[match(bootModel$Affy, PROBEFILTERDB.v3$Affy)] > .80
GSPE_accepts0.85 = PROBEFILTERDB.v3$Specificity[match(bootModel$Affy, PROBEFILTERDB.v3$Affy)] > .85
GSPE_accepts0.90 = PROBEFILTERDB.v3$Specificity[match(bootModel$Affy, PROBEFILTERDB.v3$Affy)] > .90
GSPE_accepts0.95 = PROBEFILTERDB.v3$Specificity[match(bootModel$Affy, PROBEFILTERDB.v3$Affy)] > .95
GSPE_accepts1 = PROBEFILTERDB.v3$Specificity[match(bootModel$Affy, PROBEFILTERDB.v3$Affy)] == 1


#bootMergedWithPairs <- merge(bootModel,pairs, by="Affy")

bootMergedWithPairs$GSPE_accepts0 = bootMergedWithPairs$Affy %in% (PROBEFILTERDB.v3$Affy[PROBEFILTERDB.v3$Specificity > 0])
bootMergedWithPairs$GSPE_accepts0.05 = bootMergedWithPairs$Affy %in% (PROBEFILTERDB.v3$Affy[PROBEFILTERDB.v3$Specificity > 0.05])
bootMergedWithPairs$GSPE_accepts0.10 = bootMergedWithPairs$Affy %in% (PROBEFILTERDB.v3$Affy[PROBEFILTERDB.v3$Specificity > 0.10])
bootMergedWithPairs$GSPE_accepts0.15 = bootMergedWithPairs$Affy %in% (PROBEFILTERDB.v3$Affy[PROBEFILTERDB.v3$Specificity > 0.15])
bootMergedWithPairs$GSPE_accepts0.20 = bootMergedWithPairs$Affy %in% (PROBEFILTERDB.v3$Affy[PROBEFILTERDB.v3$Specificity > 0.20])
bootMergedWithPairs$GSPE_accepts0.25 = bootMergedWithPairs$Affy %in% (PROBEFILTERDB.v3$Affy[PROBEFILTERDB.v3$Specificity > 0.25])
bootMergedWithPairs$GSPE_accepts0.30 = bootMergedWithPairs$Affy %in% (PROBEFILTERDB.v3$Affy[PROBEFILTERDB.v3$Specificity > 0.30])
bootMergedWithPairs$GSPE_accepts0.35 = bootMergedWithPairs$Affy %in% (PROBEFILTERDB.v3$Affy[PROBEFILTERDB.v3$Specificity > 0.35])
bootMergedWithPairs$GSPE_accepts0.40 = bootMergedWithPairs$Affy %in% (PROBEFILTERDB.v3$Affy[PROBEFILTERDB.v3$Specificity > 0.40])
bootMergedWithPairs$GSPE_accepts0.45 = bootMergedWithPairs$Affy %in% (PROBEFILTERDB.v3$Affy[PROBEFILTERDB.v3$Specificity > 0.45])
bootMergedWithPairs$GSPE_accepts0.50 = bootMergedWithPairs$Affy %in% (PROBEFILTERDB.v3$Affy[PROBEFILTERDB.v3$Specificity > 0.50])
bootMergedWithPairs$GSPE_accepts0.55 = bootMergedWithPairs$Affy %in% (PROBEFILTERDB.v3$Affy[PROBEFILTERDB.v3$Specificity > 0.55])
bootMergedWithPairs$GSPE_accepts0.60 = bootMergedWithPairs$Affy %in% (PROBEFILTERDB.v3$Affy[PROBEFILTERDB.v3$Specificity > 0.60])
bootMergedWithPairs$GSPE_accepts0.65 = bootMergedWithPairs$Affy %in% (PROBEFILTERDB.v3$Affy[PROBEFILTERDB.v3$Specificity > 0.65])
bootMergedWithPairs$GSPE_accepts0.70 = bootMergedWithPairs$Affy %in% (PROBEFILTERDB.v3$Affy[PROBEFILTERDB.v3$Specificity > 0.70])
bootMergedWithPairs$GSPE_accepts0.75 = bootMergedWithPairs$Affy %in% (PROBEFILTERDB.v3$Affy[PROBEFILTERDB.v3$Specificity > 0.75])
bootMergedWithPairs$GSPE_accepts0.80 = bootMergedWithPairs$Affy %in% (PROBEFILTERDB.v3$Affy[PROBEFILTERDB.v3$Specificity > 0.80])
bootMergedWithPairs$GSPE_accepts0.85 = bootMergedWithPairs$Affy %in% (PROBEFILTERDB.v3$Affy[PROBEFILTERDB.v3$Specificity > 0.85])
bootMergedWithPairs$GSPE_accepts0.90 = bootMergedWithPairs$Affy %in% (PROBEFILTERDB.v3$Affy[PROBEFILTERDB.v3$Specificity > 0.90])
bootMergedWithPairs$GSPE_accepts0.95 = bootMergedWithPairs$Affy %in% (PROBEFILTERDB.v3$Affy[PROBEFILTERDB.v3$Specificity > 0.95])
bootMergedWithPairs$GSPE_accepts1 = bootMergedWithPairs$Affy %in% (PROBEFILTERDB.v3$Affy[PROBEFILTERDB.v3$Specificity == 1])







Lfp=1
Utp=2
deltaPlus=1/2
guarantee = 1e-5

resultGSPE = rbind(
  expectedUtility(label="0", dataset=bootMergedWithPairs %where% GSPE_accepts0),   # PDBA accepts
  expectedUtility(label=".05", dataset=bootMergedWithPairs %where% GSPE_accepts0.05),
  expectedUtility(label=".10", dataset=bootMergedWithPairs %where% GSPE_accepts0.10),
  expectedUtility(label=".15", dataset=bootMergedWithPairs %where% GSPE_accepts0.15),
  expectedUtility(label=".20", dataset=bootMergedWithPairs %where% GSPE_accepts0.20),
  expectedUtility(label=".25", dataset=bootMergedWithPairs %where% GSPE_accepts0.25),
  expectedUtility(label=".30", dataset=bootMergedWithPairs %where% GSPE_accepts0.30),
  expectedUtility(label=".35", dataset=bootMergedWithPairs %where% GSPE_accepts0.35),
  expectedUtility(label=".40", dataset=bootMergedWithPairs %where% GSPE_accepts0.40),
  expectedUtility(label=".45", dataset=bootMergedWithPairs %where% GSPE_accepts0.45),
  expectedUtility(label=".50", dataset=bootMergedWithPairs %where% GSPE_accepts0.50),
  expectedUtility(label=".55", dataset=bootMergedWithPairs %where% GSPE_accepts0.55),
  expectedUtility(label=".60", dataset=bootMergedWithPairs %where% GSPE_accepts0.60),
  expectedUtility(label=".65", dataset=bootMergedWithPairs %where% GSPE_accepts0.65),
  expectedUtility(label=".70", dataset=bootMergedWithPairs %where% GSPE_accepts0.70),
  expectedUtility(label=".75", dataset=bootMergedWithPairs %where% GSPE_accepts0.75),
  expectedUtility(label=".80", dataset=bootMergedWithPairs %where% GSPE_accepts0.80),
  expectedUtility(label=".85", dataset=bootMergedWithPairs %where% GSPE_accepts0.85),
  expectedUtility(label=".90", dataset=bootMergedWithPairs %where% GSPE_accepts0.90),
  expectedUtility(label=".95", dataset=bootMergedWithPairs %where% GSPE_accepts0.95),
  expectedUtility(label="1", dataset=bootMergedWithPairs %where% GSPE_accepts1)
  )


rownames(resultGSPE) = NULL
options(digits=3)
utilitytableGSPE<- resultGSPE[order(resultGSPE$Eutility),]
utilitytableGSPE

plot(resultGSPE$label,resultGSPE$Eutility,type="p",cex=.1, las=1,main="Optimize GeneAnnot Specificity", xlab="Probe Specificity Threshold", ylab="Total Expected Utility")


plot(resultGSPE$PrPlus,resultGSPE$nPairs)


############JETSET GREEDY SEARCH

Lfp=1
Utp=2
deltaPlus=1
guarantee = 1e-5

resultAllTEU = rbind(
  expectedUtility(label="Use All", dataset=bootMergedWithPairs),  # All (no filtering)
  ###level1
  expectedUtility(label="PD", dataset=bootMergedWithPairs %where% PDBA_accepts),   # PDBA accepts
  expectedUtility(label="M", dataset=bootMergedWithPairs %where% Masker_accepts),   # “Masker accepts”
  expectedUtility(label="AT", dataset=bootMergedWithPairs %where% Affy_tag_accepts),
  expectedUtility(label="AG", dataset=bootMergedWithPairs %where% Affy_Grade_accepts),
  expectedUtility(label="E", dataset=bootMergedWithPairs %where% ENCODE_accepts),
  expectedUtility(label="GQ", dataset=bootMergedWithPairs %where% Geneannot_quality_accepts),
  expectedUtility(label="GSEN", dataset=bootMergedWithPairs %where% Geneannot_sensitivity_accepts),
  expectedUtility(label="GSPE", dataset=bootMergedWithPairs %where% Geneannot_specificity_accepts),
  expectedUtility(label="J", dataset=bootMergedWithPairs %where% Jetset_accepts),

  ###Level2
  expectedUtility(label="+int(PD)", dataset=bootMergedWithPairs %where% (PDBA_accepts & Jetset_accepts)),
  expectedUtility(label="+int(M)", dataset=bootMergedWithPairs %where% (Masker_accepts & Jetset_accepts)),
  expectedUtility(label="+int(AT)", dataset=bootMergedWithPairs %where% (Affy_tag_accepts & Jetset_accepts)),
  expectedUtility(label="+int(AG)", dataset=bootMergedWithPairs %where% (Affy_Grade_accepts & Jetset_accepts)),
  expectedUtility(label="+int(E)", dataset=bootMergedWithPairs %where% (ENCODE_accepts & Jetset_accepts)),
  expectedUtility(label="+int(GQ)", dataset=bootMergedWithPairs %where% (Geneannot_quality_accepts & Jetset_accepts)),
  expectedUtility(label="+int(GSEN)", dataset=bootMergedWithPairs %where% (Geneannot_sensitivity_accepts & Jetset_accepts)),
  expectedUtility(label="+int(GSPE)", dataset=bootMergedWithPairs %where% (Geneannot_specificity_accepts & Jetset_accepts)),
  expectedUtility(label="+union(PD)", dataset=bootMergedWithPairs %where% (PDBA_accepts | Jetset_accepts)),
  expectedUtility(label="+union(M)", dataset=bootMergedWithPairs %where% (Masker_accepts | Jetset_accepts)),
  expectedUtility(label="+union(AT)", dataset=bootMergedWithPairs %where% (Affy_tag_accepts | Jetset_accepts)),
  expectedUtility(label="+union(AG)", dataset=bootMergedWithPairs %where% (Affy_Grade_accepts | Jetset_accepts)),
  expectedUtility(label="+union(E)", dataset=bootMergedWithPairs %where% (ENCODE_accepts | Jetset_accepts)),
  expectedUtility(label="+union(GQ)", dataset=bootMergedWithPairs %where% (Geneannot_quality_accepts | Jetset_accepts)),
  expectedUtility(label="+union(GSEN)", dataset=bootMergedWithPairs %where% (Geneannot_sensitivity_accepts | Jetset_accepts)),
  expectedUtility(label="+union(GSPE)", dataset=bootMergedWithPairs %where% (Geneannot_specificity_accepts | Jetset_accepts)),

  #######Level3
  expectedUtility(label="+int(PD)", dataset=bootMergedWithPairs %where% (PDBA_accepts & Jetset_accepts & Geneannot_specificity_accepts)),
  expectedUtility(label="+int(M)", dataset=bootMergedWithPairs %where% (Masker_accepts & Jetset_accepts & Geneannot_specificity_accepts)),
  expectedUtility(label="+int(AT)", dataset=bootMergedWithPairs %where% (Affy_tag_accepts & Jetset_accepts & Geneannot_specificity_accepts)),
  expectedUtility(label="+int(AG)", dataset=bootMergedWithPairs %where% (Affy_Grade_accepts & Jetset_accepts & Geneannot_specificity_accepts)),
  expectedUtility(label="+int(E)", dataset=bootMergedWithPairs %where% (ENCODE_accepts & Jetset_accepts & Geneannot_specificity_accepts)),
  expectedUtility(label="+int(GQ)", dataset=bootMergedWithPairs %where% (Geneannot_quality_accepts & Jetset_accepts & Geneannot_specificity_accepts)),
  expectedUtility(label="+int(GSEN)", dataset=bootMergedWithPairs %where% (Geneannot_sensitivity_accepts & Jetset_accepts & Geneannot_specificity_accepts)),
  expectedUtility(label="+int(GSPE)", dataset=bootMergedWithPairs %where% (Geneannot_specificity_accepts & Jetset_accepts & Geneannot_specificity_accepts)),
  expectedUtility(label="+union(PD)", dataset=bootMergedWithPairs %where% (PDBA_accepts | Jetset_accepts & Geneannot_specificity_accepts)),
  expectedUtility(label="+union(M)", dataset=bootMergedWithPairs %where% (Masker_accepts | Jetset_accepts & Geneannot_specificity_accepts)),
  expectedUtility(label="+union(AT)", dataset=bootMergedWithPairs %where% (Affy_tag_accepts | Jetset_accepts & Geneannot_specificity_accepts)),
  expectedUtility(label="+union(AG)", dataset=bootMergedWithPairs %where% (Affy_Grade_accepts | Jetset_accepts & Geneannot_specificity_accepts)),
  expectedUtility(label="+union(E)", dataset=bootMergedWithPairs %where% (ENCODE_accepts | Jetset_accepts & Geneannot_specificity_accepts)),
  expectedUtility(label="+union(GQ)", dataset=bootMergedWithPairs %where% (Geneannot_quality_accepts | Jetset_accepts & Geneannot_specificity_accepts)),
  expectedUtility(label="+union(GSEN)", dataset=bootMergedWithPairs %where% (Geneannot_sensitivity_accepts | Jetset_accepts & Geneannot_specificity_accepts)),
  expectedUtility(label="+union(GSPE)", dataset=bootMergedWithPairs %where% (Geneannot_specificity_accepts | Jetset_accepts & Geneannot_specificity_accepts)),


  #######Level4

  expectedUtility(label="+int(PD)", dataset=bootMergedWithPairs %where% (PDBA_accepts & Jetset_accepts & Geneannot_specificity_accepts & Geneannot_sensitivity_accepts)),
  expectedUtility(label="+int(M)", dataset=bootMergedWithPairs %where% (Masker_accepts & Jetset_accepts & Geneannot_specificity_accepts & Geneannot_sensitivity_accepts)),
  expectedUtility(label="+int(AT)", dataset=bootMergedWithPairs %where% (Affy_tag_accepts & Jetset_accepts & Geneannot_specificity_accepts & Geneannot_sensitivity_accepts)),
  expectedUtility(label="+int(AG)", dataset=bootMergedWithPairs %where% (Affy_Grade_accepts & Jetset_accepts & Geneannot_specificity_accepts & Geneannot_sensitivity_accepts)),
  expectedUtility(label="+int(E)", dataset=bootMergedWithPairs %where% (ENCODE_accepts & Jetset_accepts & Geneannot_specificity_accepts & Geneannot_sensitivity_accepts)),
  expectedUtility(label="+int(GQ)", dataset=bootMergedWithPairs %where% (Geneannot_quality_accepts & Jetset_accepts & Geneannot_specificity_accepts & Geneannot_sensitivity_accepts)),
  expectedUtility(label="+int(GSEN)", dataset=bootMergedWithPairs %where% (Geneannot_sensitivity_accepts & Jetset_accepts & Geneannot_specificity_accepts & Geneannot_sensitivity_accepts)),
  expectedUtility(label="+int(GSPE)", dataset=bootMergedWithPairs %where% (Geneannot_specificity_accepts & Jetset_accepts & Geneannot_specificity_accepts & Geneannot_sensitivity_accepts)),
  expectedUtility(label="+union(PD)", dataset=bootMergedWithPairs %where% (PDBA_accepts | Jetset_accepts & Geneannot_specificity_accepts & Geneannot_sensitivity_accepts)),
  expectedUtility(label="+union(M)", dataset=bootMergedWithPairs %where% (Masker_accepts | Jetset_accepts & Geneannot_specificity_accepts & Geneannot_sensitivity_accepts)),
  expectedUtility(label="+union(AT)", dataset=bootMergedWithPairs %where% (Affy_tag_accepts | Jetset_accepts & Geneannot_specificity_accepts & Geneannot_sensitivity_accepts)),
  expectedUtility(label="+union(AG)", dataset=bootMergedWithPairs %where% (Affy_Grade_accepts | Jetset_accepts & Geneannot_specificity_accepts & Geneannot_sensitivity_accepts)),
  expectedUtility(label="+union(E)", dataset=bootMergedWithPairs %where% (ENCODE_accepts | Jetset_accepts & Geneannot_specificity_accepts & Geneannot_sensitivity_accepts)),
  expectedUtility(label="+union(GQ)", dataset=bootMergedWithPairs %where% (Geneannot_quality_accepts | Jetset_accepts & Geneannot_specificity_accepts & Geneannot_sensitivity_accepts)),
  expectedUtility(label="+union(GSEN)", dataset=bootMergedWithPairs %where% (Geneannot_sensitivity_accepts | Jetset_accepts & Geneannot_specificity_accepts & Geneannot_sensitivity_accepts)),
  expectedUtility(label="+union(GSPE)", dataset=bootMergedWithPairs %where% (Geneannot_specificity_accepts | Jetset_accepts & Geneannot_specificity_accepts & Geneannot_sensitivity_accepts)),

  ###############Level5

  expectedUtility(label="+int(PD)", dataset=bootMergedWithPairs %where% (PDBA_accepts & Jetset_accepts & Geneannot_specificity_accepts & Geneannot_sensitivity_accepts & Affy_Grade_accepts)),
  expectedUtility(label="+int(M)", dataset=bootMergedWithPairs %where% (Masker_accepts & Jetset_accepts & Geneannot_specificity_accepts & Geneannot_sensitivity_accepts & Affy_Grade_accepts)),
  expectedUtility(label="+int(AT)", dataset=bootMergedWithPairs %where% (Affy_tag_accepts & Jetset_accepts & Geneannot_specificity_accepts & Geneannot_sensitivity_accepts & Affy_Grade_accepts)),
  expectedUtility(label="+int(AG)", dataset=bootMergedWithPairs %where% (Affy_Grade_accepts & Jetset_accepts & Geneannot_specificity_accepts & Geneannot_sensitivity_accepts & Affy_Grade_accepts)),
  expectedUtility(label="+int(E)", dataset=bootMergedWithPairs %where% (ENCODE_accepts & Jetset_accepts & Geneannot_specificity_accepts & Geneannot_sensitivity_accepts & Affy_Grade_accepts)),
  expectedUtility(label="+int(GQ)", dataset=bootMergedWithPairs %where% (Geneannot_quality_accepts & Jetset_accepts & Geneannot_specificity_accepts & Geneannot_sensitivity_accepts & Affy_Grade_accepts)),
  expectedUtility(label="+int(GSEN)", dataset=bootMergedWithPairs %where% (Geneannot_sensitivity_accepts & Jetset_accepts & Geneannot_specificity_accepts & Geneannot_sensitivity_accepts & Affy_Grade_accepts)),
  expectedUtility(label="+int(GSPE)", dataset=bootMergedWithPairs %where% (Geneannot_specificity_accepts & Jetset_accepts & Geneannot_specificity_accepts & Geneannot_sensitivity_accepts & Affy_Grade_accepts)),
  expectedUtility(label="+union(PD)", dataset=bootMergedWithPairs %where% (PDBA_accepts | Jetset_accepts & Geneannot_specificity_accepts & Geneannot_sensitivity_accepts & Affy_Grade_accepts)),
  expectedUtility(label="+union(M)", dataset=bootMergedWithPairs %where% (Masker_accepts | Jetset_accepts & Geneannot_specificity_accepts & Geneannot_sensitivity_accepts & Affy_Grade_accepts)),
  expectedUtility(label="+union(AT)", dataset=bootMergedWithPairs %where% (Affy_tag_accepts | Jetset_accepts & Geneannot_specificity_accepts & Geneannot_sensitivity_accepts & Affy_Grade_accepts)),
  expectedUtility(label="+union(AG)", dataset=bootMergedWithPairs %where% (Affy_Grade_accepts | Jetset_accepts & Geneannot_specificity_accepts & Geneannot_sensitivity_accepts & Affy_Grade_accepts)),
  expectedUtility(label="+union(E)", dataset=bootMergedWithPairs %where% (ENCODE_accepts | Jetset_accepts & Geneannot_specificity_accepts & Geneannot_sensitivity_accepts & Affy_Grade_accepts)),
  expectedUtility(label="+union(GQ)", dataset=bootMergedWithPairs %where% (Geneannot_quality_accepts | Jetset_accepts & Geneannot_specificity_accepts & Geneannot_sensitivity_accepts & Affy_Grade_accepts)),
  expectedUtility(label="+union(GSEN)", dataset=bootMergedWithPairs %where% (Geneannot_sensitivity_accepts | Jetset_accepts & Geneannot_specificity_accepts & Geneannot_sensitivity_accepts & Affy_Grade_accepts)),
  expectedUtility(label="+union(GSPE)", dataset=bootMergedWithPairs %where% (Geneannot_specificity_accepts | Jetset_accepts & Geneannot_specificity_accepts & Geneannot_sensitivity_accepts & Affy_Grade_accepts)),


  ################Level6

  expectedUtility(label="+int(PD)", dataset=bootMergedWithPairs %where% (PDBA_accepts & Jetset_accepts & Geneannot_specificity_accepts & Geneannot_sensitivity_accepts & Affy_Grade_accepts & Geneannot_quality_accepts)),
  expectedUtility(label="+int(M)", dataset=bootMergedWithPairs %where% (Masker_accepts & Jetset_accepts & Geneannot_specificity_accepts & Geneannot_sensitivity_accepts & Affy_Grade_accepts & Geneannot_quality_accepts)),
  expectedUtility(label="+int(AT)", dataset=bootMergedWithPairs %where% (Affy_tag_accepts & Jetset_accepts & Geneannot_specificity_accepts & Geneannot_sensitivity_accepts & Affy_Grade_accepts & Geneannot_quality_accepts)),
  expectedUtility(label="+int(AG)", dataset=bootMergedWithPairs %where% (Affy_Grade_accepts & Jetset_accepts & Geneannot_specificity_accepts & Geneannot_sensitivity_accepts & Affy_Grade_accepts & Geneannot_quality_accepts)),
  expectedUtility(label="+int(E)", dataset=bootMergedWithPairs %where% (ENCODE_accepts & Jetset_accepts & Geneannot_specificity_accepts & Geneannot_sensitivity_accepts & Affy_Grade_accepts & Geneannot_quality_accepts)),
  expectedUtility(label="+int(GQ)", dataset=bootMergedWithPairs %where% (Geneannot_quality_accepts & Jetset_accepts & Geneannot_specificity_accepts & Geneannot_sensitivity_accepts & Affy_Grade_accepts & Geneannot_quality_accepts)),
  expectedUtility(label="+int(GSEN)", dataset=bootMergedWithPairs %where% (Geneannot_sensitivity_accepts & Jetset_accepts & Geneannot_specificity_accepts & Geneannot_sensitivity_accepts & Affy_Grade_accepts & Geneannot_quality_accepts)),
  expectedUtility(label="+int(GSPE)", dataset=bootMergedWithPairs %where% (Geneannot_specificity_accepts & Jetset_accepts & Geneannot_specificity_accepts & Geneannot_sensitivity_accepts & Affy_Grade_accepts & Geneannot_quality_accepts)),
  expectedUtility(label="+union(PD)", dataset=bootMergedWithPairs %where% (PDBA_accepts | Jetset_accepts & Geneannot_specificity_accepts & Geneannot_sensitivity_accepts & Affy_Grade_accepts & Geneannot_quality_accepts)),
  expectedUtility(label="+union(M)", dataset=bootMergedWithPairs %where% (Masker_accepts | Jetset_accepts & Geneannot_specificity_accepts & Geneannot_sensitivity_accepts & Affy_Grade_accepts & Geneannot_quality_accepts)),
  expectedUtility(label="+union(AT)", dataset=bootMergedWithPairs %where% (Affy_tag_accepts | Jetset_accepts & Geneannot_specificity_accepts & Geneannot_sensitivity_accepts & Affy_Grade_accepts & Geneannot_quality_accepts)),
  expectedUtility(label="+union(AG)", dataset=bootMergedWithPairs %where% (Affy_Grade_accepts | Jetset_accepts & Geneannot_specificity_accepts & Geneannot_sensitivity_accepts & Affy_Grade_accepts & Geneannot_quality_accepts)),
  expectedUtility(label="+union(E)", dataset=bootMergedWithPairs %where% (ENCODE_accepts | Jetset_accepts & Geneannot_specificity_accepts & Geneannot_sensitivity_accepts & Affy_Grade_accepts & Geneannot_quality_accepts)),
  expectedUtility(label="+union(GQ)", dataset=bootMergedWithPairs %where% (Geneannot_quality_accepts | Jetset_accepts & Geneannot_specificity_accepts & Geneannot_sensitivity_accepts & Affy_Grade_accepts & Geneannot_quality_accepts)),
  expectedUtility(label="+union(GSEN)", dataset=bootMergedWithPairs %where% (Geneannot_sensitivity_accepts | Jetset_accepts & Geneannot_specificity_accepts & Geneannot_sensitivity_accepts & Affy_Grade_accepts & Geneannot_quality_accepts)),
  expectedUtility(label="+union(GSPE)", dataset=bootMergedWithPairs %where% (Geneannot_specificity_accepts | Jetset_accepts & Geneannot_specificity_accepts & Geneannot_sensitivity_accepts & Affy_Grade_accepts & Geneannot_quality_accepts))


  )


########################################################


resultTEUOptimize = rbind(
  expectedUtility(label="Use All", dataset=bootMergedWithPairs),  # All (no filtering)
  expectedUtility(label="J", dataset=bootMergedWithPairs %where% Jetset_accepts),
  expectedUtility(label="+int(GSPE)", dataset=bootMergedWithPairs %where% (Jetset_accepts & Geneannot_specificity_accepts)),
  expectedUtility(label="+int(GSEN)", dataset=bootMergedWithPairs %where% (Jetset_accepts & Geneannot_specificity_accepts & Geneannot_sensitivity_accepts)),
  expectedUtility(label="+int(AG)", dataset=bootMergedWithPairs %where% (Jetset_accepts & Geneannot_sensitivity_accepts & Geneannot_specificity_accepts & Affy_Grade_accepts)),
  expectedUtility(label="+int(GQ)", dataset=bootMergedWithPairs %where% (Jetset_accepts & Geneannot_sensitivity_accepts & Geneannot_specificity_accepts & Affy_Grade_accepts & Geneannot_quality_accepts))
  )

plot(resultTEUOptimize$label,resultTEUOptimize$Eutility, type="p",cex=.05, las=1,main="Greedy Search Summary", xlab="Probeset Resource Set", ylab="Total Expected Utility")



##########################################################

resultLevel1TEU = rbind(
  expectedUtility(label="PD", dataset=bootMergedWithPairs %where% PDBA_accepts),   # PDBA accepts
  expectedUtility(label="M", dataset=bootMergedWithPairs %where% Masker_accepts),   # “Masker accepts”
  expectedUtility(label="AT", dataset=bootMergedWithPairs %where% Affy_tag_accepts),
  expectedUtility(label="AG", dataset=bootMergedWithPairs %where% Affy_Grade_accepts),
  expectedUtility(label="E", dataset=bootMergedWithPairs %where% ENCODE_accepts),
  expectedUtility(label="GQ", dataset=bootMergedWithPairs %where% Geneannot_quality_accepts),
  expectedUtility(label="GSEN", dataset=bootMergedWithPairs %where% Geneannot_sensitivity_accepts),
  expectedUtility(label="GSPE", dataset=bootMergedWithPairs %where% Geneannot_specificity_accepts),
  expectedUtility(label="J", dataset=bootMergedWithPairs %where% Jetset_accepts)
  )

###########################################################

resultLevel2TEU = rbind(
  expectedUtility(label="+int(PD)", dataset=bootMergedWithPairs %where% (PDBA_accepts & Jetset_accepts)),
  expectedUtility(label="+int(M)", dataset=bootMergedWithPairs %where% (Masker_accepts & Jetset_accepts)),
  expectedUtility(label="+int(AT)", dataset=bootMergedWithPairs %where% (Affy_tag_accepts & Jetset_accepts)),
  expectedUtility(label="+int(AG)", dataset=bootMergedWithPairs %where% (Affy_Grade_accepts & Jetset_accepts)),
  expectedUtility(label="+int(E)", dataset=bootMergedWithPairs %where% (ENCODE_accepts & Jetset_accepts)),
  expectedUtility(label="+int(GQ)", dataset=bootMergedWithPairs %where% (Geneannot_quality_accepts & Jetset_accepts)),
  expectedUtility(label="+int(GSEN)", dataset=bootMergedWithPairs %where% (Geneannot_sensitivity_accepts & Jetset_accepts)),
  expectedUtility(label="+int(GSPE)", dataset=bootMergedWithPairs %where% (Geneannot_specificity_accepts & Jetset_accepts)),
  expectedUtility(label="+union(PD)", dataset=bootMergedWithPairs %where% (PDBA_accepts | Jetset_accepts)),
  expectedUtility(label="+union(M)", dataset=bootMergedWithPairs %where% (Masker_accepts | Jetset_accepts)),
  expectedUtility(label="+union(AT)", dataset=bootMergedWithPairs %where% (Affy_tag_accepts | Jetset_accepts)),
  expectedUtility(label="+union(AG)", dataset=bootMergedWithPairs %where% (Affy_Grade_accepts | Jetset_accepts)),
  expectedUtility(label="+union(E)", dataset=bootMergedWithPairs %where% (ENCODE_accepts | Jetset_accepts)),
  expectedUtility(label="+union(GQ)", dataset=bootMergedWithPairs %where% (Geneannot_quality_accepts | Jetset_accepts)),
  expectedUtility(label="+union(GSEN)", dataset=bootMergedWithPairs %where% (Geneannot_sensitivity_accepts | Jetset_accepts)),
  expectedUtility(label="+union(GSPE)", dataset=bootMergedWithPairs %where% (Geneannot_specificity_accepts | Jetset_accepts))
  )


#################################
resultLevel3TEU =rbind(
  expectedUtility(label="+int(PD)", dataset=bootMergedWithPairs %where% (PDBA_accepts & Jetset_accepts & Geneannot_specificity_accepts)),
  expectedUtility(label="+int(M)", dataset=bootMergedWithPairs %where% (Masker_accepts & Jetset_accepts & Geneannot_specificity_accepts)),
  expectedUtility(label="+int(AT)", dataset=bootMergedWithPairs %where% (Affy_tag_accepts & Jetset_accepts & Geneannot_specificity_accepts)),
  expectedUtility(label="+int(AG)", dataset=bootMergedWithPairs %where% (Affy_Grade_accepts & Jetset_accepts & Geneannot_specificity_accepts)),
  expectedUtility(label="+int(E)", dataset=bootMergedWithPairs %where% (ENCODE_accepts & Jetset_accepts & Geneannot_specificity_accepts)),
  expectedUtility(label="+int(GQ)", dataset=bootMergedWithPairs %where% (Geneannot_quality_accepts & Jetset_accepts & Geneannot_specificity_accepts)),
  expectedUtility(label="+int(GSEN)", dataset=bootMergedWithPairs %where% (Geneannot_sensitivity_accepts & Jetset_accepts & Geneannot_specificity_accepts)),
  expectedUtility(label="+int(GSPE)", dataset=bootMergedWithPairs %where% (Geneannot_specificity_accepts & Jetset_accepts & Geneannot_specificity_accepts)),
  expectedUtility(label="+union(PD)", dataset=bootMergedWithPairs %where% (PDBA_accepts | Jetset_accepts & Geneannot_specificity_accepts)),
  expectedUtility(label="+union(M)", dataset=bootMergedWithPairs %where% (Masker_accepts | Jetset_accepts & Geneannot_specificity_accepts)),
  expectedUtility(label="+union(AT)", dataset=bootMergedWithPairs %where% (Affy_tag_accepts | Jetset_accepts & Geneannot_specificity_accepts)),
  expectedUtility(label="+union(AG)", dataset=bootMergedWithPairs %where% (Affy_Grade_accepts | Jetset_accepts & Geneannot_specificity_accepts)),
  expectedUtility(label="+union(E)", dataset=bootMergedWithPairs %where% (ENCODE_accepts | Jetset_accepts & Geneannot_specificity_accepts)),
  expectedUtility(label="+union(GQ)", dataset=bootMergedWithPairs %where% (Geneannot_quality_accepts | Jetset_accepts & Geneannot_specificity_accepts)),
  expectedUtility(label="+union(GSEN)", dataset=bootMergedWithPairs %where% (Geneannot_sensitivity_accepts | Jetset_accepts & Geneannot_specificity_accepts)),
  expectedUtility(label="+union(GSPE)", dataset=bootMergedWithPairs %where% (Geneannot_specificity_accepts | Jetset_accepts & Geneannot_specificity_accepts))

  )


##############################
resultLevel4TEU = rbind(
  expectedUtility(label="+int(PD)", dataset=bootMergedWithPairs %where% (PDBA_accepts & Jetset_accepts & Geneannot_specificity_accepts & Geneannot_sensitivity_accepts)),
  expectedUtility(label="+int(M)", dataset=bootMergedWithPairs %where% (Masker_accepts & Jetset_accepts & Geneannot_specificity_accepts & Geneannot_sensitivity_accepts)),
  expectedUtility(label="+int(AT)", dataset=bootMergedWithPairs %where% (Affy_tag_accepts & Jetset_accepts & Geneannot_specificity_accepts & Geneannot_sensitivity_accepts)),
  expectedUtility(label="+int(AG)", dataset=bootMergedWithPairs %where% (Affy_Grade_accepts & Jetset_accepts & Geneannot_specificity_accepts & Geneannot_sensitivity_accepts)),
  expectedUtility(label="+int(E)", dataset=bootMergedWithPairs %where% (ENCODE_accepts & Jetset_accepts & Geneannot_specificity_accepts & Geneannot_sensitivity_accepts)),
  expectedUtility(label="+int(GQ)", dataset=bootMergedWithPairs %where% (Geneannot_quality_accepts & Jetset_accepts & Geneannot_specificity_accepts & Geneannot_sensitivity_accepts)),
  expectedUtility(label="+int(GSEN)", dataset=bootMergedWithPairs %where% (Geneannot_sensitivity_accepts & Jetset_accepts & Geneannot_specificity_accepts & Geneannot_sensitivity_accepts)),
  expectedUtility(label="+int(GSPE)", dataset=bootMergedWithPairs %where% (Geneannot_specificity_accepts & Jetset_accepts & Geneannot_specificity_accepts & Geneannot_sensitivity_accepts)),
  expectedUtility(label="+union(PD)", dataset=bootMergedWithPairs %where% (PDBA_accepts | Jetset_accepts & Geneannot_specificity_accepts & Geneannot_sensitivity_accepts)),
  expectedUtility(label="+union(M)", dataset=bootMergedWithPairs %where% (Masker_accepts | Jetset_accepts & Geneannot_specificity_accepts & Geneannot_sensitivity_accepts)),
  expectedUtility(label="+union(AT)", dataset=bootMergedWithPairs %where% (Affy_tag_accepts | Jetset_accepts & Geneannot_specificity_accepts & Geneannot_sensitivity_accepts)),
  expectedUtility(label="+union(AG)", dataset=bootMergedWithPairs %where% (Affy_Grade_accepts | Jetset_accepts & Geneannot_specificity_accepts & Geneannot_sensitivity_accepts)),
  expectedUtility(label="+union(E)", dataset=bootMergedWithPairs %where% (ENCODE_accepts | Jetset_accepts & Geneannot_specificity_accepts & Geneannot_sensitivity_accepts)),
  expectedUtility(label="+union(GQ)", dataset=bootMergedWithPairs %where% (Geneannot_quality_accepts | Jetset_accepts & Geneannot_specificity_accepts & Geneannot_sensitivity_accepts)),
  expectedUtility(label="+union(GSEN)", dataset=bootMergedWithPairs %where% (Geneannot_sensitivity_accepts | Jetset_accepts & Geneannot_specificity_accepts & Geneannot_sensitivity_accepts)),
  expectedUtility(label="+union(GSPE)", dataset=bootMergedWithPairs %where% (Geneannot_specificity_accepts | Jetset_accepts & Geneannot_specificity_accepts & Geneannot_sensitivity_accepts))
  )

###########################


resultLevel5TEU = rbind(
  expectedUtility(label="+int(PD)", dataset=bootMergedWithPairs %where% (PDBA_accepts & Jetset_accepts & Geneannot_specificity_accepts & Geneannot_sensitivity_accepts & Affy_Grade_accepts)),
  expectedUtility(label="+int(M)", dataset=bootMergedWithPairs %where% (Masker_accepts & Jetset_accepts & Geneannot_specificity_accepts & Geneannot_sensitivity_accepts & Affy_Grade_accepts)),
  expectedUtility(label="+int(AT)", dataset=bootMergedWithPairs %where% (Affy_tag_accepts & Jetset_accepts & Geneannot_specificity_accepts & Geneannot_sensitivity_accepts & Affy_Grade_accepts)),
  expectedUtility(label="+int(AG)", dataset=bootMergedWithPairs %where% (Affy_Grade_accepts & Jetset_accepts & Geneannot_specificity_accepts & Geneannot_sensitivity_accepts & Affy_Grade_accepts)),
  expectedUtility(label="+int(E)", dataset=bootMergedWithPairs %where% (ENCODE_accepts & Jetset_accepts & Geneannot_specificity_accepts & Geneannot_sensitivity_accepts & Affy_Grade_accepts)),
  expectedUtility(label="+int(GQ)", dataset=bootMergedWithPairs %where% (Geneannot_quality_accepts & Jetset_accepts & Geneannot_specificity_accepts & Geneannot_sensitivity_accepts & Affy_Grade_accepts)),
  expectedUtility(label="+int(GSEN)", dataset=bootMergedWithPairs %where% (Geneannot_sensitivity_accepts & Jetset_accepts & Geneannot_specificity_accepts & Geneannot_sensitivity_accepts & Affy_Grade_accepts)),
  expectedUtility(label="+int(GSPE)", dataset=bootMergedWithPairs %where% (Geneannot_specificity_accepts & Jetset_accepts & Geneannot_specificity_accepts & Geneannot_sensitivity_accepts & Affy_Grade_accepts)),
  expectedUtility(label="+union(PD)", dataset=bootMergedWithPairs %where% (PDBA_accepts | Jetset_accepts & Geneannot_specificity_accepts & Geneannot_sensitivity_accepts & Affy_Grade_accepts)),
  expectedUtility(label="+union(M)", dataset=bootMergedWithPairs %where% (Masker_accepts | Jetset_accepts & Geneannot_specificity_accepts & Geneannot_sensitivity_accepts & Affy_Grade_accepts)),
  expectedUtility(label="+union(AT)", dataset=bootMergedWithPairs %where% (Affy_tag_accepts | Jetset_accepts & Geneannot_specificity_accepts & Geneannot_sensitivity_accepts & Affy_Grade_accepts)),
  expectedUtility(label="+union(AG)", dataset=bootMergedWithPairs %where% (Affy_Grade_accepts | Jetset_accepts & Geneannot_specificity_accepts & Geneannot_sensitivity_accepts & Affy_Grade_accepts)),
  expectedUtility(label="+union(E)", dataset=bootMergedWithPairs %where% (ENCODE_accepts | Jetset_accepts & Geneannot_specificity_accepts & Geneannot_sensitivity_accepts & Affy_Grade_accepts)),
  expectedUtility(label="+union(GQ)", dataset=bootMergedWithPairs %where% (Geneannot_quality_accepts | Jetset_accepts & Geneannot_specificity_accepts & Geneannot_sensitivity_accepts & Affy_Grade_accepts)),
  expectedUtility(label="+union(GSEN)", dataset=bootMergedWithPairs %where% (Geneannot_sensitivity_accepts | Jetset_accepts & Geneannot_specificity_accepts & Geneannot_sensitivity_accepts & Affy_Grade_accepts)),
  expectedUtility(label="+union(GSPE)", dataset=bootMergedWithPairs %where% (Geneannot_specificity_accepts | Jetset_accepts & Geneannot_specificity_accepts & Geneannot_sensitivity_accepts & Affy_Grade_accepts))

  )

resultLevel6TEU = rbind(
  expectedUtility(label="+int(PD)", dataset=bootMergedWithPairs %where% (PDBA_accepts & Jetset_accepts & Geneannot_specificity_accepts & Geneannot_sensitivity_accepts & Affy_Grade_accepts & Geneannot_quality_accepts)),
  expectedUtility(label="+int(M)", dataset=bootMergedWithPairs %where% (Masker_accepts & Jetset_accepts & Geneannot_specificity_accepts & Geneannot_sensitivity_accepts & Affy_Grade_accepts & Geneannot_quality_accepts)),
  expectedUtility(label="+int(AT)", dataset=bootMergedWithPairs %where% (Affy_tag_accepts & Jetset_accepts & Geneannot_specificity_accepts & Geneannot_sensitivity_accepts & Affy_Grade_accepts & Geneannot_quality_accepts)),
  expectedUtility(label="+int(AG)", dataset=bootMergedWithPairs %where% (Affy_Grade_accepts & Jetset_accepts & Geneannot_specificity_accepts & Geneannot_sensitivity_accepts & Affy_Grade_accepts & Geneannot_quality_accepts)),
  expectedUtility(label="+int(E)", dataset=bootMergedWithPairs %where% (ENCODE_accepts & Jetset_accepts & Geneannot_specificity_accepts & Geneannot_sensitivity_accepts & Affy_Grade_accepts & Geneannot_quality_accepts)),
  expectedUtility(label="+int(GQ)", dataset=bootMergedWithPairs %where% (Geneannot_quality_accepts & Jetset_accepts & Geneannot_specificity_accepts & Geneannot_sensitivity_accepts & Affy_Grade_accepts & Geneannot_quality_accepts)),
  expectedUtility(label="+int(GSEN)", dataset=bootMergedWithPairs %where% (Geneannot_sensitivity_accepts & Jetset_accepts & Geneannot_specificity_accepts & Geneannot_sensitivity_accepts & Affy_Grade_accepts & Geneannot_quality_accepts)),
  expectedUtility(label="+int(GSPE)", dataset=bootMergedWithPairs %where% (Geneannot_specificity_accepts & Jetset_accepts & Geneannot_specificity_accepts & Geneannot_sensitivity_accepts & Affy_Grade_accepts & Geneannot_quality_accepts)),
  expectedUtility(label="+union(PD)", dataset=bootMergedWithPairs %where% (PDBA_accepts | Jetset_accepts & Geneannot_specificity_accepts & Geneannot_sensitivity_accepts & Affy_Grade_accepts & Geneannot_quality_accepts)),
  expectedUtility(label="+union(M)", dataset=bootMergedWithPairs %where% (Masker_accepts | Jetset_accepts & Geneannot_specificity_accepts & Geneannot_sensitivity_accepts & Affy_Grade_accepts & Geneannot_quality_accepts)),
  expectedUtility(label="+union(AT)", dataset=bootMergedWithPairs %where% (Affy_tag_accepts | Jetset_accepts & Geneannot_specificity_accepts & Geneannot_sensitivity_accepts & Affy_Grade_accepts & Geneannot_quality_accepts)),
  expectedUtility(label="+union(AG)", dataset=bootMergedWithPairs %where% (Affy_Grade_accepts | Jetset_accepts & Geneannot_specificity_accepts & Geneannot_sensitivity_accepts & Affy_Grade_accepts & Geneannot_quality_accepts)),
  expectedUtility(label="+union(E)", dataset=bootMergedWithPairs %where% (ENCODE_accepts | Jetset_accepts & Geneannot_specificity_accepts & Geneannot_sensitivity_accepts & Affy_Grade_accepts & Geneannot_quality_accepts)),
  expectedUtility(label="+union(GQ)", dataset=bootMergedWithPairs %where% (Geneannot_quality_accepts | Jetset_accepts & Geneannot_specificity_accepts & Geneannot_sensitivity_accepts & Affy_Grade_accepts & Geneannot_quality_accepts)),
  expectedUtility(label="+union(GSEN)", dataset=bootMergedWithPairs %where% (Geneannot_sensitivity_accepts | Jetset_accepts & Geneannot_specificity_accepts & Geneannot_sensitivity_accepts & Affy_Grade_accepts & Geneannot_quality_accepts)),
  expectedUtility(label="+union(GSPE)", dataset=bootMergedWithPairs %where% (Geneannot_specificity_accepts | Jetset_accepts & Geneannot_specificity_accepts & Geneannot_sensitivity_accepts & Affy_Grade_accepts & Geneannot_quality_accepts))

  )

plot(resultAllTEU$PrPlus, resultAllTEU$nPairs, pch = "", main = "Greedy Forward Selection Across IDF Resources", xlab = "Estimated Proportion of Concordance", ylab = "Number of Pairs")
lines(resultTEUOptimize$PrPlus, resultTEUOptimize$nPairs)

points(resultLevel4TEU$PrPlus, resultLevel4TEU$nPairs, col="green", pch = "4")
points(resultLevel3TEU$PrPlus, resultLevel3TEU$nPairs, col="red", pch = "3")
points(resultLevel2TEU$PrPlus, resultLevel2TEU$nPairs, col="blue",pch = "2")
points(resultLevel1TEU$PrPlus, resultLevel1TEU$nPairs, col="black",pch = "1")

#points(resultLevel5TEU$PrPlus, resultLevel5TEU$nPairs, col="purple", pch = "5")
#points(resultLevel6TEU$PrPlus, resultLevel6TEU$nPairs, col="red", pch = "6")
lines(resultJGreedy$PrPlus, resultJGreedy$nPairs, lty = 2)
text(resultLevel1TEU$PrPlus, resultLevel1TEU$nPairs, labels = resultLevel1$label, pos = 1)
legend(locator(1),c("Total EU", "Mean EU"), lty = c(1,2))
###6leveltree###legend(locator(1), c("Tree Level 1","Tree Level 2","Tree Level 3","Tree Level 4","Tree Level 5","Tree Level 6"), text.col= c("black","blue","green","orange","purple","red"))

legend(locator(1), c("Tree Level 1","Tree Level 2","Tree Level 3","Tree Level 4"), text.col= c("black","blue","red","green"))

####Useall
draw.circle(.303,887,.015,border="gray")
draw.circle(.303,887,.0125,border="gray",lty=2)

####LEVEL1
draw.circle(.396,434,.015,border="black")
draw.circle(.396,434,.0125,border="black",lty=2)


####LEVEL2
draw.circle(.455,387,.015,border="blue")
draw.circle(.464,357,.0125,border="blue",lty=2)

####Level3
draw.circle(.494,319,.015,border="red")
draw.circle(.497,297,.0125,border="red",lty=2)






#plot(1:5,seq(1,10,length=5),type="n",xlab="",ylab="",main="Test draw.circle")
#draw.circle(2,4,c(1,0.66,0.33),border="purple",
#            col=c("#ff00ff","#ff77ff","#ffccff"),lty=1,lwd=1)
#draw.circle(2.5,8,0.6,border="red",lty=3,lwd=3)
#draw.circle(4,3,0.7,border="green",lty=1,lwd=1)
#draw.circle(3.5,7,0.8,border="blue",lty=2,lwd=2)




