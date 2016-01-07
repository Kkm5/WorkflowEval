fit2clusters.RNAseq<-function(Y, Ysigsq,
                              piStart = c(0.5, 0.5),
                              VStart = c(0.1,0.1),
                              psiStart = c(0,0.1),
                              NinnerLoop = 1,
                              nReps=500,
                              psi0Constraint,
                              V0Constraint,
                              sameV=FALSE,
                              estimatesOnly=TRUE,
                              printMe = TRUE,
                              plotMe = TRUE,
                              testMe=FALSE,
                              Ntest = 5000,
                              seed) {
  ### EM algorithm for 2 clusters,
  ### with constraints on the cluster means and variances, and known data variances
  if(testMe) {
    if(missing(seed)) .Random.seed <<- Random.seed.save
    else if(!is.na(seed)) .Random.seed <<- seed
    #  NA ==>  a new dataset.
    simPsi = c(0, 0.4)  ##
    simPi = c(2/3, 1/3)
    simData = data.frame(G = 1+rbinom(Ntest, 1, simPi[2]))
    simV = c(0.05^2, 0.05^2)
    simData$Ysigsq = rgamma(Ntest, 10, 400)
    simData$sd = sqrt(simV[simData$G] +simData$Ysigsq)
    simData = within(simData, Y <- simPsi[G] + rnorm(Ntest)*sqrt(simV[G]) + rnorm(Ntest)*sqrt(Ysigsq))
    print(summary(simData$Y))
    Y = simData$Y
    Ysigsq = simData$Ysigsq
  }
  ###############  Begining of EM algorithm  ##############
  piStar = piStart
  VStar = VStart
  psiStar = psiStart
  stopMe = FALSE
  iRep = 0
  while(1) {
    iRep = iRep + 1
    #    catn(", ", missing(V0Constraint))
    if(!missing(V0Constraint))
      VStar[1] = V0Constraint
    if(!missing(psi0Constraint))
      psiStar[1] = psi0Constraint
    #    print(psiStar)
    piStarOdds = piStar[2]/piStar[1]
    piStarOddsGK = piStarOdds *
      dnorm(Y, psiStar[2], sqrt(VStar[2] + Ysigsq)) /
      dnorm(Y, psiStar[1], sqrt(VStar[1] + Ysigsq))
    piStarGK = cbind(1/(1+piStarOddsGK), piStarOddsGK/(1+piStarOddsGK))
    EstarN = apply(piStarGK, 2, sum)
    piStar = apply(piStarGK, 2, mean)
    psiHat = psiStar
    VHat = VStar
    for(iRepInner in 1:NinnerLoop) {
      varHatTotal = colSums(outer(Y, psiHat, "-")^2 * piStarGK)
      sigsqTotal  = Ysigsq %*% piStarGK
      VHat = pmax(0, varHatTotal - sigsqTotal) / EstarN
      if(!missing(V0Constraint))
        VHat[1] = V0Constraint
      psiHat = colSums( Y %*% (piStarGK
                               / outer(Ysigsq, VHat, "+"))) /
        colSums( piStarGK
                 / outer(Ysigsq, VHat, "+"))
      if(!missing(psi0Constraint))
        psiHat[1] = psi0Constraint
      if(sameV)
        VHat[1] = VHat[2] = mean(VHat[1], VHat[2])
      if(max(abs(psiHat-psiStar), abs(VHat-VStar)) < 1e-7)
        stopMe = TRUE;
      psiHat -> psiStar
      VHat -> VStar
    }
    if(iRep >= nReps) stopMe = TRUE
    if(stopMe) break
  }
  cat(iRep, ifelse(iRep==nReps, ".  Loop exhausted.", ".  Converged."), "\n")
  
  if(plotMe) {
    options(echo=F)
    plot(col="blue", type = "l",main ="Mixture density", Ytemp<-seq(-1,1,length=100),
         xlab="Pearson's correlation", cex.lab = 1.5, ylab="Density", sub = "RPPA fold change vs RNASEQ workflow v1 and v2",
         piStar[1]*dnorm(Ytemp, psiStar[1], sqrt(VStar[1] + mean(Ysigsq)))
         +
           piStar[2]*dnorm(Ytemp, psiStar[2], sqrt(VStar[2]  + mean(Ysigsq)))
    )
    for(g in 1:2) lines(col="blue", lty=2, Ytemp<-seq(-1,1,length=100),
                        piStar[g]*dnorm(Ytemp, psiStar[g], sqrt(VStar[g]  + mean(Ysigsq))))
    lines(density(Y), lwd=2, col="black")
    abline(v = 0)
    abline(v = 0.38352)
    ###  Should we make a better choice than the means of the Ysigsq?
    if(testMe) lines(col="red", Ytemp<-seq(-1,1,length=100),
                     piStar[1]*dnorm(Ytemp, simPsi[1], sqrt(simV[1] + mean(simData$Ysigsq)))
                     +
                       piStar[2]*dnorm(Ytemp, simPsi[2], sqrt(simV[2]  + mean(simData$Ysigsq)))
    )
    legend(x=par("usr")[1], y=par("usr")[4],
           legend=c(ifelse(testMe, "truth", ""),
                    "data smooth", "estimate", " x or 0 component", " + component"),
           col=c("red", "black", "blue", "blue", "blue"),
           lty=c(ifelse(testMe, 1,0),1,1,2,2),
           lwd=c(1,2,1,1,1)
    )
    options(echo=T)
  }
  estimates = c(pi1=piStar[2], psi0=psiHat[1],
                psi1=psiHat[2], Var0=VStar[1], Var1=VStar[2])
  
  posteriorOdds =
    piStar[2]*dnorm(Y, psiHat[2], sqrt(VStar[2] + Ysigsq)) /
    piStar[1]/dnorm(Y, psiHat[1], sqrt(VStar[1] + Ysigsq))
  postProb = posteriorOdds/(1+posteriorOdds)
  postProbVar = Ysigsq * (postProb*(1-postProb))^2 *
    ((Y-psiHat[1])/(VStar[1]+Ysigsq) - (Y-psiHat[2])/(VStar[2]+Ysigsq))^2
  
  if(testMe) {
    simTruth = c(pi1=simPi[2], psi0=simPsi[1],
                 psi1=simPsi[2], Var0=simV[1], Var1=simV[2])
    estimates = data.frame(row.names=c("true","estimated"),
                           rbind(simTruth, estimates))
  }
  if(estimatesOnly) return(estimates)
  else {
    attr(x=posteriorOdds, which="estimates") = estimates
    return(data.frame(posteriorOdds,postProbVar))
  }
  
  
}
