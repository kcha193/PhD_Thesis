setwd("C:/Documents and Settings/krug001/My Documents/My Research/Proteomics/Paper 1/Decomposition example")
rm(list=ls())

# Source file containing functions to compute ESS

source("Computation of ESS v2.r")

# 4-plex iTRAQ design with 4 runs, 4 samples (independent biological reps) per run
# 4 Treatments (Ref, A, B, C), each replicated 4 times

# Define design parameters

  # Associated with physical material

  nRuns <- 4
  nSubj <- 13

  # Associated with treatments

  nTrt  <- 4
  trtLevels <- c("Ref", "A", "B", "C")
  v     <- nTrt
  nReps <- 4

  nObs  <- nReps * v


# Set up design

myDesign <- local({

   Run <- factor(rep(1:nRuns, each=nTrt))
   Tag <- factor(rep(114:117, times=nRuns))
   Treatment <- factor( c(circ(trtLevels, 4)), levels=trtLevels)
   Control <- factor(ifelse(Treatment=="Ref", "Ref", "Other"))
   Subject <- rep(0, nObs)
   Subject[Treatment != "Ref"] <- 1:12
   Subject[!Subject] <- 13
   Subject <- factor(Subject)
   
   data.frame(Subject, Control, Treatment, Run, Tag)

})
myDesign

# Phase 1 design

Z  <- makeBlkDesMatrix(design.df=myDesign, blkOrder=c(1,2,3,4))
Pb <- makeBlockProjectors(Z)
Pb <- makeBlockProjectors2(Z)

myEffects <- factorIncidenceMatrix(trtCols=2:3, factorNames=names(myDesign)[2:3])
   # X <- makeTreatDesignMatrix(design.df=myDesign, trtCols=2:3, effectsMatrix=myEffects[-3,])
   # T <- makeTreatProjectors(design.df=myDesign, trtCols=2:3, effectsMatrix=myEffects[-3,])

X <- NULL
X$Control <- matrix(0, nrow=nTrt, ncol=nlevels(myDesign$Control), dimnames=list(levels(myDesign$Treatment), levels(myDesign$Control)))
X$Control[cbind(1:4, c(2,1,1,1))] <- 1
X$Treatment <- diag(4)
dimnames(X$Treatment) <- list(levels(myDesign$Treatment), levels(myDesign$Treatment))

T <- NULL
T$Control <- projMat(X$Control) - K(4)
T$Treatment <- projMat(X$Treatment) - T$Control - K(4)

N <- matrix(0, nrow=nObs, ncol=nlevels(myDesign$Treatment), dimnames=list(1:nObs, levels(myDesign$Treatment)))
N[cbind(1:nObs, match(myDesign$Treatment, levels(myDesign$Treatment)))] <- 1

# Get treatment information matrices in each stratum
A <- lapply(Pb, function(x) lapply(T, function(y) InfMat(C=x, N=N, T=y)))

VCs <- getVCs(Z, Pb, N, T)
lapply(VCs, function(x) round(x, 4))

# Phase 2 ignoring a) peptide derivatization and b) tags

#-- Fit Runs then Subjects, exclude intrablock stratum
Z  <- makeBlkDesMatrix(design.df=myDesign, blkOrder=c(1,4))
Pb <- makeBlockProjectors(Z)
# Get treatment information matrices in each stratum
A <- lapply(Pb, function(x) lapply(T, function(y) InfMat(C=x, N=N, T=y)))

VCs <- getVCs2(Z, Pb, N, T)
lapply(VCs, function(x) round(x, 4))

# Let's try using an artificial stratum construct, called "column"

myDesign$Column <- factor( as.numeric(as.character(myDesign$Tag)) - 113 )
data.frame(names(myDesign),1:ncol(myDesign))
Z  <- makeBlkDesMatrix(design.df=myDesign, blkOrder=c(1,4))
Z2 <- makeBlkDesMatrix(design.df=myDesign, blkOrder=c(6,4))
Pb <- local({

   run <- projMat(Z2$Run) - K(nObs)
   col <- projMat(Z2$Column) - K(nObs)
   run.col <- projMat(Z2$e) - run - col - K(nObs)

   list(e=run.col, Column=col, Run=run)

})

# Get treatment information matrices in each stratum
A <- lapply(Pb, function(x) lapply(T, function(y) InfMat(C=x, N=N, T=y)))

VCs <- getVCs2(Z, Pb, N, T)
lapply(VCs, function(x) round(x, 4))

# Phase 2 including tags, but still ignoring peptide derivatization

#-- Now include the tag effect

myEffects <- factorIncidenceMatrix(trtCols=c(2:3,5), factorNames=names(myDesign)[c(2:3,5)])
nTag <- 4
nTrtTag <- nTrt * nTag
trtTagNames <- with(myDesign, paste(rep(levels(Tag), each=nTrt), levels(Treatment), sep="."))

X <- NULL
nCol <- nlevels(myDesign$Control) * nTags
X$Control <- matrix(0, nrow=nTrtTag, ncol=nCol, dimnames=list(trtTagNames, levels(myDesign$Control)))
X$Control[cbind(1:4, c(2,1,1,1))] <- 1
X$Treatment <- diag(8)
dimnames(X$Treatment) <- list(levels(myDesign$Treatment), levels(myDesign$Treatment))

T <- NULL
T$Control <- projMat(X$Control) - K(4)
T$Treatment <- projMat(X$Treatment) - T$Control - K(4)

# Get treatment information matrices in each stratum
A <- lapply(Pb, function(x) lapply(T, function(y) InfMat(C=x, N=N, T=y)))

VCs <- getVCs2(Z, Pb, N, T)
lapply(VCs, function(x) round(x, 4))


#-- Fit Subjects then Runs, exclude intrablock stratum
Z  <- makeBlkDesMatrix(design.df=myDesign, blkOrder=c(4,1))
Pb <- makeBlockProjectors(Z)
Pb <- Pb[-1]
# Get treatment information matrices in each stratum
A <- lapply(Pb, function(x) lapply(T, function(y) InfMat(C=x, N=N, T=y)))

VCs <- getVCs2(Z, Pb, N, T)
lapply(VCs, function(x) round(x, 4))

#-- Fit Subjects, exclude Runs
Z  <- makeBlkDesMatrix(design.df=myDesign, blkOrder=c(1,4))
Pb <- makeBlockProjectors(Z)
Pb <- Pb[-1]
# Get treatment information matrices in each stratum
A <- lapply(Pb, function(x) lapply(T, function(y) InfMat(C=x, N=N, T=y)))

VCs <- getVCs2(Z, Pb, N, T)
lapply(VCs, function(x) round(x, 4))



#-- Fit Subjects then Runs
Z  <- makeBlkDesMatrix(design.df=myDesign, blkOrder=c(1,4))
names(Z)
Pb <- local({

   run  <- projMat(Z$Run)-K(nObs)
   subj <- projMat(Z$Subject)-K(nObs)
   run.subj <- projMat(Z$Subject)-K(nObs) - run - subj - K(nObs)
   list(Run=run, Subject=subj, Run.Subject=run.subj)

})
# Get treatment information matrices in each stratum
A <- lapply(Pb, function(x) lapply(T, function(y) InfMat(C=x, N=N, T=y)))

VCs <- getVCs2(Z, Pb, N, T)
lapply(VCs, function(x) round(x, 4))

