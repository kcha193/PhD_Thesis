setwd("C:/Documents and Settings/krug001/My Documents/My Research/Proteomics/Paper 1/Decomposition example")
rm(list=ls())

# Source file containing functions to compute ESS

source("Computation of ESS v2.r")

# 4-plex iTRAQ design with 4 runs, 4 samples (independent biological reps) per run
# 4 Treatments (A, B, C, D), each replicated 4 times
# 3 peptide Sequences per protein

# Define design parameters

  # Associated with physical material

  nRuns <- 4
  nSubj <- 16
  nPeps <- 3

  # Associated with treatments

  nTrt  <- 4
  nSeq  <- 3
  v     <- nTrt * nSeq
  nReps <- 4

  nObs  <- nReps * v

# Set up design

myDesign <- local({

   Run <- factor(rep(1:nRuns, each=nTrt*nPeps))
   Peptide <- factor(rep(paste("P", 1:3, sep=""), times=nSubj))
   Subject <- factor(rep(1:nSubj, each=nPeps))
   Treatment <- factor(rep(LETTERS[1:nTrt], each=nPeps, times=nRuns))
   Sequence <- factor(rep(paste("Q", 1:3, sep=""), times=nSubj))
   
   data.frame(Run, Peptide, Subject, Treatment, Sequence)

})
myDesign

Z <- makeBlkDesMatrix(design.df=myDesign, blkOrder=c(3,1))
Pb <- makeBlockProjectors(Z)
myEffects <- factorIncidenceMatrix(trtCols=c(4,5), factorNames=names(myDesign)[c(4,5)])
X <- makeTreatDesignMatrix(design.df=myDesign, trtCols=c(4,5), effectsMatrix=myEffects)
T <- makeTreatProjectors(design.df=myDesign, trtCols=c(4,5), effectsMatrix=myEffects)
N <- getIncidenceMatrix(nRows=nObs, design.df=myDesign, trtCols=4:5)

# Get treatment information matrices in each stratum
A <- lapply(Pb, function(x) lapply(T, function(y) InfMat(C=x, N=N, T=y)))

VCs <- getVCs(Z, Pb, N, T)
lapply(VCs, function(x) round(x, 4))



