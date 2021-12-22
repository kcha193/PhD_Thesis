setwd("C:/Documents and Settings/krug001/My Documents/My Research/Proteomics/Paper 1/Decomposition example")
rm(list=ls())

# Source file containing functions to compute ESS

source("Computation of ESS v2.r")

# 4-plex iTRAQ design with 4 runs, 4 samples (independent biological reps) per run
# 4 Treatments (A, B, C, D), each replicated 4 times
# Samples within runs differentially labelled with 4-plex iTRAQ reagents
# 3 peptide Sequences per protein

# This design now includes 2 spectra for each peptide species

# Define design parameters

  # Associated with physical material

  nRuns <- 4
  nSubj <- 16
  nPeps <- 3
  nSpec <- 2

  # Associated with treatments

  nTrt  <- 4
  nTags <- 4
  nSeq  <- 3
  v     <- nTrt * nTags * nSeq
  nReps <- 4

  nObs  <- nReps * nTrt * nSeq * nSpec

# Set up design

circ <- function(x)
{
   l <- length(x)
   y <- 1:l
   for(i in 2:l){
      y <- (y %% l) +1
      x <- c(x, x[y])
   }
   return(x)
}

myDesign <- local({

   Run <- factor(rep(1:nRuns, each=nTrt*nPeps))
   Subject <- factor(rep(1:nSubj, each=nPeps))
   Peptide <- factor(rep(paste("P", 1:3, sep=""), times=nSubj))
   Peptide <- interaction(Subject, Peptide)
   Treatment <- factor(rep(LETTERS[1:nTrt], each=nPeps, times=nRuns))
   Tag <- factor(rep(circ(114:117), each=nPeps))
   Sequence <- factor(rep(paste("Q", 1:3, sep=""), times=nSubj))
   
   dat <- data.frame(Run, Peptide, Subject, Treatment, Tag, Sequence)
   n <- nrow(dat)
   dat <- dat[rep(1:n, each=nSpec),]
   dat$Spectra <- factor(rep(1:nSpec, times=n))
   dat$Spectra <- interaction(Peptide, dat$Spectra)
   
   dat

})
myDesign 

Z <- makeBlkDesMatrix(design.df=myDesign, blkOrder=c(2, 3,1))
Pb <- makeBlockProjectors(Z)
myEffects <- factorIncidenceMatrix(trtCols=c(5,4,6), factorNames=names(myDesign)[c(5,4,6)])
X <- makeTreatDesignMatrix(design.df=myDesign, trtCols=c(5,4,6), effectsMatrix=myEffects[-c(3,7),])
T <- makeTreatProjectors(design.df=myDesign, trtCols=c(5,4,6), effectsMatrix=myEffects[-c(3,7),])
N <- getIncidenceMatrix(nRows=nObs, design.df=myDesign, trtCols=4:6)

# Get treatment information matrices in each stratum
A <- lapply(Pb, function(x) lapply(T, function(y) InfMat(C=x, N=N, T=y)))
VCs <- getVCs(Z, Pb, N, T)
lapply(VCs, function(x) round(x, 4))

