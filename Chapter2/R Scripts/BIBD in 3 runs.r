setwd("C:/Documents and Settings/krug001/My Documents/My Research/Proteomics/Paper 1/Decomposition example")
rm(list=ls())

# Source file containing functions to compute ESS

source("Computation of ESS v2.r")

# 4-plex iTRAQ design
# $ reps of 3 Treatments (A, B, C) arranged in a BIBD in 3 runs

# Define design parameters

  # Associated with physical material

  nRuns <- 3
  nSubj <- 12

  # Associated with treatments

  nTrt  <- 3
  nTags <- 4
  trtLevels <- c("A", "B", "C")
  v     <- nTrt * nTags
  nReps <- 4

  nObs  <- nReps * nTrt

# Set up design

myDesign <- local({

   Run <- factor(rep(1:nRuns, each=nTags))
   Tag <- factor(rep(114:117, times=nRuns))
   Treatment <- factor(c( t(cbind(circ(LETTERS[1:3], 3), circ(LETTERS[1:3], 3))) ), levels=LETTERS[1:3])
   #Subject <- factor(makeTechRepIDs(Treatment))
   Subject <- factor(1:12)

   data.frame(Run, Subject, Tag, Treatment)

})
myDesign

Z <- makeBlkDesMatrix(design.df=myDesign, blkOrder=2)
Pb <- makeBlockProjectors(Z)
myEffects <- factorIncidenceMatrix(trtCols=4, factorNames=names(myDesign)[4])
X <- makeTreatDesignMatrix(design.df=myDesign, trtCols=4, effectsMatrix=myEffects)
T <- makeTreatProjectors(design.df=myDesign, trtCols=4, effectsMatrix=myEffects)
N <- getIncidenceMatrix(nRows=nObs, design.df=myDesign, trtCols=4)
# Get treatment information matrices in each stratum
A <- lapply(Pb, function(x) lapply(T, function(y) InfMat(C=x, N=N, T=y)))
VCs <- getVCs(Z, Pb, N, T)
lapply(VCs, function(x) round(x, 4))

# Add peptides

  # Associated with physical material

  nRuns <- 3
  nSubj <- 12
  nPeps <- 3
  #nSpec <- 2

  # Associated with treatments

  nTrt  <- 3
  nTags <- 4
  nSeq  <- 3
  v     <- nTrt * nSeq
  #v     <- nTrt * nTags * nSeq
  nReps <- 4

  nObs  <- nReps * nTrt * nSeq
#  nObs  <- nReps * nTrt * nSeq * nSpec

# Set up design

myDesign <- local({

   Run <- factor(rep(1:nRuns, each=nTags))
   Subject <- factor(1:12)
   Tag <- factor(rep(114:117, times=nRuns))
   Treatment <- factor(c( t(cbind(circ(LETTERS[1:3], 3), circ(LETTERS[1:3], 1))) ), levels=LETTERS[1:3])
   #Subject <- factor(makeTechRepIDs(Treatment)) # Use this when get function to check for nested structure

   df <- data.frame(Run, Subject, Treatment, Tag)
   df <- df[rep(1:nrow(df), each=nPeps), ]

   Peptide <- factor(rep(paste("P", 1:3, sep=""), times=nSubj))
   Sequence <- factor(rep(paste("Q", 1:3, sep=""), times=nSubj))

   df <- data.frame(df, Peptide, Sequence)
   rownames(df) <- 1:nrow(df)
   return(df)

})
myDesign

Z <- makeBlkDesMatrix(design.df=myDesign, blkOrder=2)
Pb <- makeBlockProjectors(Z)
myTrtCols <- c(3,6)
myEffects <- factorIncidenceMatrix(trtCols=myTrtCols, factorNames=names(myDesign)[myTrtCols])
X <- makeTreatDesignMatrix(design.df=myDesign, trtCols=myTrtCols, effectsMatrix=myEffects)
T <- makeTreatProjectors(design.df=myDesign, trtCols=myTrtCols, effectsMatrix=myEffects)
N <- getIncidenceMatrix(nRows=nObs, design.df=myDesign, trtCols=myTrtCols)
# Get treatment information matrices in each stratum
A <- lapply(Pb, function(x) lapply(T, function(y) InfMat(C=x, N=N, T=y)))
VCs <- getVCs(Z, Pb, N, T)
lapply(VCs, function(x) round(x, 4))

# After peptide derivatization and labelling, but before MudPIT analysis

Z <- makeBlkDesMatrix(design.df=myDesign, blkOrder=2)
Pb <- makeBlockProjectors(Z)
myTrtCols <- c(4,3,6)
myEffects <- factorIncidenceMatrix(trtCols=myTrtCols, factorNames=names(myDesign)[myTrtCols])
rmEffects <- c(3,7)
X <- makeTreatDesignMatrix(design.df=myDesign, trtCols=myTrtCols, effectsMatrix=myEffects[-rmEffects,])
T <- makeTreatProjectors(design.df=myDesign, trtCols=myTrtCols, effectsMatrix=myEffects[-rmEffects,])
N <- getIncidenceMatrix(nRows=nObs, design.df=myDesign, trtCols=myTrtCols)
# Get treatment information matrices in each stratum
A <- lapply(Pb, function(x) lapply(T, function(y) InfMat(C=x, N=N, T=y)))
VCs <- getVCs(Z, Pb, N, T)
lapply(VCs, function(x) round(x, 4))

# Now add 2 spectra per peptide

  # Associated with physical material

  nRuns <- 3
  nSubj <- 12
  nPeps <- 3
  nSpec <- 2

  # Associated with treatments

  nTrt  <- 3
  nTags <- 4
  nSeq  <- 3
  v     <- nTrt * nTags * nSeq
  nReps <- 4

  nObs  <- nReps * nTrt * nSeq * nSpec

# Set up design

myDesign <- local({

   df <- myDesign
   n <- nrow(df)
   df <- df[rep(1:n, each=nSpec), ]
   df$Peptide <- with(df, interaction(Subject, Peptide))
   df$Spectra <- factor(rep(1:nSpec, times=n))
   df$Spectra <- with(df, interaction(Peptide, Spectra))

   rownames(df) <- 1:nrow(df)
   return(df)

})
myDesign

# Full MudPIT

Z <- makeBlkDesMatrix(design.df=myDesign, blkOrder=c(5, 2,1))
Pb <- makeBlockProjectors(Z)
myTrtCols <- c(4,3,6)
myEffects <- factorIncidenceMatrix(trtCols=myTrtCols, factorNames=names(myDesign)[myTrtCols])
rmEffects <- c(3,7)
X <- makeTreatDesignMatrix(design.df=myDesign, trtCols=myTrtCols, effectsMatrix=myEffects[-rmEffects,])
T <- makeTreatProjectors(design.df=myDesign, trtCols=myTrtCols, effectsMatrix=myEffects[-rmEffects,])
N <- getIncidenceMatrix(nRows=nObs, design.df=myDesign, trtCols=myTrtCols)
# Get treatment information matrices in each stratum
A <- lapply(Pb, function(x) lapply(T, function(y) InfMat(C=x, N=N, T=y)))
VCs <- getVCs(Z, Pb, N, T)
lapply(VCs, function(x) round(x, 4))
