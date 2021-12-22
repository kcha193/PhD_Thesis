
#setwd("C:/Documents and Settings/krug001/My Documents/My Research/Proteomics/Paper 1/Decomposition example")
#rm(list=ls())

# Set up functions for constructing projection matrices, the averaging matrix, and the trace

projMat <- function(X) X %*% solve(t(X) %*% X) %*% t(X)
K <- function(n) matrix(1/n, nrow=n, ncol=n)
J <- function(n) diag(n) - K(n)
tr <- function(X) sum(diag(X))
m1 <- function(n) matrix(rep(1,n), nrow = n, ncol = 1)


circ <- function(x, blkSize)
{
   l <- length(x)
   y <- 1:l
   for(i in 2:l){
      y <- (y %% l) +1
      x <- c(x, x[y])
   }
   circDes <- matrix(x, ncol=l, nrow=blkSize, byrow=TRUE)
   return(circDes)
}


makeTechRepIDs <- function(Factor)
{

   # Use this function when subsamples are taken from experimental units

   uniqueLevels <- make.unique( as.character(Factor) )
   tmp <- unlist(lapply( strsplit(uniqueLevels, "\\."), function(x) x[2]))
   tmp[is.na(tmp)] <- 0
   techRep <- as.numeric(tmp) + 1

   return(techRep)

}

isFactorNameNumeric <- function(levels) !as.logical( length(grep("[A-Z]|[a-z]", levels)) )

makeDesignMatrix <- function(nRows, design.df, col)
{
   factor <- design.df[,col]
   facName <- col
   nCols <- nlevels(factor)
   
   Z <- matrix(0, nrow=nRows, ncol=nCols)
   Z[cbind(1:nRows, match(c(factor), 1:nCols))] <- 1
   if(isFactorNameNumeric(levels(factor))) colNames <- paste(substr(facName,1,1), 1:nCols,sep="")
   else                                    colNames <- levels(factor)
   dimnames(Z) <- list(1:nRows, colNames)
   return(Z)
}

# Set up the block design matrices

makeBlkDesMatrix <- function(design.df, blkOrder){

   # design.df = data.frame containing design
   # blkOrder = order in which to set up block design matrices

   n <- length(blkOrder)
   nRows <- nrow(design.df)
   Z <- list(NULL)

   Z[[1]] <- diag(nrow(design.df))
   for(i in 2:(n+1)) Z[[i]] <- makeDesignMatrix(nRows=nRows, design.df=design.df, col=blkOrder[i-1]) 
   names(Z) <- c("e", blkOrder)

   return(Z)
}   

# Set up the Treatment design matrices (Richard)
makeTrtDesMatrix <- function(design.df, trtOrder, interactions = TRUE){
   # design.df = data.frame containing design
   # trtOrder = order in which to set up treatment design matrices
   # interactions= compute the matrices that represents the interaction
   
   n <- length(trtOrder)
   nRows <- nrow(design.df)
   
   G <- rep(1, nRows)
   T <- list()
   inter <-list()
   interCounter <- 1 
        
   for(i in 1:n){       
    nCols = nlevels(design.df[,trtOrder[i]])
    T[[i]] <- matrix(-1/nCols, ncol = nCols, nrow = nRows)
    T[[i]][cbind(1:nRows, match(design.df[,trtOrder[i]], 1:nCols))] <-  T[[i]][cbind(1:nRows, match(design.df[,trtOrder[i]], 1:nCols))]+1
    colnames(T[[i]]) <- paste(trtOrder[i], 1:nCols, sep ="")
    G <- cbind(G,T[[i]])     
    if( i > 1 && interactions==TRUE) { 
      for(j in 1:(i-1)){      
        Cols = sort(levels(interaction(levels(design.df[,trtOrder[j]]), levels(design.df[,trtOrder[i]]))))
        sepCols = strsplit(Cols, "\\.")
        inter[[interCounter]] <- sapply(sepCols, function(x) T[[j]][,as.numeric(x[1])] * T[[i]][,as.numeric(x[2])]) 
        colnames(inter[[interCounter]]) = sapply(sepCols, function(x) paste(trtOrder[j], x[1], trtOrder[i],x[2], sep = "")  )
        G <- cbind(G,inter[[interCounter]])
        interCounter <- interCounter + 1
      }
    }    
   }   
          
   return(G)
}   


# Make block projection matrices

makeBlockProjectors <- function(BlkDesignMatrixList)
{
   n <- nrow(BlkDesignMatrixList$e)
   Q <- lapply(BlkDesignMatrixList, function(z) projMat(z) - mK(n))
   Q <- Q[sort(1:length(Q), decreasing=TRUE)]
   cusumQ <- P <- NULL
   P[[1]] <- Q[[1]]
   cusumQ[[1]] <- matrix(0, nrow=n, ncol=n)
   for(i in 2:length(Q)){
      cusumQ[[i]] <- P[[i-1]] + cusumQ[[i-1]]
      P[[i]] <- Q[[i]] - cusumQ[[i]]
   }
   P <- P[sort(1:length(P), decreasing=TRUE)]
   names(P) <- names(BlkDesignMatrixList)
   return(P)

}

makeBlockProjectors2 <- function(BlkDesignMatrixList, initial = diag(nrow(BlkDesignMatrixList[[1]])) - K(nrow(BlkDesignMatrixList[[1]])))
{
  
   n <- nrow(BlkDesignMatrixList$e)
   Q <- lapply(BlkDesignMatrixList, function(z) projMat(z))
   Q <- Q[sort(1:length(Q), decreasing=TRUE)]
   elementToRemove <- numeric()
   
   cusumQ <- P <- NULL
   P[[1]] <- Q[[1]] %*% initial
      
   if(all(P[[1]] <0.000001)){
     elementToRemove = 1
     P[[1]] <- matrix(0, nrow = n, ncol = n)
   } else{
     P[[1]] <- (t(P[[1]]) %*% P[[1]])/svd(t(P[[1]]) %*% P[[1]])$d[1]
   }
   
   cusumQ[[1]] <- initial - P[[1]]

   for(i in 2:length(Q)){
      P[[i]] <- Q[[i]] %*% cusumQ[[i-1]]
        
      if(all(P[[i]] <0.000001)){
        elementToRemove <- c(elementToRemove, i)
        P[[i]] <- matrix(0, nrow = n, ncol = n)
      }else{
        P[[i]] <- (t(P[[i]]) %*% P[[i]])/svd(t(P[[i]]) %*% P[[i]])$d[1]

      }
      cusumQ[[i]] <- cusumQ[[i-1]] - P[[i]]
   }
   
   P <- P[sort(1:length(P), decreasing=TRUE)]
   names(P) <- names(BlkDesignMatrixList)

   if(length(elementToRemove) > 0){
    elementToRemove <- length(P) - elementToRemove + 1
    P <- P[-elementToRemove]
   }
   return(P)
}


factorIncidenceMatrix <- function(trtCols=trtOrder, factorNames)
{

   nFactors <- length(trtCols)
   facPointer <- diag(nFactors)
   facIncidMat <- facPointer[1,,drop=FALSE]
   for(i in 2:nFactors){
      if(nFactors==1) break
      facIncidMat <- rbind(facIncidMat, facPointer[i,])
      nRow <- nrow(facIncidMat)
      j <- 1
      while(j < nRow){
         facIncidMat <- rbind(facIncidMat, facIncidMat[j,] + facPointer[i,])
         j <- j + 1
      }
   }
   dimnames(facIncidMat) <- list(1:nrow(facIncidMat),factorNames)
   return(facIncidMat)
}

makeTreatProjectors <- function(design.df, trtCols, effectsMatrix)
{

   nLevels <- sapply(design.df[,trtCols], nlevels)
   effectNames <- apply(effectsMatrix, 1, function(x) paste(colnames(effectsMatrix)[as.logical(x)], collapse="."))
   nEffects <- nrow(effectsMatrix)
   X <- as.list(rep(1,nEffects))
   names(X) <- effectNames
   indMatrix <- function(x, n){
      if(x == 1) X <- J(n)
      else       X <- K(n)
      return(X)
   }
   for(i in 1:length(trtCols)){
      matList <- lapply(effectsMatrix[,i], function(y) indMatrix(y, nLevels[i]))
      for(j in 1:nEffects) X[[j]] <- X[[j]] %x% matList[[j]]
   }
   return(X)

}

makeTreatDesignMatrix <- function(design.df, trtCols, effectsMatrix)
{

   nLevels <- sapply(design.df[,trtCols], nlevels)
   effectNames <- apply(effectsMatrix, 1, function(x) paste(colnames(effectsMatrix)[as.logical(x)], collapse="."))
   nEffects <- nrow(effectsMatrix)
   X <- as.list(rep(1,nEffects))
   names(X) <- effectNames
   indMatrix <- function(x, n){
      if(x == 1) X <- diag(n)
      else       X <- rep(1, n)
      return(X)
   }
   for(i in 1:length(trtCols)){
      matList <- lapply(effectsMatrix[,i], function(y) indMatrix(y, nLevels[i]))
      for(j in 1:nEffects) X[[j]] <- X[[j]] %x% matList[[j]]
   }
   return(X)

}

getIncidenceMatrix <- function(nRows, design.df, trtCols)
{
   X <- list(NULL)
   nFacs <- length(trtCols)

   for(i in 1:nFacs) X[[i]] <- makeDesignMatrix(nRows=nObs, design.df=design.df, col=trtCols[i])
   names(X) <- names(design.df)[trtCols]

   N <- rep(1, nObs)
   index <- nFacs
   while(index){
      N <- as.matrix(data.frame(apply(data.frame(X[[index]]), 2, function(x) x * data.frame(N))))
      if(index == nFacs) colnames(N) <- colnames(X[[index]])
      index <- index - 1
   }
   rownames(N) <- 1:nrow(N)

   return(as.matrix(N))
}

# Compute information matrix

InfMat <- function(C,N,T){

   ei <- eigen(t(T) %*% t(N) %*% C %*% N %*% T)
   nn <- length(ei$values)
   L <- matrix(0, nrow=nn, ncol=nn)
   for(i in 1:(nn)) {
      if(ei$values[i]<1e-6)next
      L <- L + ei$values[i]*ei$vectors[,i]%*%t(ei$vectors[,i])
   }
   return(L)
}

# Compute g-inverse of information matrix

invInfMat <- function(C,N,T){

   ei <- eigen(t(T) %*% t(N) %*% C %*% N %*% T)
   nn <- length(ei$values)
   L <- matrix(0, nrow=nn, ncol=nn)
   for(i in 1:(nn)) {
    if( Re(ei$values[i]) <1e-6)next
    L <- L + (1/Re(ei$values[i]))*Re(ei$vectors[,i])%*%t(Re(ei$vectors[,i]))
   }
   return(L)
}

eigenValue <- function(C,N,T){
   fractions(svd(t(T) %*% t(N) %*% C %*% N %*% T)$d[1])
}


# And the coefficients of their associated VCs

getVCs <- function(blkDesMat, blkProj, incMat, trtProj)
{

   nStrata <- length(blkProj)
   nEffects <- length(trtProj)
   nBlkFactors <- length(blkDesMat)
   nObs <- nrow(blkDesMat$e)
   ginvA <- lapply(blkProj, function(x) lapply(trtProj, function(y) invInfMat(C=x, N=incMat, T=y)))
   
   # Now we want to get t(TN) %*% ginvA %*% TN for each treatment effect
   TN <- lapply(trtProj, function(x) x %*% t(incMat))
   # Intialize list
   NTginvATN <- rep( list( rep(list(NULL), nEffects) ), nStrata )
   for(i in 1:nStrata){
      for(j in 1:nEffects) NTginvATN[[i]][[j]] <- t(TN[[j]]) %*% ginvA[[i]][[j]] %*% TN[[j]]
      names(NTginvATN[[i]]) <- names(TN)
   }

   # Now, pre- and post-multiply NTginvATN by block projection matrices
   # Initialize list
   PNTginvATNP <- rep( list( rep(list(NULL), nEffects) ), nStrata )
   for(i in 1:nStrata){
      sumMat <- matrix(0, nrow=nObs, ncol=nObs)
      for(j in 1:nEffects){
         PNTginvATNP[[i]][[j]] <- t(blkProj[[i]]) %*% NTginvATN[[i]][[j]] %*% blkProj[[i]]
         sumMat <- PNTginvATNP[[i]][[j]] + sumMat
      }
      PNTginvATNP[[i]][[j+1]] <- blkProj[[i]] - sumMat
      names(PNTginvATNP[[i]]) <- c(names(NTginvATN[[i]]), "Error")
   }
   names(PNTginvATNP) <- names(NTginvATN)
   
   # Now construct variance matrices
   V <- lapply(blkDesMat, function(x) x %*% t(x))

   # And finally get coefficients of e, S and R VCs for each treatment effect in each stratum

   VC <- rep( list(NULL), nBlkFactors)
   names(VC) <- names(V)
   tmp <- matrix(NA, nrow=nEffects+1, ncol=nBlkFactors, dimnames=list(names(PNTginvATNP[[1]]), names(V)))
   for(i in 1:nStrata){
      for(j in 1:(nEffects + 1))
         for(k in 1:nBlkFactors) tmp[j,k] <- tr(PNTginvATNP[[i]][[j]] %*% V[[k]])
      VC[[i]] <- tmp
   }
   VC <- VC[!unlist(lapply(VC, function(x) is.null(x)))]  # Need this when want fewere strata than there are block structuresV
   VC <- VC[sort(1:length(VC), decreasing=TRUE)]
   return(VC)
}

getVCs2 <- function(blkDesMat, blkProj, incMat, trtProj)
{

   nStrata <- length(blkProj)
   nEffects <- length(trtProj)
   nBlkFactors <- length(blkDesMat)
   ginvA <- lapply(blkProj, function(x) lapply(trtProj, function(y) invInfMat(C=x, N=incMat, T=y)))
   
   # Now we want to get t(TN) %*% ginvA %*% TN for each treatment effect
   TN <- lapply(trtProj, function(x) x %*% t(incMat))
   # Intialize list
   NTginvATN <- rep( list( rep(list(NULL), nEffects) ), nStrata )
   for(i in 1:nStrata){
      for(j in 1:nEffects) NTginvATN[[i]][[j]] <- t(TN[[j]]) %*% ginvA[[i]][[j]] %*% TN[[j]]
      names(NTginvATN[[i]]) <- names(TN)
   }

   # Now, pre- and post-multiply NTginvATN by block projection matrices
   # Initialize list
   PNTginvATNP <- rep( list( rep(list(NULL), nEffects) ), nStrata )
   for(i in 1:nStrata){
      sumMat <- matrix(0, nrow=nObs, ncol=nObs)
      for(j in 1:nEffects){
         PNTginvATNP[[i]][[j]] <- t(blkProj[[i]]) %*% NTginvATN[[i]][[j]] %*% blkProj[[i]]
         sumMat <- PNTginvATNP[[i]][[j]] + sumMat
      }
      PNTginvATNP[[i]][[j+1]] <- blkProj[[i]] - sumMat
      names(PNTginvATNP[[i]]) <- c(names(NTginvATN[[i]]), "Error")
   }
   names(PNTginvATNP) <- names(NTginvATN)

   
   # Now construct variance matrices
   V <- lapply(blkDesMat, function(x) x %*% t(x))


   # And finally get coefficients of e, S and R VCs for each treatment effect in each stratum

   stratumVC <- rep( list(NULL), nStrata)
   names(stratumVC) <- names(blkProj)
   tmp <- matrix(NA, nrow=nEffects+1, ncol=nBlkFactors, dimnames=list(names(PNTginvATNP[[1]]), names(V)))
   for(i in 1:nStrata){
      for(j in 1:(nEffects + 1))
         for(k in 1:nBlkFactors) tmp[j,k] <- tr(PNTginvATNP[[i]][[j]] %*% V[[k]])
      stratumVC[[i]] <- tmp
   }

#   stratumVC <- stratumVC[!unlist(lapply(strata, function(x) is.null(x)))]  # Need this when want fewer strata than there are block structuresV
   stratumVC <- stratumVC[sort(1:length(stratumVC), decreasing=TRUE)]

   return(stratumVC)
}

is.nested <- function (factor1,factor2) 
  { 
    # only one positive number per line in the f1 * f2 crosstable 
    all(apply(table(factor1,factor2)>0,1,sum) == 1) 
  } 

are.crossed <- function (factor1,factor2) 
  { all(table(factor1,factor2) > 0 ) } 
