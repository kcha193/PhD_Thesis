getVCs <- function(blkDesMat, blkProj, incMat, trtProj, nTrtLevel)
{

   blkDesMat = Z; blkProj = Pb1; incMat = N; trtProj = Treat; nTrtLevel = c(nTre, nPos, nlevels(designPhase1$Inter))
   nStrata <- length(blkProj)
   nEffects <- length(trtProj)
   nBlkFactors <- length(blkDesMat)

   nObs <- nrow(blkDesMat$e)

   ginvA <- lapply(blkProj, function(x) lapply(trtProj, function(y) invInfMat(C=x, N=incMat, T=y)))

   effFactor <- t(sapply(blkProj, function(x) sapply(trtProj, function(y) eigenValue(C=x, N=incMat, T=y))))
   
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
   V <- lapply(Z1, function(x) x %*% t(x))

   # And finally get coefficients of e, S and R VCs for each treatment effect in each stratum 
   VC <- rep( list(NULL), length(V))
   names(VC) <- names(V)
   for(i in 1:nStrata){
      tmp <- matrix(0, nrow=nEffects+1, ncol=length(V),
            dimnames=list(names(PNTginvATNP[[1]]), names(V)))
      for(j in 1:(nEffects + 1)){
         for(k in 1:length(V)){
           tmp[j,k] <- tr(PNTginvATNP[[i]][[j]] %*% V[[k]])
         }
                  }
      VC[[i]] <- tmp[-which(apply(tmp, 1, function(x) all(x < 0.000001))),]
   }
   VC <- VC[!unlist(lapply(VC, function(x) is.null(x)))]  # Need this when want fewere strata than there are block structuresV
   VC <- VC[sort(1:length(VC), decreasing=TRUE)]
   return(VC)
}


bdeff <- function(blocks, varieties) {
         blocks <- as.factor(blocks)             # minor safety move
         b <- length(levels(blocks))
         varieties <- as.factor(varieties)       # minor safety move
         v <- length(levels(varieties))
         K <- as.vector(table(blocks))           # remove dim attr
         R <- as.vector(table(varieties))        # remove dim attr
         N <- table(blocks, varieties)
         A <- 1/sqrt(K) * N * rep(1/sqrt(R), rep(b, v))
         sv <- svd(A)
         list(eff=1 - sv$d^2, blockcv=sv$u, varietycv=sv$v)
     }
     
getVCs2 <- function(blkDesMat, blkProj, incMat, trtProj, tInc, nTrtLevel)
{
   blkDesMat = Z; blkProj = PbPhase1; incMat = N; trtProj = Treat; tInc = T; 
   nTrtLevel = c(nTre, nMet, nlevels(designPhase2$Inter))
     
   nEffects <- length(trtProj)
   nBlkFactors <- length(blkDesMat)

   nObs <- nrow(blkDesMat$e)

   effFactor <- lapply(blkProj,  function(x) lapply(x, function(z)
   fractions(sapply(trtProj, function(y) eigenValue(C=z, N=incMat, T=y)))))
   
   #a function that computes the sum of the matrices for computing the PNTginvATNP
   matrixSum =
     function(matrixList){
      matrixSum = matrix(0, nrow = dim(matrixList[[1]])[1], ncol = dim(matrixList[[1]])[2])
      for(i in 1:length(matrixList)){      
        matrixSum = matrixSum + matrixList[[i]]
      }      
      return(matrixSum)   
   }
   
  
   #pre- and post-multiply NTginvATN by block projection matrices
   PNTginvATNP <- lapply(blkProj, function(x) lapply(x, function(z){
      blkProkMat <- lapply(trtProj, function(y) t(z) %*% incMat %*% t(y) %*% invInfMat(C=z, N=incMat, T=y) %*% y %*% t(incMat) %*% z)      
      blkProkMat$Error = (t(z) %*% z) - matrixSum(blkProkMat)
      return(blkProkMat)      
   }))
  
   
   lapply(PNTginvATNP, function(x) lapply(x, function(y) lapply(y,  tr)))
   
   #Now construct variance matrices
   V <- lapply(blkDesMat, function(x) x %*% t(x))
   
   #And finally get coefficients of e, S and R VCs for each treatment effect in each stratum
   VC <- PNTginvATNP #to keep the structure
   contrT <- PNTginvATNP                        
   
   for(i in 1:length(PNTginvATNP)){
      for(k in 1:length(PNTginvATNP[[i]])){
        tmp <- matrix(0, nrow=nEffects+1, ncol=nBlkFactors, dimnames=list(names(PNTginvATNP[[i]][[k]]), names(V)))
        for(j in 1:(nEffects + 1)){
         for(z in 1:nBlkFactors){
           tmp[j,z] <- tr(PNTginvATNP[[i]][[k]][[j]] %*% V[[z]])
         }
         
         contrT[[i]][[k]][[j]] <- t(tInc) %*% PNTginvATNP[[i]][[k]][[j]] %*% tInc 
         
         #if(j != (nEffects + 1) && tmp[j,1] > 0.0000001){
         #  tmp[j, nBlkFactors+j] <- (nObs/nTrtLevel[j]) * (effFactor[i,j]/tmp[j,1])
         #}
        }
      VC[[i]][[k]] <- tmp[-which(apply(tmp, 1, function(x) all(x < 0.000001))),-which(apply(tmp, 2, function(x) all(x < 0.000001)))]   
      }
   }
  

   #VC <- VC[!unlist(lapply(VC, function(x) is.null(x)))]  # Need this when want fewere strata than there are block structuresV
   VC <- VC[sort(1:length(VC), decreasing=TRUE)]
   
   contrT <- contrT[sort(1:length(contrT), decreasing=TRUE)]
      
   return(VC)
}
