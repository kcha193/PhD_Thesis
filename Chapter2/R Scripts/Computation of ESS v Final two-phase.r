#2011/2/20 01:35:30
#Start from scratch
# input should be the design of data frame
#                 specify the random and fixed terms
#                  using the "terms" function to help to indicate cross and next relationships
#using the Nelder's method i.e. eignevalue decomposition to determine the fixed part of the EMS
# 
#The results should be in two: 1)ANOVA table structure with DF and coeffients of the variances compoenents of EMS for the random effects
#                              2)Treatment contribution matrices which may get complicated when the numbers of levels of treatment is large.
#2011/2/21 02:13:39
# next thing to do:
#         a) add the block strcutere for two-phase experiment: need to be carefull with the confounding between the block structures!!!
#         b) add the add the effeciency factors: careful with different contrasts!!


getVCs.twoPhase  <- 
function(design.df, random.terms1, random.terms2, fixed.terms){
  library(MASS)
 
#######################################################################################  
  #All the pre-initiated functions   
  projMat <- function(X) X %*% ginv(t(X) %*% X) %*% t(X) #projection matrix
  mK <- function(n) matrix(1/n, nrow=n, ncol=n)           #Averaging matrix
  mJ <- function(n) diag(n) - mK(n)                        #J matrix 
  tr <- function(X) sum(diag(X))                          #Trace operation
  m1 <- function(n) matrix(rep(1,n), nrow = n, ncol = 1)  #M1 matrix
  mI <- function(n) diag(n)                               #identity matrix
  
  
  #check for factor names 
  isFactorNameNumeric <- function(levels) !as.logical( length(grep("[A-Z]|[a-z]", levels)) )

  #make design matrix
  makeDesignMatrix <- function(nRows, design.df, col){    
    if(grepl(":", col)){ 
      factor <-as.factor( apply(design.df[,unlist(strsplit(col, ":"))], 1, function(x) paste(x, collapse =".")))
    }else{
      factor <- as.factor(design.df[,col])
    }  
        
    facName <- col
    nCols <- nlevels(factor)
   
    Z <- matrix(0, nrow=nRows, ncol=nCols)
    Z[cbind(1:nRows, match(c(factor), 1:nCols))] <- 1
    if(isFactorNameNumeric(levels(factor))){ 
      colNames <- paste(substr(facName,1,1), 1:nCols,sep="")
    }else{                                    
      colNames <- levels(factor)
    }
    
    dimnames(Z) <- list(1:nRows, colNames)
    return(Z)
  } 
 
 #make block design matrix 
 makeBlkDesMatrix <- function(design.df, blkTerm){

   # design.df = data.frame containing design
   # blkTerm = block terms

   n <- length(blkTerm)
   nRows <- nrow(design.df)
   Z <- list(NULL)
   Z[[1]] <- diag(nrow(design.df))
   
   for(i in 2:(n+1)){   
    Z[[i]] <- makeDesignMatrix(nRows=nRows, design.df=design.df, col=blkTerm[i-1]) 
   }   
        
   names(Z) <- c("e", blkTerm)  
   return(Z)
  }   
 
 #make the block projection matrix
  makeBlockProjectors <- function(BlkDesignMatrixList, initial = diag(nrow(BlkDesignMatrixList[[1]])) - mK(nrow(BlkDesignMatrixList[[1]]))){
 n <- nrow(BlkDesignMatrixList$e)
    Q <- lapply(BlkDesignMatrixList, function(z) projMat(z))
    
    Q <- Q[sort(1:length(Q), decreasing=TRUE)]
    
    cusumQ <- P <- NULL
    P[[1]] <- Q[[1]] %*% initial
    
    cusumQ[[1]] <- initial - P[[1]]    
    
    for(i in 2:length(Q)){
      P[[i]] <- Q[[i]] %*% cusumQ[[i-1]]
      cusumQ[[i]] <- cusumQ[[i-1]] - P[[i]]
       
    }   
    
    P <- P[sort(1:length(P), decreasing=TRUE)]
    names(P) <- names(BlkDesignMatrixList)           
    
    elementToRemove = numeric(0)
    for(i in 1:length(Q)){        
      if(all(P[[i]] <1e-6)){
        elementToRemove = c(elementToRemove, i)
      } 
    }
    
    if(length(elementToRemove)>0){
      P= P[-elementToRemove]    
    }                               
                             
    return(P)                   
  }
  
  #make the treatment C matrix
   makeTreatProjectors <- function(design.df, trtCols, effectsMatrix){       
   
    if(length(trtCols) == 1){
      nLevels = nlevels(design.df[,trtCols])
      names(nLevels) =  trtCols   
    }else if(any(grepl(":", trtCols))){
      uniqueTrtCols = unique(unlist(strsplit(trtCols, "\\:")))  
      nLevels <- sapply(design.df[,uniqueTrtCols], function(x) nlevels(as.factor(x)))
    }else{
      nLevels <- sapply(design.df[,trtCols], function(x) nlevels(as.factor(x)))
    }  
  
   nLevels = nLevels[rownames(effectsMatrix)]
   nEffects = ncol(effectsMatrix)  
   
   X <- as.list(rep(1,nEffects))
   names(X)  <- colnames(effectsMatrix)
   
   indMatrix <- function(x,n){
      if(x == 1) X <- mJ(n)
      else if(x == 2)  X <- mI(n)     
      else       X <- mK(n)
      return(X)
   } 
     
   for(i in 1:nrow(effectsMatrix)){
      
      matList <- lapply(effectsMatrix[i,], function(y) indMatrix(y, nLevels[i]))        
      for(j in 1:nEffects) X[[j]] <- X[[j]] %x% matList[[j]]    
   }
   
   return(X)  
  }
   
        
  #get the incidence matrix N
  getIncidenceMatrix <- function(design.df, trtCols){
   
   if(length(trtCols) == 1){
      incident = design.df[,trtCols]
      nLevels = levels(design.df[,trtCols])    
    }else if(any(grepl(":", trtCols))){ 
      uniqueTrtCols = unique(unlist(strsplit(trtCols, "\\:")))  

      incident = as.factor(apply(design.df[,uniqueTrtCols], 1, function(x) paste(x, collapse =".")))
      nLevels = sort(levels(interaction(design.df[,uniqueTrtCols])))

    }else{
     incident = as.factor(apply(design.df[,trtCols], 1, function(x) paste(x, collapse =".")))
     nLevels = sort(levels(interaction(design.df[,trtCols])))
    }       
   
   N <- matrix(0, nrow=nrow(design.df), ncol=length(nLevels))
   N[cbind(1:nrow(design.df), match(incident, nLevels))] <- 1
  
   return(N)
  }

  #mutiply mutiple rows or column of the matrix
  mutiplyRows = function(x){  
     unlist(lapply(strsplit(names(x), "\\."), function(y) cumprod(x[y])[length(y)]))  
  }

  #make the treatment ocntribution matrix 
  makeTreatContributMatrix <- function(design.df, trtCols){
    
    n <- length(trtCols)
    nRows <- nrow(design.df)
    X <- as.matrix(rep(1, nRows), nrow = nRows, ncol = 1)
    colnames(X) = "mu"
   
    for(i in 2:(n+1)){ 
      newX = makeDesignMatrix(nRows=nRows, design.df=design.df, col=trtCols[i-1]) 
      X <- cbind(X, newX - 1/ncol(newX))
    }   

    newX = t(apply(X, 1,  function(x)  unlist(lapply(strsplit(names(x), "\\."), function(y) cumprod(x[y])[length(y)]))))  
    colnames(newX) = colnames(X) 
    return(newX)
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

  #get eigenvalue using single value decomposition 
  eigenValue <- function(C,N,T){
   fractions(eigen(t(T) %*% t(N) %*% C %*% N %*% T)$va[1])
  }
  
  #pre- and post-multiply NTginvATN by block projection matrices  
  blkProkMat = 
    function(z, T, N){
    
      nEffect = length(T)
      PNTginvATNP = T
      PNTginvATNP[[1]] = z %*% N %*% T[[1]] %*% invInfMat(C=z, N=N, T=T[[1]]) %*% T[[1]] %*% t(N) %*% t(z) 

      newZ = (z %*% t(z)) - PNTginvATNP[[1]]
      
      if(nEffect !=1){
            
        for(i in 2:nEffect){
                      
          PNTginvATNP[[i]] = newZ %*% N %*% t(T[[i]]) %*% invInfMat(C=newZ, N=N, T=T[[i]]) %*% T[[i]] %*% t(N) %*% t(newZ) 
         
          newZ = (newZ %*% t(newZ)) - PNTginvATNP[[i]]
        }
      }

      PNTginvATNP$Error = newZ
        
      elementToRemove = numeric(0)
      for(i in 1:length(PNTginvATNP)){        
        if(all(PNTginvATNP[[i]] <1e-6)){
          elementToRemove = c(elementToRemove, i)
        } 
      }
    
      if(length(elementToRemove)>0){
        PNTginvATNP= PNTginvATNP[-elementToRemove]    
      }                               

      return(PNTginvATNP) 
    }

  #a function that computes the sum of the matrices for computing the PNTginvATNP
   matrixSum =
     function(matrixList){
      matrixSum = matrix(0, nrow = dim(matrixList[[1]])[1], ncol = dim(matrixList[[1]])[2])
      for(i in 1:length(matrixList)){      
        matrixSum = matrixSum + matrixList[[i]]
      }      
      return(matrixSum)   
   }
  
#########################################################################################  
#Main methods starts here->
#Extract the fixed and random terms   
  
  rT1 = terms(as.formula(paste("~", random.terms1, sep = "")), keep.order = TRUE) #random terms phase 1
  rT2 = terms(as.formula(paste("~", random.terms2, sep = "")), keep.order = TRUE) #random terms phase 2
  fT = terms(as.formula(paste("~", fixed.terms, sep = "")), keep.order = TRUE)  #fixed terms 
  
#########################################################################################  
#Preperating the block structures    
  print("Preperating the block structure.")

  blkTerm1 = attr(rT1,"term.labels")
  blkTerm2 = attr(rT2,"term.labels")
  
  Z1 = makeBlkDesMatrix(design.df, rev(blkTerm1))
  Z2 = makeBlkDesMatrix(design.df, rev(blkTerm2))
  
  print("Defining the block structures of second Phase.")
  Pb <- makeBlockProjectors(Z2)
  
  print("Defining the block structures of first phase within second Phase.")
  Pb1 <- lapply(Pb,function(z) makeBlockProjectors(Z1, z) )
  

#########################################################################################  
#Prepating the treatment structures
  print("Preperating the treatment structure.")

  trtTerm = attr(fT,"term.labels")
  effectsMatrix = attr(fT,"factor")

  T =  makeTreatProjectors(design.df, trtTerm, effectsMatrix)
  N =  getIncidenceMatrix(design.df, trtTerm)
  
  trtContr = makeTreatContributMatrix(design.df, trtTerm) 
#########################################################################################  
#Start calculating the VCs
#2-phase experiment
print("Start calculating the variance components.")
print("Pre- and post-multiply NTginvATN by block projection matrices.")  
   PNTginvATNP <- lapply(Pb1, function(y) lapply(y, function(z) blkProkMat(z, T, N)))
  
   #Now construct variance matrices
   V <- lapply(c(Z1, Z2[-1]), function(x) x %*% t(x))
  
   VC <- PNTginvATNP #to keep the structure
   contrT <- PNTginvATNP
   
   # And finally get coefficients of e, S and R VCs for each treatment effect in each stratum
print("Finally get coefficients for each treatment effect in each stratum.")

    for(i in 1:length(PNTginvATNP)){
          for(k in 1:length(PNTginvATNP[[i]])){ 
          tmp <- matrix(0, nrow=length(names(PNTginvATNP[[i]][[k]])), ncol=length(blkTerm1)+ length(blkTerm2)+ 1, dimnames=list(names(PNTginvATNP[[i]][[k]]), names(V)))
       
        for(j in 1:(length(names(PNTginvATNP[[i]][[k]])))){
          for(z in 1:(length(blkTerm1)+ length(blkTerm2)+ 1)){
            tmp[j,z] <- tr(PNTginvATNP[[i]][[k]][[j]] %*% V[[z]])                                   
          }         
          contrT[[i]][[k]][[j]] <- t(trtContr) %*% PNTginvATNP[[i]][[k]][[j]] %*% trtContr           
        }
      contrT[[i]][[k]] = lapply(contrT[[i]][[k]], function(x) if(all(x<1e-6)){ x = NULL} else{fractions(x)})              
      VC[[i]][[k]] <- fractions(tmp) 
      
      if(any(apply(tmp, 1, function(x) all(x < 1e-6)))){
       VC[[i]][[k]]<-  VC[[i]][[k]][-which(apply(tmp, 1, function(x) all(x < 1e-6))),,drop = FALSE]
      }
      }      
    }
  
  
   VC <- VC[sort(1:length(VC), decreasing=TRUE)]
   VC <- lapply(VC, function(x) x[sort(1:length(x), decreasing=TRUE)])
    
   contrT <- contrT[sort(1:length(contrT), decreasing=TRUE)]
   contrT <- lapply(contrT, function(x) x[sort(1:length(x), decreasing=TRUE)])
  

   return(list(random = VC, fixed = contrT))
}





#Our iTRAQ experiment
VC = getVCs.twoPhase.complex(design ,random.terms1 ="Cage/Animal/Sample/Subsample", random.terms2 = "Run",  fixed.terms = "Tag + Treat*Position")
VC$random
VC$fixed

VC = getVCs.twoPhase(design[-which(design$Run == "1"),] ,random.terms1 ="Cage/Animal/Sample/Subsample", random.terms2 = "Run",  fixed.terms = "Tag + Treat*Position")
VC$random
VC$fixed


VC = getVCs.twoPhase(design[-c(which(design$Run =="1"), which(design$Run =="8")),] ,random.terms1 ="Cage/Animal/Sample/Subsample", random.terms2 = "Run",  fixed.terms = "Tag + Treat*Position")
VC$random
VC$fixed

VC = getVCs.twoPhase(design[-c(which(design$Run =="1"), which(design$Run =="8")),] ,random.terms1 ="Cage/Animal/Sample/Subsample", random.terms2 = "Run",  fixed.terms = "Treat*Position + Tag")
VC$random
VC$fixed

VC1 = getVCs.twoPhase(design[-c(which(design$Run =="1"), which(design$Run =="7"), which(design$Run =="8")),] ,random.terms1 ="Cage/Animal/Sample/Subsample", random.terms2 = "Run",  fixed.terms = "Treat*Position + Tag")
VC$random
VC$fixed


