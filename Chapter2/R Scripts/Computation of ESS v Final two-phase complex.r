#2011/2/20 01:35:30
#Start from scratch
# input should be the design of data frame
#                 specify the random and fixed terms
#                  using the "terms" function to help to indicate cross and next relationships
#using the Nelder's method i.e. eignevalue decomposition to determine the fixed part of the EMS
# 
#The results should be in two: 1)ANOVA table structure with DF and coeffients of the variances compoenents of EMS for the random effects
#                              2)Treatment contribution matrices which may get complicated when the numbers of levels of treatment is large.
#25/02/2011 12:32:29 p.m.
# next thing to do:
#         a) add the coeffient of the treatment effects using the replication.



getVCs.twoPhase <- 
function(design.df, random.terms1, random.terms2, fixed.terms, var.comp = NA, trt.contr = NA, table.legend = FALSE){
  #design.df = data.frame contain the design of the experiment 
  #random.terms1 = a single string of characters contains the tier structrue for the random factors of the 1st phase using the W-R syntax. 
  #random.terms1 = a single string of characters contains the tier structrue for the random factors of the 2nd phase using the W-R syntax. 
  #fixed.terms = a single string of characters contains the tier structrue for the fixed factors using the W-R syntax. 
  #var.comps = a vector of characters contains the varaince componenets of interestm this allows the user to specify the variance componenets to be shown on the ANOVA table. Default is NA, which uses every random factors as the variance components with the 1st phase's variance components appear before the 2nd phase's variance components. 
  #trt.contr = a list of treatment contrast matrices, this allows the user to specify the contrasts for each treatment factor. Note that if the user uses this parameter, they need to specify the contrasts for every treatment factor. Default is NA, which uses the C matrix proposed by John and Williams.
  #table.legend = a logical, allow the users to using the lengend for the variance compoenents of the ANOVA table for a large design. Default is FALSE, which uses the orignal names.  

  
  library(MASS)   #need this package to perform ginv and factions functions
 
#######################################################################################  
  #All the pre-initiated functions   
  projMat <- function(X) X %*% ginv(t(X) %*% X) %*% t(X)  #projection matrix
  mK <- function(n) matrix(1/n, nrow=n, ncol=n)           #Averaging matrix
  mJ <- function(n) diag(n) - mK(n)                       #J matrix 
  tr <- function(X) sum(diag(X))                          #Trace operation
  m1 <- function(n) matrix(rep(1,n), nrow = n, ncol = 1)  #M1 matrix
  mI <- function(n) diag(n)                               #identity matrix
  
#######################################################################################  
  #check for factor names 
  isFactorNameNumeric <- function(levels) !as.logical( length(grep("[A-Z]|[a-z]", levels)) )

#######################################################################################  
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
      colNames <- paste(facName, 1:nCols,sep="")
    }else{                                    
      colNames <- levels(factor)
    }
    
    dimnames(Z) <- list(1:nRows, colNames)
    return(Z)
  } 
 
#######################################################################################  
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
 
#######################################################################################  
 #make the block projection matrix
  makeBlockProjectors <- function(BlkDesignMatrixList, initial = diag(nrow(BlkDesignMatrixList[[1]])) - mK(nrow(BlkDesignMatrixList[[1]]))){
   
    n <- nrow(BlkDesignMatrixList$e)
    Q <- lapply(BlkDesignMatrixList, function(z) projMat(z))
    
    Q <- Q[sort(1:length(Q), decreasing=TRUE)]
    
    elementToRemove = numeric(0)

    cusumQ <- P <- NULL
    P[[1]] <- Q[[1]] %*% initial     
                                                                                          
    if(all(P[[1]] <1e-6)){
     elementToRemove = 1                       #revord the elements that are all zeros
     P[[1]] <- matrix(0, nrow = n, ncol = n)   #make the projectors which has less than 1e-6 to zero, to avoid rounding error
    } 
        
    cusumQ[[1]] <- initial - P[[1]]    
    
    for(i in 2:(length(Q))){
      P[[i]] <- Q[[i]] %*% cusumQ[[i-1]]

      if((all(P[[i]] <1e-6))|| tr(P[[i]]) <= 0){
        elementToRemove = c(elementToRemove,i )   #revord the elements that are all zeros
        P[[i]] <- matrix(0, nrow = n, ncol = n)   #make the projectors which has less than 1e-6 to zero, to avoid rounding error
      }                 
      cusumQ[[i]] <- cusumQ[[i-1]] - P[[i]]
    }     
       
    P <- P[sort(1:length(P), decreasing=TRUE)]
    names(P) <- names(BlkDesignMatrixList)   
   
    if(length(elementToRemove)>0)
      P= P[-(length(Q) - elementToRemove+1)]
                                    
    return(P)                   
  }


#######################################################################################  
  #transfrom each treatment contrast metrix to C matrix
transContrToT = function(fact, contr){
  T = matrix(0, nrow = nlevels(fact), ncol= nlevels(fact))
  if(is.matrix(contr)){
    rownames(contr) = fact
    for(i in 1:ncol(contr)){         
      x = contr[levels(fact),i]
      T = T + (x %*% t(x))/as.numeric(t(x) %*% x)
    }
  
  } else{
     names(contr) = fact
     x = contr[levels(fact)]
      T = T + (x %*% t(x))/as.numeric(t(x) %*% x)      
  }                                                     
  return(T)
}                                                                        

  
#######################################################################################  
  #make the treatment C matrix
  makeTreatProjectors <- function(design.df, trtCols, effectsMatrix, trt.contr){       
   
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
   
   if(all(is.na(trt.contr))){
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
   }else{
      indMatrix <- function(x,n, trtContr){
        if(x == 1) X <- trtContr
        else if(x == 2)  X <- mI(n)     
        else       X <- mK(n)
        return(X)
     } 
       
     for(i in 1:nrow(effectsMatrix)){
        trtContr = transContrToT(design.df[,rownames(effectsMatrix)[i]], trt.contr[[rownames(effectsMatrix)[i]]])        
        matList <- lapply(effectsMatrix[i,], function(y) indMatrix(y, nLevels[i], trtContr))     
        for(j in 1:nEffects) X[[j]] <- X[[j]] %x% matList[[j]]    
     }   
   
   }
   
   return(X)  
  }
  
        
#######################################################################################  
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

#######################################################################################  
  #get a list of replication for each treatment
  getReplicationList <- function(design.df, trtCols){
    
    if(length(trtCols) == 1){
      return(mean(table(design.df[,trtCols])))
      
    }else if(any(grepl(":", trtCols))){ 
     
      trtColsList = lapply(strsplit(trtCols, "\\:"), function(x) design.df[,x])
      repList = sapply(trtColsList, function(y) if(is.factor(y)) {mean(table(y))} else{ mean(table(apply(y, 1, function(x) paste(x, collapse ="."))))} )
        
      levelList =  sapply(trtColsList, function(y) if(is.factor(y)) {nlevels(y)} else{ nlevels(as.factor(apply(y, 1, function(x) paste(x, collapse ="."))))} )
        
      repList = repList * levelList/nlevels(interaction(design.df[,unique(unlist(strsplit(trtCols, "\\:")))  ]))
     
      return(repList)      
    
    }else{
       repList = apply(design.df[,trtCols], 2, function(x)mean(table(x)))
       levelList = apply(design.df[,trtCols], 2, function(x) nlevels(as.factor(x)))
       repList = repList* levelList/ nlevels(interaction(design.df[,trtCols]))
     return(repList)      
        
    }                                                                 
  }    

#######################################################################################  
  #get a list of coeffient for each treatment
  getTrtCoef <- function(design.df, trtCols){
    
    if(length(trtCols) == 1){
      return(length(design.df[,trtCols])/nlevels(design.df[,trtCols]))
    }else if(any(grepl(":", trtCols))){ 
     
      trtColsList = lapply(strsplit(trtCols, "\\:"), function(x) design.df[,x])
      repList = sapply(trtColsList, function(y) if(is.factor(y)) {mean(table(y))} else{ mean(table(apply(y, 1, function(x) paste(x, collapse ="."))))} )
       
      return(repList)      
    }else{
       repList = apply(design.df[,trtCols], 2, function(x)mean(table(x)))
      return(repList)      
    }                                                                 
  }    

#######################################################################################  
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
#######################################################################################  
  #get eigenvalue using single value decomposition 
  eigenValue <- function(C,N,T){
   fractions(eigen(t(T) %*% t(N) %*% C %*% N %*% T)$va[1])
  }
  

#######################################################################################  
  #pre- and post-multiply NTginvATN by block projection matrices  
  blkProkMat = function(z, T, N, Rep){
      if(!is.matrix(z)) return(z)
      
      nEffect = length(T)
      PNTginvATNP = T
      
      PNTginvATNP[[1]] = z %*% N %*% T[[1]] %*% invInfMat(C=z, N=N, T=T[[1]]) %*% T[[1]] %*% t(N) %*% t(z) 

      va =Re( eigen(t(T[[1]]) %*% t(N) %*% z %*% N %*% T[[1]])$va)

      effFactor = mean(va[which(va>1e-6)]/Rep[names(T[1])])
      
      newZ = (z %*% t(z)) - PNTginvATNP[[1]]
      
      if(nEffect !=1){            
        for(i in 2:nEffect){                      
          PNTginvATNP[[i]] = newZ %*% N %*% t(T[[i]]) %*% invInfMat(C=newZ, N=N, T=T[[i]]) %*% T[[i]] %*% t(N) %*% t(newZ)    
          va =Re(eigen(t(T[[i]]) %*% t(N) %*% z %*% N %*% T[[i]])$va)
          effFactor = c(effFactor, mean(va[va>1e-6]/Rep[names(T[i])]))            
          newZ = (newZ %*% t(newZ)) - PNTginvATNP[[i]]
        }
      }

      PNTginvATNP$Residual = newZ
      elementToRemove = numeric(0)
      for(i in 1:length(PNTginvATNP)){        
        if(all(PNTginvATNP[[i]] <1e-6))
          elementToRemove = c(elementToRemove, i)
      }
    
      if(length(elementToRemove)>0)        
        PNTginvATNP= PNTginvATNP[-elementToRemove] 
                                     
      PNTginvATNP$effFactor = fractions(effFactor)
      
      names(PNTginvATNP$effFactor) = names(T)
      return(PNTginvATNP) 
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
  if(names(Pb)[1] == "e")
     names(Pb)[1] = "Within"
  
  print("Defining the block structures of first phase within second Phase.")
  trtTerm = attr(rT1,"term.labels")
  effectsMatrix = attr(rT1,"factor")

  T =  makeTreatProjectors(design.df, trtTerm, effectsMatrix, trt.contr = NA)
  N =  getIncidenceMatrix(design.df, trtTerm)
  Rep = getReplicationList(design.df, trtTerm)
  names(Rep) = names(T)
  
  
  Pb1 <- lapply(Pb,function(z) blkProkMat(z, T, N, Rep))

#########################################################################################  
#Prepating the treatment structures
  print("Preperating the treatment structure.")

  trtTerm = attr(fT,"term.labels")
  effectsMatrix = attr(fT,"factor")
                                                    
  T =  makeTreatProjectors(design.df, trtTerm, effectsMatrix, trt.contr)
  N =  getIncidenceMatrix(design.df, trtTerm)
  Rep = getReplicationList(design.df, trtTerm)
  trt.Coef = getTrtCoef(design.df, trtTerm)  
    
  names(Rep) = names(T)
  names(trt.Coef) = names(T)
  
 #########################################################################################  
#Start calculating the VCs
#2-phase experiment
print("Start calculating the variance components.")
print("Pre- and post-multiply NTginvATN by block projection matrices.")  
  
  
  PNTginvATNP <- lapply(Pb1, function(y) lapply(y, function(z) blkProkMat(z, T, N, Rep)))
   
  PNTginvATNP<- PNTginvATNP[sort(1:length(PNTginvATNP), decreasing=TRUE)]
      
  #Now construct variance matrices
  if(all(is.na(var.comp))){ 
    V <- lapply(c(Z1, Z2[-1]), function(x) x %*% t(x)) 
  }else{
    V <- lapply(makeBlkDesMatrix(design.df, var.comp), function(x) x %*% t(x))
  }
   
   VC <-character(length(V))
   names(VC)  = names(V)
   VC = t(data.frame(VC))
   
   trt  = numeric(length(trt.Coef) + length(Rep))
   names(trt) = c( names(T), paste("eff",  names(T), sep = ".") )
    
    
print("Finally get coefficients for each treatment effect in each stratum.")

    for(i in 1:length(PNTginvATNP)){
      VC = rbind(VC, character(length = length(V)))
      rownames(VC)[nrow(VC)] = paste("Between", names(PNTginvATNP[i]), sep = " ")
      
      for(k in 1:(length(PNTginvATNP[[i]])-1)){ 
        tmp <- matrix(0, nrow=length(names(PNTginvATNP[[i]][[k]])) -1, ncol=length(V), 
                  dimnames=list(names(PNTginvATNP[[i]][[k]])[- length(names(PNTginvATNP[[i]][[k]]))], names(V)))
       
        for(j in 1:((length(names(PNTginvATNP[[i]][[k]])))-1)){
          for(z in 1:length(V)){
              tmp[j,z] <- tr(PNTginvATNP[[i]][[k]][[j]] %*% V[[z]])                                   
          }                       
        }
      
 
      if(nrow(tmp) == 1 && rownames(tmp) == "Residual"){
        tmp = c(tmp[1], tmp[-1]/tmp[1])  
        VC = rbind(VC, attr(fractions(tmp),"fracs")) 
        if(names(PNTginvATNP[[i]][k]) == "Residual"  && length(PNTginvATNP[[i]]) > 2){
          rownames(VC)[nrow(VC)] =  paste("  ", names(PNTginvATNP[[i]][k]), sep = " ")
        }else if(names(PNTginvATNP[[i]][k]) == "Residual" && length(PNTginvATNP[[i]]) == 2){
          rownames(VC)[nrow(VC)] = paste("Between", names(PNTginvATNP[i]), sep = " ")
          VC = VC[-(nrow(VC)-1),]
        }else{
          rownames(VC)[nrow(VC)] = paste("   Between", names(PNTginvATNP[[i]][k]), sep = " ")
        }
      }else{
        VC = rbind(VC, character(length = length(V)))
        
        if(names(PNTginvATNP[[i]][k]) == "Residual"){
          rownames(VC)[nrow(VC)] =  paste("  ", names(PNTginvATNP[[i]][k]), sep = " ")
        }else{
          rownames(VC)[nrow(VC)] = paste("   Between", names(PNTginvATNP[[i]][k]), sep = " ")
        }         
        rownames(tmp) = paste("     ", rownames(tmp), sep = " ")
        
        tmp = t(apply(tmp, 1, function(x) attr(fractions(c(x[1], x[-1]/x[1])),"fracs")))
        VC = rbind(VC, tmp)
      }     
       trt.temp =  attr(fractions(c(trt.Coef*PNTginvATNP[[i]][[k]][[length(names(PNTginvATNP[[i]][[k]]))]],
                    PNTginvATNP[[i]][[k]][[length(names(PNTginvATNP[[i]][[k]]))]])),"fracs") 
      trt = rbind(trt,trt.temp)
      if( names(PNTginvATNP[i]) =="Within"){
         rownames(trt)[nrow(trt)] =  paste(names(PNTginvATNP[i]), names(PNTginvATNP[[i]][k]), sep = ".")
      }else{
        rownames(trt)[nrow(trt)] =  paste("Between", names(PNTginvATNP[i]), names(PNTginvATNP[[i]][k]), sep = ".")
      }
        
    }
   }
   
   
   trt = trt[-1,]

   trt = noquote(ifelse(trt == "NaN", "", trt))
   
   colnames(VC)[1] = "DF"
   VC = noquote(VC[-1,])   
    
   if(table.legend){
    Legend = paste(paste( letters[1:(length(colnames(VC))-1)],colnames(VC)[-1], sep = "="), collapse = ", ")
    colnames(VC)[-1] = letters[1:(length(colnames(VC))-1)]
    VC = list(VC = VC, legend = Legend)
   }
 
 
   return(list(random = VC, fixed = trt))
}

design.df = design

#Brien Payne's design
 design[,1] =as.factor(paste("T", design[,1], sep = ""))
 design[,2] =as.factor(paste("M", design[,2], sep = ""))
 
VC = getVCs.twoPhase.complex(design ,random.terms1 ="(Row*(Squ/Col))/Hal", random.terms2 = "((Oc/In/St)*Ju)/Pos",  fixed.terms = "Tre*Met")


VC$random
VC$effFactor
VC$fixed



                 && rownames(tmp) == "Residual"
