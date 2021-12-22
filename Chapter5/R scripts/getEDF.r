
#Transform the fraction to the numeric
fracToNum =
function(x){
  sapply(strsplit(x, "/"), function(s) ifelse(length(s) ==1, as.numeric(s),
            as.numeric(s)[1]/as.numeric(s)[2]))
  }

#Exact and organize the row names of the ANOVA table from inforDecompuTE
extractName= 
function(r.names){
   new.r.names = character(length(r.names))
  stratum.names = character(2)  
  for(i in 1:length(r.names)){
     if(grepl("^      ", r.names[i])){
      temp = c(stratum.names, gsub(" ", "",r.names[i]))
      new.r.names[i] = paste(temp, collapse = "")

      } else if(grepl("^   ", r.names[i])){
      temp = c(stratum.names[1], gsub(" ", "",r.names[i]))
      new.r.names[i] = paste(temp, collapse = "")
      stratum.names[2] = gsub(" ", "",r.names[i])
     } else {
       stratum.names[1] = r.names[i]
        new.r.names[i] = r.names[i]
     }

  }

   return(new.r.names)
}

#Pooling method
ptmv  =
function(VC.numeric){
  MS = VC.numeric[,"MS"]
  DF = VC.numeric[,"DF"]
  
  if(MS[4]>MS[3]){
    MS[4] = MS[3] =((DF[4] * MS[4]) + (DF[3] * MS[3]))/(DF[3] + DF[4])
  
    if(MS[4]>MS[2]){
      MS[4] = MS[3] = MS[2] =((DF[2] * MS[2]) + (DF[4] * MS[4]))/(DF[2] + DF[4])
      
      if(MS[2]>MS[1]){
        MS[4] = MS[3] = MS[2] = MS[1] =((DF[2] * MS[2]) + (DF[1] * MS[1]))/(DF[2] + DF[1])
  
    }
  
  }
 
  } else if(MS[4]>MS[2]){
    MS[4] = MS[2] =((DF[2] * MS[2]) + (DF[4] * MS[4]))/(DF[2] + DF[4])
    
      if(MS[2]>MS[1]){
        MS[4] = MS[2] = MS[1] =((DF[2] * MS[2]) + (DF[1] * MS[1]))/(DF[2] + DF[1])
  
    }

  }
  
  VC.numeric[,"MS"] = MS

  return(VC.numeric)
}

################################################################################
#Step 1:
#Extact the mean squares and degrees of freedom form the aov function
#These that means squares and degrees of freedoms are used to construct the initial
#information matrix
#Extact the variance components structure for inforDecompuTE to contrsut matrix G


################################################################################
#Step 3:
#REML to obtain the better estimates of the parameters of interest
#Use some random inital esitmates of variance componenets
#Get the score functions and then use the inverse of information matrix
#Start the fisher's scoring algorithm

#This function is to find the independant row set of linear combination of the
#variance componenets
getVcEDF = function(aov.table, MS = NA, row.MS = NA, true.VC = rep(NA, ncol(aov.table$A) - 2), 
    neg.VC = TRUE, ptmv = FALSE) {
    
    VC.numeric = t(apply(aov.table$A, 1, fracToNum))
    rownames(VC.numeric) = extractName(rownames(VC.numeric))
    #browser()
    if(!all(is.na(MS))){
      VC.numeric[,"MS"] = MS
    }         
    
    # Extract the mean squares for estimating variance components and EDF
    if (any(is.na(row.MS))) {
        EF.numeric = t(apply(aov.table$E, 1, fracToNum))
        
        rownames(EF.numeric) = extractName(rownames(EF.numeric))
        
        VC.numeric = VC.numeric[-which(apply(VC.numeric, 1, function(x) all(is.na(x)))), ]
        
        VC.numeric = VC.numeric[-match(names(which(!apply(EF.numeric, 1, function(x) all(is.na(x))))), 
            rownames(VC.numeric)), ]
    } else {
        VC.numeric = VC.numeric[row.MS, ]
    }
    
    # Contruct the Fisher's information matrix with respected to the mean squares
    m = VC.numeric[, "DF"]/(2 * VC.numeric[, "MS"]^2)
    m.mat = diag(m)
    
    # Extract G matrix, for transformation
    G.mat = VC.numeric[, -match(c("DF", "MS"), colnames(VC.numeric))]
    # G.mat[G.mat != 0 ] = 1

    # Fisher's information matrix
    infor.mat = t(G.mat) %*% m.mat %*% G.mat
    
    
    # Compute the variance componenets based on the linear combination
    LC.VC.numeric = VC.numeric[order(rowSums(G.mat)), ]
    LC.var.comp = numeric(ncol(G.mat))
    LC.G.mat = G.mat[order(rowSums(G.mat)), ]
    
    
    if (ncol(G.mat) == nrow(G.mat)) {
        LC.var.comp = ginv(LC.G.mat) %*% LC.VC.numeric[, "MS"]
    } else {
        # Use the qr decompositon to find the linear indenpendence
        qr.rank = qr(t(LC.G.mat))$rank
        qr.index = qr(t(LC.G.mat))$pivot[1:qr.rank]
        
        LC.var.comp = ginv(LC.G.mat[qr.index, ]) %*% LC.VC.numeric[qr.index, "MS"]
    }
    
    
    if (neg.VC && any(LC.var.comp < 0)) {
        LC.var.comp[which(LC.var.comp < 0)] = 0
    }   
  
    # Compute the variance componenets based on the Fisher's scoring of the REML method
    var.comp =  LC.var.comp # true.VC    #initial guess
    # var.comp[which(var.comp == 0)] = 1e-10
    
    old.var.comp = rep(0, ncol(G.mat))
    
    MSS = as.numeric(VC.numeric[, "MS"])
    DFF = as.numeric(VC.numeric[, "DF"])
    
    conv = matrix(0, nrow = 2, ncol = ncol(G.mat))
    conv[nrow(conv), ] = var.comp
    counter = 0
    # print(var.comp) Check over 100 iterations browser()
    while (!all(apply(conv, 1, function(x) isTRUE(all.equal(x, conv[nrow(conv), ], tolerance = 1e-7))))) {
        old.var.comp = var.comp
        old.MSS.hat = MSS.hat
        MSS.hat = G.mat %*% old.var.comp
        score.fun = t(G.mat) %*% (DFF * (MSS - MSS.hat)/(2 * MSS.hat^2))
        
        m.mat = diag(DFF/(2 * as.numeric(MSS.hat)^2))
        
        #n.infor.mat = t(sqrt(G.mat)) %*% m.mat %*% sqrt(G.mat)
        n.infor.mat = t(G.mat) %*% m.mat %*% G.mat
        inv.infor.mat = try(solve(n.infor.mat), TRUE)  
        #if the matrix is singluar
        if(class(inv.infor.mat) =="try-error"){
         inv.infor.mat = ginv(n.infor.mat) 
        }
         inv.infor.mat[which(abs(inv.infor.mat) <1e-10)]<-0
       
        var.comp = old.var.comp + inv.infor.mat %*% score.fun
         
        print(var.comp) 
        conv = rbind(conv, t(var.comp))
        conv = conv[-1, ]
        counter = counter + 1
        
        #force the variance compoenent estimates to be at zero, if it is below zero
        #This is done after 500 iterations, because the negative variance can sometimes 
        #fail the convergence during REML 
        var.comp[which(var.comp < 0)] = 0

       if (counter > 1000) {
            #if still fail to converge, thne last set of VC estimate is picked 
            break
        }
    }
    
    if (neg.VC && any(var.comp < 0)) {
        var.comp[which(var.comp < 0)] = 0
    }
    
    MSS.hat = G.mat %*% var.comp
    m.mat = diag(DFF/(2 * as.numeric(MSS.hat)^2))
    n.infor.mat = t(G.mat) %*% m.mat %*% G.mat
    inv.infor.mat = try(solve(n.infor.mat), TRUE)  
    #if the matrix is singluar
    if(class(inv.infor.mat) =="try-error"){
      inv.infor.mat = ginv(n.infor.mat) 
     }      
    inv.infor.mat[which(abs(inv.infor.mat) <1e-10)]<-0
    
    REML.EDF = numeric(length(m))
    
    rownames(var.comp) = colnames(G.mat)
    REML.var.comp = var.comp
    
    
    # Variance components from three methods
    var.comp = cbind(LC.var.comp, REML.var.comp, true.VC)
    
    colnames(var.comp) = c("LC.var.comp", "REML.var.comp", "TRUE.var.comp")
    
    lc.MSS.hat = G.mat %*% LC.var.comp
    lc.m.mat = diag(DFF/(2 * as.numeric(lc.MSS.hat)^2))
    lc.infor.mat = t(G.mat) %*% lc.m.mat %*% G.mat
    inv.lc.infor.mat = try(solve(lc.infor.mat), TRUE)  
    #if the matrix is singluar
    if(class(inv.lc.infor.mat) =="try-error"){
       inv.lc.infor.mat = ginv(lc.infor.mat) 
    }      
     inv.lc.infor.mat[which(abs(inv.lc.infor.mat) <1e-10)]<-0
    
    LC.EDF = numeric(length(m))
    
    if (!all(is.na(true.VC))) {
        true.MSS.hat = G.mat %*% true.VC
        true.m.mat = diag(DFF/(2 * as.numeric(true.MSS.hat)^2))
        true.infor.mat = t(G.mat) %*% true.m.mat %*% G.mat
        inv.true.infor.mat = solve(true.infor.mat)
         inv.true.infor.mat = try(solve(true.infor.mat), TRUE)  
    #if the matrix is singluar
    if(class(inv.true.infor.mat) =="try-error"){
        inv.true.infor.mat = ginv(true.infor.mat) 
      #inv.true.infor.mat = as.matrix(solve(nearPD(true.infor.mat)$mat)) 
    }      
       
        # inv.true.infor.mat[which(abs(inv.true.infor.mat) <1e-10)]<-0
        TRUE.EDF = numeric(length(m))
    }
    
    # Obtain the variances of all the paramters, than compute the EDF
    
    
    
    for (i in 1:nrow(G.mat)) {
        index = which(!G.mat[i, ] == 0)
        
        coef.mat = outer(G.mat[i, index], G.mat[i, index], "*")
        
        var.mat = inv.infor.mat[index, index]
        
        # Cmpute the EDF based on the VC from the REML
        REML.EDF[i] = (2 * MSS.hat[i]^2)/sum(coef.mat * var.mat)
        
        # Compute the EDF based on the VC from the linear combination
        lc.var.mat = inv.lc.infor.mat[index, index]
        
        LC.EDF[i] = (2 * lc.MSS.hat[i]^2)/sum(coef.mat * lc.var.mat)
 
        
        # Compute the EDF based on the true variance compoenents
        if (!all(is.na(true.VC))) {
            true.var.mat = inv.true.infor.mat[index, index]
            
            TRUE.EDF[i] = (2 * true.MSS.hat[i]^2)/sum(coef.mat * true.var.mat)
            
        }          
    }
    
    list(Stratum = cbind(VC.numeric, REML.EDF, LC.EDF, TRUE.EDF), Var.comp = var.comp)
}

 
##############o##############################################################################
getVcEDF.modified = function(aov.table, row.MS = NA, true.VC = rep(NA, ncol(aov.table$A) - 2), 
    neg.VC = TRUE, ptmv = FALSE) {
    
    VC.numeric = t(apply(aov.table$A, 1, fracToNum))
    rownames(VC.numeric) = extractName(rownames(VC.numeric))
    
    # Extract the mean squares for estimating variance components and EDF
    if (any(is.na(row.MS))) {
        EF.numeric = t(apply(aov.table$E, 1, fracToNum))
        
        rownames(EF.numeric) = extractName(rownames(EF.numeric))
        
        VC.numeric = VC.numeric[-which(apply(VC.numeric, 1, function(x) all(is.na(x)))), ]
        
        VC.numeric = VC.numeric[-match(names(which(!apply(EF.numeric, 1, function(x) all(is.na(x))))), 
            rownames(VC.numeric)), ]
    } else {
        VC.numeric = VC.numeric[row.MS, ]
    }
    
    # Pooling method
    if (ptmv) {
        VC.numeric = ptmv(VC.numeric)
    }
    
    # Contruct the Fisher's information matrix with respected to the mean squares
    m = VC.numeric[, "DF"]/(2 * VC.numeric[, "MS"]^2)
    m.mat = diag(m)
    
    # Extract G matrix, for transformation
    G.mat = VC.numeric[, -match(c("DF", "MS"), colnames(VC.numeric))]
    
    # Fisher's information matrix
    #infor.mat = t(G.mat) %*% m.mat %*% G.mat
    
    # infor.mat = ginv(sqrt(G.mat)) %*% m.mat %*% t(ginv(sqrt(G.mat)))
  
    # Compute the variance componenets based on the linear combination
    LC.VC.numeric = VC.numeric[order(rowSums(G.mat)), ]
    LC.var.comp = numeric(length(ncol(G.mat)))
    LC.G.mat = G.mat[order(rowSums(G.mat)), ]
    
    
    if (ncol(G.mat) == nrow(G.mat)) {
        LC.var.comp = ginv(LC.G.mat) %*% LC.VC.numeric[, "MS"]
    } else {
        # Use the qr decompositon to find the linear indenpendence
        qr.rank = qr(t(LC.G.mat))$rank
        qr.index = qr(t(LC.G.mat))$pivot[1:qr.rank]
        
        LC.var.comp = ginv(LC.G.mat[qr.index, ]) %*% LC.VC.numeric[qr.index, "MS"]
        
    }
    
    
    
    if (neg.VC && any(LC.var.comp < 0)) {
        LC.var.comp[which(LC.var.comp < 0)] = 0
    }
    
   
    # Compute the variance componenets based on the Fisher's scoring of the REML method
    var.comp = true.VC #LC.var.comp # true.VC    #initial guess
    # var.comp[which(var.comp == 0)] = 1e-10
    
    old.var.comp = rep(0, ncol(G.mat))
    
    MSS = as.numeric(VC.numeric[, "MS"])
    DFF = as.numeric(VC.numeric[, "DF"])
    
    conv = matrix(0, nrow = 2, ncol = ncol(G.mat))
    conv[nrow(conv), ] = var.comp
    counter = 0
    # print(var.comp) Check over 100 iterations browser()
    while (!all(apply(conv, 1, function(x) isTRUE(all.equal(x, conv[nrow(conv), ], tolerance = 1e-10))))) {
        old.var.comp = var.comp
        MSS.hat = G.mat %*% old.var.comp
        score.fun = ginv(G.mat) %*% (DFF * (MSS - MSS.hat)/(2 * MSS.hat^2))
  
        m.mat = diag(DFF/(2 * as.numeric(MSS.hat)^2))
        
        n.infor.mat = t(sqrt(G.mat)) %*% m.mat %*% sqrt(G.mat)
        #n.infor.mat = ginv(sqrt(G.mat)) %*% m.mat %*% t(ginv(sqrt(G.mat)))
        
        inv.infor.mat = ginv(n.infor.mat)
        
        var.comp = old.var.comp + inv.infor.mat %*% score.fun
        
        #print(var.comp) 
        conv = rbind(conv, t(var.comp))
        conv = conv[-1, ]
        counter = counter + 1
        if (counter > 1000) {
            stop("too many iterations")
        }
    }
    
    if (neg.VC && any(var.comp < 0)) {
        var.comp[which(var.comp < 0)] = 0
    }
    
    #browser()
    
    MSS.hat = G.mat %*% var.comp
    m.mat = diag(DFF/(2 * as.numeric(MSS.hat)^2))
    n.infor.mat = t(G.mat) %*% m.mat %*% G.mat

    
    # inv.infor.mat[which(abs(inv.infor.mat) <1e-10)]<-0
    
    REML.EDF = numeric(length(m))
    
    rownames(var.comp) = colnames(G.mat)
    REML.var.comp = var.comp
    
    
    # Variance components from three methods
    var.comp = cbind(LC.var.comp, REML.var.comp, true.VC)
    
    colnames(var.comp) = c("LC.var.comp", "REML.var.comp", "TRUE.var.comp")
    
    lc.MSS.hat = G.mat %*% LC.var.comp
    lc.m.mat = diag(DFF/(2 * as.numeric(lc.MSS.hat)^2))
    lc.infor.mat = t(G.mat) %*% lc.m.mat %*% G.mat
    
    
    # inv.lc.infor.mat[which(abs(inv.lc.infor.mat) <1e-10)]<-0
    
    LC.EDF = numeric(length(m))
    
    
    if (!all(is.na(true.VC))) {
        true.MSS.hat = G.mat %*% true.VC
        true.m.mat = diag(DFF/(2 * as.numeric(true.MSS.hat)^2))
        true.infor.mat = t(G.mat) %*% true.m.mat %*% G.mat
        
        # inv.true.infor.mat[which(abs(inv.true.infor.mat) <1e-10)]<-0
        TRUE.EDF = numeric(length(m))
    }
    
    # Obtain the variances of all the paramters, than compute the EDF
    
    
    
    for (i in 1:nrow(G.mat)) {
        index = which(!G.mat[i, ] == 0)
        
        coef.mat = outer(G.mat[i, index], G.mat[i, index], "*")
        
        #var.mat = inv.infor.mat[index, index]
        
        # Compute the EDF based on the VC from the REML
        #REML.EDF[i] = (2 * MSS.hat[i]^2)/sum(coef.mat * var.mat)
        
        REML.EDF[i] = (2 * MSS.hat[i]^2)/sum(coef.mat * solve(n.infor.mat[index,index]))
        
        # Compute the EDF based on the VC from the linear combination
        #lc.var.mat = inv.lc.infor.mat[index, index]
        
        #LC.EDF[i] = (2 * lc.MSS.hat[i]^2)/sum(coef.mat * lc.var.mat)
        LC.EDF[i] = (2 * lc.MSS.hat[i]^2)/sum(coef.mat * solve(lc.infor.mat[index,index]))
 
        
        # Compute the EDF based on the true variance compoenents
        if (!all(is.na(true.VC))) {
            #true.var.mat = inv.true.infor.mat[index, index]
            
            TRUE.EDF[i] =  (2 * true.MSS.hat[i]^2)/sum(coef.mat * solve(true.infor.mat[index,index]))

            
        }
        
    }
    
    list(Stratum = cbind(VC.numeric, REML.EDF, LC.EDF, TRUE.EDF), Var.comp = var.comp)
}
 




getFailVc=
function(aov.table, row.MS = NA, true.VC = rep(NA, ncol(aov.table$A) - 2), neg.VC = TRUE, ptmv = FALSE){

  VC.numeric = t(apply(aov.table$A, 1, fracToNum))
  rownames(VC.numeric) = extractName(rownames(VC.numeric))

  #Extract the mean squares for estimating variance components and EDF
  if(any(is.na(row.MS))){
    EF.numeric = t(apply(aov.table$E, 1, fracToNum))
       
    rownames(EF.numeric) = extractName(rownames(EF.numeric))
  
    VC.numeric = VC.numeric[-which(apply(VC.numeric, 1, function(x) all(is.na(x)))),]
  
    VC.numeric = VC.numeric[-match(names(which(!apply(EF.numeric, 1, function(x) all(is.na(x))))), 
                                        rownames(VC.numeric)),]
  }else{
    VC.numeric= VC.numeric[row.MS,]
  }
  
  #Pooling method
   if(ptmv){
   VC.numeric = ptmv(VC.numeric)
   }
 
  #Contruct the Fisher's information matrix with respected to the mean squares  
  m = VC.numeric[,"DF"]/(2*VC.numeric[,"MS"]^2) 
  m.mat = diag(m)   
  
  #Extract G matrix, for transformation 
  G.mat = VC.numeric[, -match(c("DF", "MS"), colnames(VC.numeric))]

  #Fisher's information matrix
  infor.mat = t(G.mat) %*% m.mat %*% G.mat


  #Compute the variance componenets based on the linear combination
    LC.VC.numeric = VC.numeric[order(rowSums(G.mat)),]
    LC.var.comp = numeric(length(ncol(G.mat)))
    LC.G.mat = G.mat[order(rowSums(G.mat)),]
     
    
    if(ncol(G.mat) == nrow(G.mat)){     
      LC.var.comp =  solve(LC.G.mat) %*%LC.VC.numeric[,"MS"]       
    } else{
      #Use the qr decompositon to find the linear idenpendence    
      qr.rank = qr(t(LC.G.mat))$rank
      qr.index = qr(t(LC.G.mat))$pivot[1:qr.rank]
      
      LC.var.comp =  solve(LC.G.mat[qr.index,])  %*% LC.VC.numeric[qr.index ,"MS"]  
     }

          
      if(neg.VC){
        LC.var.comp[which(LC.var.comp<0)] = 0 
      }


  #Compute the variance componenets based on the Fisher's scoring of the REML method 
  var.comp = true.VC # LC.var.comp    #initial guess
  old.var.comp = rep(0, ncol(G.mat))          
   
  MSS = as.numeric(VC.numeric[,"MS"])
  DFF = as.numeric(VC.numeric[,"DF"])
    
  conv = matrix(0, nrow = 2, ncol = ncol(G.mat))
  conv[nrow(conv), ] = var.comp
  counter = 0
  #print(var.comp)
 #Check over 100 iterations 
  while(! all(apply(conv, 1, function(x) isTRUE(all.equal(x, conv[nrow(conv),], tolerance = 1e-7))))){
    old.var.comp = var.comp
    MSS.hat = G.mat %*%  old.var.comp 
    score.fun = t(G.mat) %*% (DFF * (MSS - MSS.hat) /(2 * MSS.hat^2)) 
    m.mat = diag(DFF/(2*as.numeric(MSS.hat)^2))
    n.infor.mat =  t(G.mat) %*% m.mat %*% G.mat
    var.comp  = old.var.comp  +  solve(n.infor.mat) %*% score.fun
    #print(var.comp)
   
    #print(var.comp)
    conv = rbind(conv, t(var.comp))
    conv = conv[-1,]
    counter = counter + 1
    if(counter>100){
      var.comp = t(conv)
      break
    }
  }

    if(neg.VC){
      var.comp[which(var.comp<0)] = 0 
    }

  
   rownames(var.comp) = colnames(G.mat)
   REML.var.comp = var.comp
   
      
    #Variance components from three methods
    var.comp =  cbind(LC.var.comp, 
                      REML.var.comp,
                      true.VC) 
    
    colnames(var.comp) = c("LC.var.comp","REML.var.comp1","REML.var.comp2","TRUE.var.comp")
    
  
  list(Var.comp = var.comp)
}
