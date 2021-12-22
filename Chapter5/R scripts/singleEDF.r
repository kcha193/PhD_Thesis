


# sourceDir <- function(path, trace = TRUE, ...) {
#   library(MASS)
#   library(inline)
#   library(compiler) 
#   library(formatR)
#   library(Rcpp)
#   library(RcppArmadillo)
#   library(ggplot2)
#   library(gridExtra)
#   library(reshape)
#   library(snowfall)
#   for (nm in list.files(path, pattern = "\\.[RrSsQq]$")) {
#     if(trace) cat(nm,":")
#       source(file.path(path, nm), ...)
#       if(trace) cat("\n")
#   }
# }
# 
# sourceDir(path = "C:/Users/kcha193/My Dropbox/Packaging tools/packaging/infoDecompuTE/R")
# sourceDir(path = "C:/Users/kcha193/My Dropbox/Packaging tools/packaging/optimTE/R")

library(optimTE)

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



getVcEDF = 
function(aov.table, MS = NA, row.MS = NA, true.VC = rep(NA, ncol(aov.table$A) - 2), 
neg.VC = TRUE) {
  
  
    VC.numeric = t(apply(aov.table$A, 1, fracToNum))
    rownames(VC.numeric) = extractName(rownames(VC.numeric))
    
    if (!all(is.na(MS))) {
        VC.numeric[, "MS"] = MS
    }
    
   # Extract the mean squares for estimating variance components and EDF
    if (any(is.na(row.MS))) {
        EF.numeric = t(apply(aov.table$Fixed$EF, 1, fracToNum))
        
        rownames(EF.numeric) = extractName(rownames(EF.numeric))
        
        VC.numeric = VC.numeric[-which(apply(VC.numeric, 1, function(x) all(is.na(x)))), ]
        
        VC.numeric = VC.numeric[-match(names(which(!apply(EF.numeric, 1, function(x) all(is.na(x))))), 
            rownames(VC.numeric)), ]
    } else {
        VC.numeric = VC.numeric[row.MS, ]
    }
     
     if(any(apply(VC.numeric, 2, function(x) all(x==0))))
      VC.numeric = VC.numeric[,-which(apply(VC.numeric, 2, function(x) all(x==0)))]

    
    
    # Extract G matrix, for transformation
    G.mat = VC.numeric[, -match(c("DF", "MS"), colnames(VC.numeric))]
    # G.mat[G.mat != 0 ] = 1
    
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
    
    LC.var.comp
    
    if (neg.VC && any(LC.var.comp < 0)) {
        LC.var.comp[which(LC.var.comp < 0)] = 0
    }
    
    # Compute the variance componenets based on the Fisher's scoring of the REML method
    
    var.comp =  as.numeric(if( all(is.na(true.VC))){ LC.var.comp} else{true.VC}) #runif(ncol(G.mat), 0, 10) #initial guess
    
    old.var.comp = rep(0, ncol(G.mat))
    neg.var.comp = rep(0, ncol(G.mat))
    
    MSS = as.numeric(VC.numeric[, "MS"])
    DFF = as.numeric(VC.numeric[, "DF"])
    
    conv = matrix(0, nrow = 2, ncol = ncol(G.mat))
    conv[nrow(conv), ] = var.comp
    counter = 0
    # print(var.comp) Check over 100 iterations browser()
    while (!all(apply(conv, 1, function(x) isTRUE(all.equal(x, conv[nrow(conv), ], tolerance = 1e-07))))) {
        old.var.comp = var.comp
        
        MSS.hat = as.numeric(G.mat %*% old.var.comp)
        score.fun = t(G.mat) %*% (DFF * (MSS - MSS.hat)/(2 * MSS.hat^2))
        m.mat = diag(DFF/(2 * MSS.hat^2))
        n.infor.mat = t(G.mat) %*% m.mat %*% G.mat
        var.comp = try(old.var.comp + solve(n.infor.mat) %*% score.fun, TRUE)
        
        if (class(var.comp) == "try-error") {
          #stop()
          var.comp = rep(NA, length( old.var.comp))
            break
          
        }
        
        if (neg.VC && any(var.comp < 0)) {
          var.comp[which(var.comp < 0)] =  1e-7
        }
        #print(var.comp)
        # force the variance compoenent estimates to be at zero, if it is below zero This is done after 500
        # iterations, because the negative variance can sometimes fail the convergence during REML
        if ( counter > 1000) {
          var.comp = rep(NA, length( old.var.comp))
                   break
        }
        
        conv = rbind(conv, t(var.comp))
        conv = conv[-1, ]
        counter = counter + 1
      
         }
    
    #if(counter > 1000){
    #  var.comp = apply(conv, 2, mean)
    #}
    
    if ((!all(is.na(var.comp))) && neg.VC && any(var.comp == 1e-7)) {
        var.comp[which(var.comp == 1e-7)] = 0
    }
    
    
    names(var.comp) = colnames(G.mat)
    REML.var.comp = var.comp
    
    if(!all(is.na(REML.var.comp))){
      MSS.hat = as.numeric(G.mat %*% REML.var.comp)
      m.mat = diag(DFF/(2 * (MSS.hat)^2))
      n.infor.mat = t(G.mat) %*% m.mat %*% G.mat
      inv.infor.mat = try(solve(n.infor.mat), TRUE)
      # if the matrix is singluar
    
      if (class(inv.infor.mat) == "try-error") {
      stop()
          REML.var.comp[which(REML.var.comp < 0)] = 0
          MSS.hat = as.numeric(G.mat %*% REML.var.comp)
          m.mat = diag(DFF/2 * (MSS.hat)^2)
          infor.mat = t(G.mat) %*% m.mat %*% G.mat
          inv.infor.mat = solve(n.infor.mat)
      }
    }
    
    lc.MSS.hat = as.numeric(G.mat %*% LC.var.comp)
    lc.m.mat = diag(DFF/(2 * (lc.MSS.hat)^2))
    lc.infor.mat = t(G.mat) %*% lc.m.mat %*% G.mat
    inv.lc.infor.mat = try(solve(lc.infor.mat), TRUE)
    if (class( inv.lc.infor.mat) == "try-error") {
          stop()
      }    

    true.MSS.hat = as.numeric(G.mat %*% true.VC)
    true.m.mat = diag(DFF/(2 * (true.MSS.hat)^2))
    true.infor.mat = t(G.mat) %*% true.m.mat %*% G.mat
    inv.true.infor.mat = solve(true.infor.mat)
    
    
    # Variance components from three methods
    var.comp = cbind(REML.var.comp, LC.var.comp, true.VC)
    
    colnames(var.comp) = c("REML.var.comp", "LC.var.comp",  "TRUE.var.comp")
    var.comp                                      
    
    #Obtain the variances of all the paramters, than compute the EDF
    
    REML.EDF = numeric(nrow(G.mat))
    LC.EDF = numeric(nrow(G.mat))
    TRUE.EDF = numeric(nrow(G.mat))
    
    for (i in 1:nrow(G.mat)) {
        index = which(!G.mat[i, ] == 0)
        
        coef.mat = outer(G.mat[i, index], G.mat[i, index], "*")

        # Cmpute the EDF based on the VC from the REML
        if(all(is.na(REML.var.comp))){
          REML.EDF[i] = NA
        } else {
          var.mat = inv.infor.mat[index, index]

          REML.EDF[i] = (2 * MSS.hat[i]^2)/sum(coef.mat * var.mat)
        }
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
 
  
getVcEDFNew = 
  function(aov.table, MS = NA, row.MS = NA, true.VC = rep(NA, ncol(aov.table$A) - 2), 
           neg.VC = TRUE) {
    
    
    VC.numeric = t(apply(aov.table$A, 1, fracToNum))
    rownames(VC.numeric) = extractName(rownames(VC.numeric))
    
    if (!all(is.na(MS))) {
      VC.numeric[, "MS"] = MS
    }
    
    # Extract the mean squares for estimating variance components and EDF
    if (any(is.na(row.MS))) {
      EF.numeric = t(apply(aov.table$Fixed$EF, 1, fracToNum))
      
      rownames(EF.numeric) = extractName(rownames(EF.numeric))
      
      VC.numeric = VC.numeric[-which(apply(VC.numeric, 1, function(x) all(is.na(x)))), ]
      
      VC.numeric = VC.numeric[-match(names(which(!apply(EF.numeric, 1, function(x) all(is.na(x))))), 
                                     rownames(VC.numeric)), ]
    } else {
      VC.numeric = VC.numeric[row.MS, ]
    }
    
    if(any(apply(VC.numeric, 2, function(x) all(x==0))))
      VC.numeric = VC.numeric[,-which(apply(VC.numeric, 2, function(x) all(x==0)))]
    
    
    
    # Extract G matrix, for transformation
    G.mat = VC.numeric[, -match(c("DF", "MS"), colnames(VC.numeric))]
    # G.mat[G.mat != 0 ] = 1
    
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
    
    LC.var.comp
    
    if (neg.VC && any(LC.var.comp < 0)) {
      LC.var.comp[which(LC.var.comp < 0)] = 0
    }
    
    # Compute the variance componenets based on the Fisher's scoring of the REML method
    
    var.comp =  as.numeric(if( all(is.na(true.VC))){ LC.var.comp} else{true.VC}) #runif(ncol(G.mat), 0, 10) #initial guess
    
    old.var.comp = rep(0, ncol(G.mat))
    neg.var.comp = rep(0, ncol(G.mat))
    
    MSS = as.numeric(VC.numeric[, "MS"])
    DFF = as.numeric(VC.numeric[, "DF"])
    
    conv = matrix(0, nrow = 2, ncol = ncol(G.mat))
    conv[nrow(conv), ] = var.comp
    counter = 0
    # print(var.comp) Check over 100 iterations browser()
    while (!all(apply(conv, 1, function(x) isTRUE(all.equal(x, conv[nrow(conv), ], tolerance = 1e-07))))) {
      old.var.comp = var.comp
      
      MSS.hat = as.numeric(G.mat %*% old.var.comp)
      score.fun = t(G.mat) %*% (DFF * (MSS - MSS.hat)/(2 * MSS.hat^2))
      m.mat = diag(DFF/(2 * MSS.hat^2))
      n.infor.mat = t(G.mat) %*% m.mat %*% G.mat
      var.comp = try(old.var.comp + solve(n.infor.mat) %*% score.fun, TRUE)
      
      if (class(var.comp) == "try-error") {
        #stop()
        var.comp = rep(NA, length( old.var.comp))
        break
        
      }
      
      if (neg.VC && any(var.comp < 0)) {
        var.comp[which(var.comp < 0)] =  1e-7
      }
      #print(var.comp)
      # force the variance compoenent estimates to be at zero, if it is below zero This is done after 500
      # iterations, because the negative variance can sometimes fail the convergence during REML
      if ( counter > 1000) {
        var.comp = rep(NA, length( old.var.comp))
        break
      }
      
      conv = rbind(conv, t(var.comp))
      conv = conv[-1, ]
      counter = counter + 1
      
    }
    
    #if(counter > 1000){
    #  var.comp = apply(conv, 2, mean)
    #}
    
    if ((!all(is.na(var.comp))) && neg.VC && any(var.comp == 1e-7)) {
      var.comp[which(var.comp == 1e-7)] = 0
    }
    
    
    names(var.comp) = colnames(G.mat)
    REML.var.comp = var.comp
    
    if(!all(is.na(REML.var.comp))){
      MSS.hat = as.numeric(G.mat %*% REML.var.comp)
      m.mat = diag(DFF/(2 * (MSS.hat)^2))
      n.infor.mat = t(G.mat) %*% m.mat %*% G.mat
      inv.infor.mat = try(solve(n.infor.mat), TRUE)
      # if the matrix is singluar
      
      if (class(inv.infor.mat) == "try-error") {
        stop()
        REML.var.comp[which(REML.var.comp < 0)] = 0
        MSS.hat = as.numeric(G.mat %*% REML.var.comp)
        m.mat = diag(DFF/2 * (MSS.hat)^2)
        infor.mat = t(G.mat) %*% m.mat %*% G.mat
        inv.infor.mat = solve(n.infor.mat)
      }
    }
    
    lc.MSS.hat = as.numeric(G.mat %*% LC.var.comp)
    lc.m.mat = diag(DFF/(2 * (lc.MSS.hat)^2))
    lc.infor.mat = t(G.mat) %*% lc.m.mat %*% G.mat
    inv.lc.infor.mat = try(solve(lc.infor.mat), TRUE)
    if (class( inv.lc.infor.mat) == "try-error") {
      stop()
    }    
    
    true.MSS.hat = as.numeric(G.mat %*% true.VC)
    true.m.mat = diag(DFF/(2 * (true.MSS.hat)^2))
    true.infor.mat = t(G.mat) %*% true.m.mat %*% G.mat
    inv.true.infor.mat = solve(true.infor.mat)
    
    
    # Variance components from three methods
    var.comp = cbind(REML.var.comp, LC.var.comp, true.VC)
    
    colnames(var.comp) = c("REML.var.comp", "LC.var.comp",  "TRUE.var.comp")
    var.comp                                      
    
    #Obtain the variances of all the paramters, than compute the EDF
    
    REML.EDF = numeric(nrow(G.mat))
    LC.EDF = numeric(nrow(G.mat))
    TRUE.EDF = numeric(nrow(G.mat))

    for (i in 1:nrow(G.mat)) {
      index = which(!G.mat[i, ] == 0)
      
      coef.mat = outer(G.mat[i, index], G.mat[i, index], "*")
      
      # Cmpute the EDF based on the VC from the REML
      if(all(is.na(REML.var.comp))){
        REML.EDF[i] = NA
      } else {
        var.mat = inv.infor.mat[index, index]
        
        REML.EDF[i] = (2 * MSS.hat[i]^2)/sum(coef.mat * var.mat)
      }
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



edfPlot =
  function(sim,
           rowNames = "Within RunBetweenAniResidual",
           pTitle = "EDF plots") {
    
    tempS.bwAni = sim$tempS[grep(rowNames, rownames(sim$tempS)), ]
    
    mm = melt(tempS.bwAni)
    
    names(mm)[3] = "Method"
    
    mm$gamma.ani =   as.numeric(as.character(mm$gamma.ani))
    mm$gamma.run =   ifelse(
      is.na(mm$gamma.run),
      "plain(\"Runs Fixed\")",
      paste("sigma[r]^2/sigma^2 ==", mm$gamma.run)
    )
   
    mm$gamma.run = factor(mm$gamma.run, levels =  unique(mm$gamma.run))
    
    pg <- ggplot(mm, aes(x = log(gamma.ani), y = value)) + ylab("EDF") +
      xlab(expression(paste(sigma[a] ^ 2, " / ", sigma ^ 2))) +
      ggtitle(pTitle) +
      scale_x_continuous(breaks = log(unique(mm$gamma.ani))[seq(1, 17, 2)],
                         labels = as.character(unique(round(mm$gamma.ani, 4)))[seq(1, 17, 2)]) +
      geom_line(aes(colour = Method), size = 1, alpha = I(.5)) +
      facet_wrap(~  gamma.run, labeller = label_parsed) +
      scale_colour_discrete(name = "Method") + theme(legend.position = "top") +
      theme_bw()
    
    pg
  }

edfPlot.old =
  function(sim,
           rowNames = "Within RunBetweenAniResidual",
           pTitle = "EDF plots") {
    
    tempS.bwAni = sim$tempS[grep(rowNames, rownames(sim$tempS)), ]
    
    mm = melt(tempS.bwAni)
    
    names(mm)[3] = "Method"
    
    mm$gamma.ani =   as.numeric(as.character(mm$gamma.ani))
    #mm$gamma.run =  as.factor( paste("gamma.run =", mm$gamma.run))
    #levels(mm$gamma.run ) = sort(levels(mm$gamma.run ))
    

    if (any(is.na(mm$gamma.run))) {
      level.temp =  levels(mm$gamma.run)
      temp =  as.character(mm$gamma.run)
      temp[which(is.na(temp))] = "Fixed"
      mm$gamma.run = factor(temp,  levels =  c(level.temp, "Fixed"))
    }
    
    pg <-
      ggplot(mm, aes(x = log(gamma.ani), y = value), group = gamma.run) + ylab("EDF") +
      xlab(expression(paste(sigma[A] ^ 2, " / ", sigma ^ 2))) +
      ggtitle(pTitle) +
      scale_x_continuous(breaks = log(unique(mm$gamma.ani))[seq(1, 17, 2)],
                         labels = as.character(unique(round(mm$gamma.ani, 4)))[seq(1, 17, 2)]) +
      geom_line(aes(colour = Method), size = 1, alpha = I(.5)) +
      facet_wrap( ~  gamma.run) +
      scale_colour_discrete(name = "Method")    + theme(legend.position = "top")
    
    pg
  }

edfPlotCompare =
  function(sim,
           sim1,
           compare,
           rowNames = "Within RunBetweenAniResidual",
           pTitle = "EDF plots") {
    tempS.bwAni = sim$tempS[grep(rowNames, rownames(sim$tempS)), c(1, 2)]
    tempS.bwAni = cbind(tempS.bwAni, sim$tempS[grep(rowNames, rownames(sim$tempS)), 3])
    tempS.bwAni = cbind(tempS.bwAni, sim1$tempS[grep(rowNames, rownames(sim1$tempS)), 3])
    
    colnames(tempS.bwAni) [c(3, 4)] =  compare
    mm = melt(tempS.bwAni)
    
    colnames(mm) [3] = "Design"
    
    mm$gamma.ani =   as.numeric(as.character(mm$gamma.ani))
    mm$gamma.run =   ifelse(
      is.na(mm$gamma.run),
      "plain(\"Runs Fixed\")",
      paste("sigma[r]^2/sigma^2 ==", mm$gamma.run)
    )
    #levels(mm$gamma.run ) = sort(levels(mm$gamma.run ))
    
    mm$gamma.run = factor(mm$gamma.run, levels =  unique(mm$gamma.run))
    
    pg <- ggplot(mm, aes(x = log(gamma.ani), y = value)) + ylab("EDF") +
      xlab(expression(paste(sigma[a] ^ 2, " / ", sigma ^ 2))) +
      ggtitle(pTitle) +
      scale_x_continuous(breaks = log(unique(mm$gamma.ani))[seq(1, 17, 2)],
                         labels = as.character(unique(round(mm$gamma.ani, 4)))[seq(1, 17, 2)]) +
      geom_line(aes(colour = Design), size = 1, alpha = I(.5)) +
      facet_wrap(~  gamma.run, labeller = label_parsed) +
      #facet_wrap( ~  gamma.run, labeller = label_bquote(sigma[r]^2/sigma^2 == .(as.character(gamma.run)))) +
      scale_colour_discrete(name = "Design") + theme(legend.position = "top") +
      theme_bw()
    pg
  }



edfPlotCompare.REML =
  function(sim, sim1, compare, rowNames = "Within RunBetweenAniResidual", pTitle = "EDF plots"){
    
    tempS.bwAni = sim$tempS[grep(rowNames, rownames(sim$tempS)), c(1,2)]
    tempS.bwAni = cbind(tempS.bwAni, sim$tempS[grep(rowNames, rownames(sim$tempS)), c(3,4)])
    tempS.bwAni = cbind(tempS.bwAni, sim1$tempS[grep(rowNames, rownames(sim1$tempS)), c(3,4)])
    
    colnames(tempS.bwAni) [3:6] = paste0(rep(compare, each = 2), "/", colnames(tempS.bwAni) [3:6])
    mm = melt(tempS.bwAni)  
    
    colnames(mm) [3] = "Design"
    
   mm$gamma.ani =   as.numeric(as.character(mm$gamma.ani))
    mm$gamma.run =   ifelse(
      is.na(mm$gamma.run),
      "plain(\"Runs Fixed\")",
      paste("sigma[r]^2/sigma^2 ==", mm$gamma.run)
    )
    
    mm$gamma.run = factor(mm$gamma.run, levels =  unique(mm$gamma.run))

    
    pg <- ggplot(mm, aes(x = log(gamma.ani), y=value))+ ylab("EDF") +
      xlab(expression(paste(sigma[a]^2," / ",sigma^2))) +
      ggtitle(pTitle) +
      scale_x_continuous( breaks = log(unique(mm$gamma.ani))[seq(1, 17, 2)],
                          labels = as.character(unique(round(mm$gamma.ani, 4)))[seq(1, 17, 2)]) +
      geom_line(aes(colour = Design), size = 1, alpha=I(.5)) +
      facet_wrap(~  gamma.run, labeller = label_parsed) +
       scale_colour_discrete(name = "Design/Method")    + theme(legend.position = "top") + 
      theme_bw()     
    
    pg
  }



edfPlotCompareCag.REML =
  function(sim, sim1, compare, rowNames = 
             "Within RunBetweenCag(Ani)Residual", pTitle = "EDF plots"){
    
    tempS.bwAni = sim[, c(1,2,3)]
    tempS.bwAni = cbind(tempS.bwAni, sim[, c(4,5)])
    tempS.bwAni = cbind(tempS.bwAni, sim1[, c(4,5)])
    
    colnames(tempS.bwAni) [4:7] = paste0(rep(compare, each = 2), "/", colnames(tempS.bwAni) [4:7])
    mm = melt(tempS.bwAni)  
    
    colnames(mm) [4] = "Design"
    
   mm$gamma.ani =   as.numeric(as.character(mm$gamma.ani))
    mm$gamma.run =   ifelse(
      is.na(mm$gamma.run),
      "plain(\"Runs Fixed\")",
      paste("sigma[r]^2/sigma^2 ==", mm$gamma.run)
    )
    
    mm$gamma.run = factor(mm$gamma.run, levels =  c( "sigma[r]^2/sigma^2 == 0", "sigma[r]^2/sigma^2 == 0.25", "sigma[r]^2/sigma^2 == 1", 
       "sigma[r]^2/sigma^2 == 4", "sigma[r]^2/sigma^2 == 100" , "plain(\"Runs Fixed\")"))

      mm$gamma.cag =   ifelse(
      is.na(mm$gamma.cag),
      "plain(\"Trays Fixed\")",
      paste("sigma[b]^2/sigma^2 ==", mm$gamma.cag)
    )
    
    mm$gamma.cag = factor(mm$gamma.cag, levels =  c( "sigma[b]^2/sigma^2 == 0", 
                                                     "sigma[b]^2/sigma^2 == 0.25", 
                                                     "sigma[b]^2/sigma^2 == 1", 
       "sigma[b]^2/sigma^2 == 4", "sigma[b]^2/sigma^2 == 100" , "plain(\"Trays Fixed\")"))
    
    pg <- ggplot(mm, aes(x = log(gamma.ani), y=value))+ ylab("EDF") +
      xlab(expression(paste(sigma[p]^2," / ",sigma^2))) +
      ggtitle(pTitle) +
      scale_x_continuous( breaks = log(unique(mm$gamma.ani))[seq(1, 17, 2)],
                          labels = as.character(unique(round(mm$gamma.ani, 4)))[seq(1, 17, 2)]) +
      geom_line(aes(colour = Design), size = 1, alpha=I(.5)) +
      facet_wrap(c("gamma.run", "gamma.cag"), labeller = label_parsed) +
       scale_colour_discrete(name = "Design/Method")    + theme(legend.position = "top") + 
      theme_bw()  
    
      pg
  }




edfPlotCompare.LC =
  function(sim, sim1, compare, rowNames = "Within RunBetweenAniResidual", pTitle = "EDF plots"){
    
    tempS.bwAni = sim$tempS[grep(rowNames, rownames(sim$tempS)), c(1,2)]
    tempS.bwAni = cbind(tempS.bwAni, apply(sim$tempS[grep(rowNames, rownames(sim$tempS)), 4], 1, max))
    tempS.bwAni = cbind(tempS.bwAni, apply(sim1$tempS[grep(rowNames, rownames(sim1$tempS)), 4], 1, max))
    
    colnames(tempS.bwAni) [c(3,4)]=  compare
    mm = melt(tempS.bwAni)  
    
    colnames(mm) [3] = "Design"
    
    mm$gamma.ani =   as.numeric(as.character(mm$gamma.ani))
    #mm$gamma.run =  as.factor( paste("gamma.run =", mm$gamma.run))
    #levels(mm$gamma.run ) = sort(levels(mm$gamma.run ))
    if(any(is.na(mm$gamma.run))){
      level.temp =  levels(mm$gamma.run)
      temp =  as.character(mm$gamma.run)
      temp[which(is.na(temp))] = "Fixed"
      mm$gamma.run = factor(temp,  levels =  c(level.temp, "Fixed"))
    }
    
    
    pg <- ggplot(mm, aes(x = log(gamma.ani), y=value))+ ylab("EDF") +
      xlab(expression(paste(sigma[a]^2," / ",sigma^2))) +
      ggtitle(pTitle) +
      scale_x_continuous( breaks = log(unique(mm$gamma.ani))[seq(1, 17, 2)],
                          labels = as.character(unique(round(mm$gamma.ani, 4)))[seq(1, 17, 2)]) +
      geom_line(aes(colour = Design), size = 1, alpha=I(.5)) +
      facet_wrap(~ gamma.run) +
      scale_colour_discrete(name = "Design") + theme(legend.position="top") + 
      theme_bw()    
    pg
  }

edfPlot1.comCag = 
function(sim, rowNames = "Within RunBetweenCag:AniResidual", method = "REML", pTitle = "EDF plots") {
    
        tempS.bwAni = sim$tempS[grep(rowNames, rownames(sim$tempS)), 1:3]
  tempS.bwAni = cbind(tempS.bwAni, apply(sim$tempS[grep(rowNames, rownames(sim$tempS)), c(4,5)], 1, max))

    mm = melt(tempS.bwAni)
    
    #names(mm)[4] = "Method"
    
    mm$gamma.ani = as.numeric(as.character(mm$gamma.ani))
    # mm$gamma.run = as.factor( paste('gamma.run =', mm$gamma.run)) levels(mm$gamma.run ) = sort(levels(mm$gamma.run ))
    
    if(any(is.na(mm$gamma.run))){
    level.temp =  levels(mm$gamma.run)
    temp =  as.character(mm$gamma.run)
    temp[which(is.na(temp))] = "Fixed"
    mm$gamma.run = factor(temp,  levels =  c(level.temp, "Fixed"))
   }
   
   if(any(is.na(mm$gamma.cag))){
    level.temp =  levels(mm$gamma.cag)
    temp =  as.character(mm$gamma.cag)
    temp[which(is.na(temp))] = "Fixed"
    mm$gamma.cag = factor(temp,  levels =  c(level.temp, "Fixed"))
   }
 
    pg <- ggplot(mm, aes(x = log(gamma.ani), y = value)) + 
    ylab("EDF") + xlab(expression(paste(sigma[A]^2, " / ", sigma^2))) + 
    ggtitle(pTitle) + 
    scale_x_continuous(breaks = log(unique(mm$gamma.ani))[seq(1, 17, 2)], 
            labels = as.character(unique(round(mm$gamma.ani, 4)))[seq(1, 17, 2)]) + 
    geom_line(aes(colour = gamma.cag), size = 1, alpha = I(0.5)) + 
    facet_wrap(~gamma.run) + 
    scale_colour_discrete(name = "gamma.cag")     + theme(legend.position="top")
    pg
}
 edfPlot2.comCag = function(sim, rowNames = "Within RunBetweenCag:AniResidual", method = "REML", pTitle = "EDF plots") {
    
        tempS.bwAni = sim$tempS[grep(rowNames, rownames(sim$tempS)), 1:3]
  tempS.bwAni = cbind(tempS.bwAni, apply(sim$tempS[grep(rowNames, rownames(sim$tempS)), c(4,5)], 1, max))

    mm = melt(tempS.bwAni)
    
    #names(mm)[4] = "Method"
    
    mm$gamma.ani = as.numeric(as.character(mm$gamma.ani))
    # mm$gamma.run = as.factor( paste('gamma.run =', mm$gamma.run)) levels(mm$gamma.run ) = sort(levels(mm$gamma.run ))
    
    if(any(is.na(mm$gamma.run))){
    level.temp =  levels(mm$gamma.run)
    temp =  as.character(mm$gamma.run)
    temp[which(is.na(temp))] = "Fixed"
    mm$gamma.run =  factor(temp,  levels =  c(level.temp, "Fixed"))
   }
   
   if(any(is.na(mm$gamma.cag))){
    level.temp =  levels(mm$gamma.cag)
    temp =  as.character(mm$gamma.cag)
    temp[which(is.na(temp))] = "Fixed"
    mm$gamma.cag = factor(temp,  levels =  c(level.temp[1:5], "Fixed", level.temp[6]))
   }
 
    pg <- ggplot(mm, aes(x = log(gamma.ani), y = value)) + 
    ylab("EDF") + xlab(expression(paste(sigma[A]^2, " / ", sigma^2))) + 
    ggtitle(pTitle) + 
    scale_x_continuous(breaks = log(unique(mm$gamma.ani))[seq(1, 17, 2)], 
            labels = as.character(unique(round(mm$gamma.ani, 4)))[seq(1, 17, 2)]) + 
    geom_line(aes(colour = gamma.cag), size = 1, alpha = I(0.5)) + 
    facet_wrap(~gamma.run) + 
    scale_colour_discrete(name = "gamma.cag")     + theme(legend.position="top")
  
    
    pg
}
 
 
 
 
 
 edfPlotCompare1 =
function(sim, sim1, sim2, sim3, compare, rowNames = "Within RunBetweenAniResidual", pTitle = "EDF plots"){

 
   tempS.bwAni =
    cbind(sim[is.na(sim$gamma.cag), c(1, 3)],
      rbind(sim[is.na(sim$gamma.cag) & !is.na(sim$gamma.run), c(4,5)], 
            sim[is.na(sim$gamma.cag) & is.na(sim$gamma.run), c(4,5)]),
      rbind(sim1[is.na(sim1$gamma.cag) & !is.na(sim1$gamma.run), c(4,5)], 
            sim1[is.na(sim1$gamma.cag) & is.na(sim1$gamma.run), c(4,5)]),
      rbind(sim2[is.na(sim2$gamma.cag) & !is.na(sim2$gamma.run), c(4,5)], 
            sim2[is.na(sim2$gamma.cag) & is.na(sim2$gamma.run), c(4,5)]),
      rbind(sim3[is.na(sim3$gamma.cag) & !is.na(sim3$gamma.run), c(4,5)], 
            sim3[is.na(sim3$gamma.cag) & is.na(sim3$gamma.run), c(4,5)]))
    

  
    colnames(tempS.bwAni) [3:10] = paste0(rep(compare, each = 2), "/", colnames(tempS.bwAni) [3:10])
    mm = melt(tempS.bwAni)  
    
    colnames(mm) [3] = "Design"
    
   mm$gamma.ani =   as.numeric(as.character(mm$gamma.ani))
    mm$gamma.run =   ifelse(
      is.na(mm$gamma.run),
      "plain(\"Runs Fixed\")",
      paste("sigma[r]^2/sigma^2 ==", mm$gamma.run)
    )
    
    mm$gamma.run = factor(mm$gamma.run, levels =  c( "sigma[r]^2/sigma^2 == 0", 
                                                     "sigma[r]^2/sigma^2 == 0.25",
                                                     "sigma[r]^2/sigma^2 == 1", 
       "sigma[r]^2/sigma^2 == 4", "sigma[r]^2/sigma^2 == 100" , "plain(\"Runs Fixed\")"))

 
    
    pg <- ggplot(mm, aes(x = log(gamma.ani), y=value))+ ylab("EDF") +
      xlab(expression(paste(sigma[p]^2," / ",sigma^2))) +
      ggtitle(pTitle) +
      scale_x_continuous( breaks = log(unique(mm$gamma.ani))[seq(1, 17, 2)],
                          labels = as.character(unique(round(mm$gamma.ani, 4)))[seq(1, 17, 2)]) +
      geom_line(aes(colour = Design), size = 1, alpha=I(.5)) +
      facet_wrap(c("gamma.run"), labeller = label_parsed) +
       scale_colour_discrete(name = "Design/Method")    + theme(legend.position = "top") + 
      theme_bw()  
    
      pg
 
}
   
   
   
    