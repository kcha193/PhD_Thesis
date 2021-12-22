

################################################################################
#Step 1:
#Extact the mean squares and degrees of freedom form the aov function
#These that means squares and degrees of freedoms are used to construct the initial
#information matrix


getMSEst =
function(response, trt.str,  blk.str1, blk.str2, design ){
  #Need all the mean squares and degrees of freedom from this aonva model,
  #because only the random effects are fitted here


  #First Phase, extract ALL information apart from the "Within" statum"
  formula1 = as.formula(paste("response ~", blk.str1, " + Error(", blk.str2, ")", sep = ""))

  aov.fit1 = aov(formula1, data = design)

  termsNotNeed = c("(Intercept)", "Within")

  index = match(termsNotNeed, names(aov.fit1 ))

  termsNeed = names(aov.fit1 )[-index]

  MSS = numeric(0)
  DFF = numeric(0)
  sv.name = character(0)
  m = numeric(0)
  for(i in 1:length(termsNeed)){
    MS = summary(aov.fit1[[termsNeed[i]]])[[1]]$'Mean Sq'
    MS = MS[length(MS)]
    DF = summary(aov.fit1[[termsNeed[i]]])[[1]]$'Df'
    DF = DF[length(DF)]
     
    m = c(m, DF/(2*MS^2))
    MSS = c(MSS, MS)
    DFF = c(DFF, DF)
    svname = paste("Between ", names(aov.fit1[termsNeed[i]]), "Between" ,
                  rownames((summary(aov.fit1[[termsNeed[i]]]))[[1]]), sep = "")

    svname = gsub("BetweenResiduals", "Residuals", svname)
    svname = gsub("      ", "", svname)

    if(length(svname) == 1){
     svname = gsub("Residuals", "", svname)

    }
    
      sv.name = c(sv.name, svname)
  }

  #Second phase, extract all the information
  #apart from the information of the strata deinfed in the first phase
  formula2 = as.formula(paste("response ~", trt.str, " + Error(", blk.str2, "+", blk.str1,")", sep = ""))

  aov.fit2 = suppressWarnings(aov(formula2, data = design))

  termsNotNeed =  c("(Intercept)", names(aov.fit1)[-index])

  index = match(termsNotNeed, names(aov.fit2))

  termsNeed =  names(aov.fit2)[-index]

  #Exact the mean squares and DF from the last sources of variation of every stratum
  for(i in 1:length(termsNeed)){
    MS = summary(aov.fit2[[termsNeed[i]]])[[1]]$'Mean Sq'

    MS = MS[length(MS)]

    DF = summary(aov.fit2[[termsNeed[i]]])[[1]]$'Df'
    DF = DF[length(DF)]

    svnames = rownames((summary(aov.fit2[[termsNeed[i]]]))[[1]])
    if(names(aov.fit2[termsNeed[i]]) == "Within"){
       svnames = paste("Residual", svnames[length(svnames)], sep = "")
    }  else {
      svnames = paste("Between", names(aov.fit2[termsNeed[i]]), svnames[length(svnames)], sep = "")
    }
    m = c(m, DF/(2*MS^2))
    MSS = c(MSS, MS)
    DFF = c(DFF, DF)
    sv.name = c(sv.name, paste("Within", svnames, sep=""))

  }

  sv.name = gsub("Residuals", "Residual", sv.name)

  #information matrix with the mean squares

  return(data.frame(sv.name, DFF, MSS, m))
}

################################################################################
#Step 2:
#Extact the variance components structure for inforDecompuTE to contrsut matrix G
#Need to modify to avoid the use of the package

fracToNum =
function(x){
  sapply(strsplit(x, "/"), function(s) ifelse(length(s) ==1, as.numeric(s),
            as.numeric(s)[1]/as.numeric(s)[2]))
}

getGMat =
function(MS.likelihood, trt.str, blk.str1, blk.str2, design, ...){

  VC = getVCs.twoPhase(design, blk.str1 , blk.str2, trt.str, ...)$random

  VC.numeric = t(apply(VC, 1, fracToNum))

  r.names = rownames(VC.numeric)
  new.r.names = character(length(r.names))
  stratum.names = character(2)

  #Exact and organize the row names of the ANOVA table from inforDecompuTE
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

  rownames(VC.numeric) = new.r.names

  VC.numeric[as.character(MS.likelihood$sv.name), -1]

}


G.mat = getGMat(MS.likelihood = MS.likelihood, trt.str = "Trt + Tag", blk.str1 = "Ani",
blk.str2 = "Run", design = design)


################################################################################
#Step 3:

################################################################################
#Step 4:
#REML to obtain the better estimates of the parameters of interest
#Use some random inital esitmates of variance componenets
#Get the score functions and then use the inverse of information matrix
#Start the fisher's scoring algorithm

#This function is to find the independant row set of linear combination of the
#variance componenets

getVcEDF =
function(MS.likelihood, G.mat, real.VC = rep(NA, ncol(G.mat))){

  m.mat = diag(MS.likelihood$m)

  infor.mat = t(G.mat) %*% m.mat %*% G.mat

  n.infor.mat = infor.mat
  var.comp = rep(0.001, ncol(G.mat))      #initial guess
  old.var.comp = rep(0, ncol(G.mat))

  MSS = MS.likelihood$MSS
  DFF = MS.likelihood$DFF
  
  while(! isTRUE(all.equal(var.comp,old.var.comp))){
    old.var.comp = var.comp

    MSS.hat = G.mat %*% var.comp

    score.fun = t(G.mat) %*% (DFF * (MSS - MSS.hat) /(2 * MSS.hat^2))
    
    
    m.mat = diag(DFF/(2*as.numeric(MSS.hat)^2))

    n.infor.mat =  t(G.mat) %*% m.mat %*% G.mat

    var.comp  = old.var.comp  +  ginv(n.infor.mat) %*% score.fun
  
    
    #MSS = MSS.hat
    #print(var.comp)
    #var.comp[which(var.comp<0)] = 0 
  }

   rownames(var.comp) = colnames(G.mat)

  if(!all(is.na(real.VC))){
    real.MSS.hat = G.mat %*% real.VC  
    real.m.mat = diag(DFF/(2*as.numeric(real.MSS.hat)^2))
    real.infor.mat =  t(G.mat) %*% real.m.mat %*% G.mat
    inv.real.infor.mat = solve(real.infor.mat)
    r.EDF = numeric(length(MS.likelihood$m))
  }
  
  #Obtain the variances of all the paramters, than compute the EDF
  inv.infor.mat = solve(n.infor.mat)

  EDF = numeric(length(MS.likelihood$m))
  
  for(i in 1:nrow(G.mat)){
    index = which(!G.mat[i,]==0)

    var.mat = inv.infor.mat[index, index]

    coef.mat = outer(G.mat[i,index], G.mat[i,index], "*")

    coef.mat[upper.tri(coef.mat)] =  coef.mat[upper.tri(coef.mat)]
    coef.mat[lower.tri(coef.mat)] =  coef.mat[lower.tri(coef.mat)]

    EDF[i] = 2 * MSS.hat[i]^2/sum(coef.mat * var.mat)
    
    if(!all(is.na(real.VC))){
      real.var.mat = inv.real.infor.mat[index, index]
    
      r.EDF[i] = 2 * real.MSS.hat[i]^2/sum(coef.mat * real.var.mat)
    }

  }

  list(Stratum = cbind(MS.likelihood, EDF, r.EDF), Var.comp = var.comp)
}


################################################################################
library(MASS)
library(asreml)
library(infoDecompuTE)

design = phase2designEX4

getVCs.twoPhase(phase2designEX4, blk.str1 = "Ani", blk.str2 = "Run",
trt.str = "Trt + Tag")

EDF = numeric(3)
r.EDF = numeric(3)
VC = numeric(3)


x = c(0.0001, 0.001, 0.01, 0.1, 1, 10, 100, 1000, 10000)
n = length(x)

total <- n
# create progress bar
pb <- txtProgressBar(min = 0, max = total, style = 3)

for( i in 1:n){
  setTxtProgressBar(pb, i)
  
  gamma.run = 0.1
  gamma.ani = x[i]

  run.eff = rnorm(4, mean = 0, sd = sqrt(gamma.run * 1))
  ani.eff = rnorm(4, mean = 0, sd = sqrt(gamma.ani * 1))
  trt.eff = c(1, 2)
  tag.eff = c(0,0,0,0)
  res.eff = rnorm(16, mean = 0, sd = 1)

  real.VC = c(1, (gamma.ani * 1),(gamma.run * 1))
   
  y = with(design, run.eff[Run] + ani.eff[Ani] + tag.eff[Tag] + trt.eff[Trt]) + res.eff
  
  MS.likelihood =
  getMSEst(response = y, trt.str = "Trt + Tag", blk.str1 = "Ani",
  blk.str2 = "Run", design = design)
 
  MS.likelihood
  
  G.mat = getGMat(MS.likelihood = MS.likelihood, trt.str = "Trt + Tag", 
  blk.str1 = "Ani", blk.str2 = "Run", design = design)

 
  temp = getVcEDF(MS.likelihood = MS.likelihood, G.mat = G.mat, real.VC = real.VC)
  
  EDF = rbind(EDF, temp$S$E)
  r.EDF = rbind(r.EDF, temp$S$"r.EDF")
  VC =  rbind(VC, t(temp$V))

}

close(pb)

apply(EDF[-1,], 2, mean)
apply(r.EDF[-1,], 2, mean)
apply(VC[-1,], 2, mean)


opar = par(mfrow = c(1,3))
  hist(EDF[-1,1])
  hist(EDF[-1,2])
  hist(EDF[-1,3])
par(opar)


opar = par(mfrow = c(1,3))
  hist(VC[-1,1])
  hist(VC[-1,2])
  hist(VC[-1,3])
par(opar)







asreml.fit = asreml(y ~ Trt + Tag, random = ~ Run + Ani, data = design)

summary(asreml.fit)
svc.asreml(asreml.fit)



#With losing runs




################################################################################


design.df = data.frame(Set = factor(rep(1:4, each = 4)),
                       Array = factor(rep(1:2, time = 4, each = 2)),
                       Dye = factor(c("Gre", "Red")),
                       Treat = factor(LETTERS[c(1,2,2,1,
                                                1,2,2,1,
                                                1,2,2,1,
                                                1,2,2,1)]),
                       Plant = interaction(factor(rep(1:4, each = 4)),
                                           factor(LETTERS[c(1,2,2,1,
                                                            1,2,2,1,
                                                            1,2,2,1,
                                                            1,2,2,1)])))



array.eff = rnorm(2, mean = 0, sd = 9)
trt.eff = rep(0,2)
dye.eff = c(0,0)
res.eff = rnorm(16, mean = 0, sd = 1)
pla.eff = rnorm(8, mean = 0, sd = 4)

y = with(design.df , array.eff[Array] + trt.eff[Treat] +
dye.eff[Dye] + pla.eff[Plant]) + res.eff


summary.aov.modified(response = y, blk.str2 = "Set/Array", blk.str1 = "Plant",
                trt.str = "Treat + Dye", design = design.df)


MS.likelihood =
getMSLikelihoodEst(response = y, blk.str2 = "Set/Array", blk.str1 = "Plant",
                trt.str = "Treat + Dye", design = design.df)

G.mat = getGMat(MS.likelihood = MS.likelihood, blk.str2 = "Set/Array", blk.str1 = "Plant",
                trt.str = "Treat + Dye", var.comp = c("Plant", "Set:Array"), design = design.df)

getVcEDF(MS.likelihood = MS.likelihood, G.mat = G.mat)
