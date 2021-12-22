fracToNum =
function(x){
  sapply(strsplit(x, "/"), function(s) ifelse(length(s) ==1, as.numeric(s),
            as.numeric(s)[1]/as.numeric(s)[2]))
  }


obj.fun.ave = function(ind) {


  Rep = Z1.rep
    newZ1 = Z1.mat[ind[-length(ind)], ]
    info.mat = t(newZ1) %*% blk.proj %*% newZ1
    e.va = eigen(info.mat, only.values = TRUE)$va
    can.eff = e.va[-which(e.va < 1e-07)]/Rep

    aveEff = 100/mean(1/can.eff)


    return(aveEff)
}

obj.fun.ms = function(ind) {

    newZ1 = Z1.mat[ind[-length(ind)], ]
    NaR = t(newZ1) %*% Zb
    NaT = t(newZ1) %*% Zt

    W =  NaR %*% t(NaR)/(nPlot * ani.Rep) + NaT %*% t(NaT)/(nBlk * ani.Rep)

    aveEff = tr(W %*% W)

    return(-aveEff)
}

obj.fun.penalised.ms  = function(ind, w) {

    newZ1 = Z1.mat[ind[-length(ind)], ]
   info.mat = t(newZ1) %*% blk.proj %*% newZ1

    #aveEff = tr(info.mat) + w/tr(info.mat %*% info.mat)
    
    aveEff = tr(info.mat) -  sum((info.mat[upper.tri(matrix(0,ncol = ncol(info.mat),nrow = nrow(info.mat)))])^2)/w
    
    return(aveEff)
}

obj.fun.m = function(ind) {

    newZ1 = Z1.mat[ind[-length(ind)], ]
    info.mat = t(newZ1) %*%(mI(n) - Pb) %*% newZ1

    aveEff = tr(info.mat)

    return(aveEff)
}

obj.fun.m = function(ind,w) {

    newZ1 = Z1.mat[ind[-length(ind)], ]
    info.mat = t(newZ1) %*% blk.proj %*% newZ1

    aveEff = tr(info.mat)

    return(aveEff)
}


obj.fun.s = function(ind,w) {

    #newZ1 = Z1.mat[ind[-length(ind)], ]
    #info.mat = t(newZ1) %*% blk.proj %*% newZ1

    #aveEff = tr(info.mat %*% info.mat)

    newZ1 = Z1.mat[ind[-length(ind)], ]
    
    info.mat = t(newZ1) %*% blk.proj %*% newZ1
    
    return(-tr(info.mat %*% info.mat))
}


swap <- function(ind) {

    idx <- seq(2, length(ind) - 2)
    changepoints <- sample(idx, size = 2, replace = FALSE)
    # Two checks for omitting the swap of two indtical animals and animals in the same block
    check1 = all(Z1.mat[ind[changepoints[1]], ] == Z1.mat[ind[changepoints[2]], ])

    while (check1) {
        changepoints <- sample(idx, size = 2, replace = FALSE)
        check1 = all(Z1.mat[ind[changepoints[1]], ] == Z1.mat[ind[changepoints[2]], ])
    }
    tmp <- ind[changepoints[1]]
    ind[changepoints[1]] <- ind[changepoints[2]]
    ind[changepoints[2]] <- tmp
    return(ind)
}


multiSA = function(newInit, temp, tmax, obj.fun) {

    top = 250
    n = 100000
    best250 = matrix(0, nrow = top, ncol = length(newInit))
    best250.obj = numeric(top)
    best250.obj = rep(-Inf,top)

    #best250[1,] = newInit
    #best250.obj[1] = obj.fun(newInit,temp)
    
    prevInit = newInit
    
    for (i in 1:n) {
        currInit = swap(prevInit)
        currTemp = temp/log(((i - 1)%/%tmax) * tmax + exp(1))
        previous = obj.fun(prevInit,currTemp)
        current = obj.fun(currInit,currTemp)
        
        if(current > previous){
          prevInit = currInit

        } else{
         compare = exp((current - previous)/currTemp)

         if(compare > runif(1,0,1)){
           prevInit = currInit
         } else{
           prevInit = prevInit
         }
        }
        
        #if(i > (n - tmax)){
          if(any(current >  best250.obj)){  #saving the top250 designs
              best250[which(current >  best250.obj)[1],] = currInit
              best250.obj[which(current >  best250.obj)[1]] = current
          }
        #}
    }
    
    return(best250)
}

multiSA = cmpfun(multiSA)

old.time = proc.time()

init = sample(1:n)

newInit = c(init, init[1])

obj.fun.penalised.ms(newInit,1)

U = numeric(100000)
temp = newInit
for(i in 1:100000){
  temp = swap(temp)
  U[i] = obj.fun.m(temp)
}

max(U) - min(U)

optimised.m = multiSA(newInit, 100, 1000)

optimised.m = multiSA(newInit, 1, 10000, obj.fun.penalised.ms)
currTemp = 1/log(((100000 - 1)%/%10000) * 10000 + exp(1))

optimised.m = multiSA(newInit, 100, 1000, obj.fun.m)

proc.time() - old.time

obj.fun.penalised.ms(optimised.m[1,], 1)

optimised.m[which(apply(optimised.m, 1, obj.fun.s) == max(apply(optimised.m, 1, obj.fun.s))),]

bestInd = optimised.m[which(apply(optimised.m, 1, obj.fun.s) == max(apply(optimised.m, 1, obj.fun.s)))[1],-13]


apply(optimised.m[which(apply(optimised.m, 1, obj.fun.s) == max(apply(optimised.m, 1, obj.fun.s))),], 1, obj.fun.ave )

index = which((apply(optimised.m, 1, obj.fun.m ) - currTemp/apply(optimised.m, 1, obj.fun.s)) == max(apply(optimised.m, 1, obj.fun.m ) - currTemp/apply(optimised.m, 1, obj.fun.s)))
index = 1:250

apply(optimised.m[index,], 1, obj.fun.ave )
apply(optimised.m[index,], 1, obj.fun.m )
apply(optimised.m[index,], 1, obj.fun.s)
apply(optimised.m[index,], 1, obj.fun.penalised.ms, w = currTemp)


which(apply(optimised.m[index,], 1, obj.fun.penalised.ms, w = 1) == 
max(apply(optimised.m[index,], 1, obj.fun.penalised.ms, w = 1)))

t(optimised.m[,-13])
bestInd = optimised.m[54,-13]


optimised.ms = multiSA(newInit, 100, 1000, obj.fun = obj.fun.ms)
apply(optimised.ms, 1, obj.fun.ave )
apply(optimised.ms, 1, obj.fun.m )
apply(optimised.ms, 1, obj.fun.s)


bestInd = optimised.m[1,-13]

write.csv(t(optimised.ms[,-13]), file= "RowColumnMSOpitmal.csv")


MSoptimal = read.csv(file.choose())
MSoptimal = t(optimised.m[index,-13])

MSoptimal = sapply(strsplit(unique(apply(MSoptimal[1:12,-1], 2, function(x) paste(x, sep ="", collapse = "."))), "\\."), function(x) x)

apply(MSoptimal, 2, function(x) obj.fun.ave(as.numeric(c(x,1))) )
apply(MSoptimal, 2, function(x) obj.fun.m(as.numeric(c(x,1))) )
apply(MSoptimal, 2, function(x) obj.fun.s(as.numeric(c(x,1))) )


##################################################################################
bestInd = as.numeric(MSoptimal[,13])

fractions(apply(MSoptimal, 2, function(x) test(X.trt = Z1.mat[as.numeric(x),], mI(n) - Pb, Rep=Z1.rep)$ave.eff)

fractions(apply(MSoptimal, 2, function(x) test(X.trt = Z1.mat[as.numeric(x),], Pb, Rep=Z1.rep)$ave.eff))

################################################################################
#Test tag test
################################################################################

testTag = function(bestInd){

  new.Z1.mat =  Z1.mat[bestInd,]
  #new.Z1.mat=  Z1.mat

################################################################################
  #Set the design from the search

  colnames(new.Z1.mat) = sort(levels(interaction(newLETTERS[1:nZ1])))
  new.Z1.des = apply(new.Z1.mat, 1,  function(x)
  colnames(new.Z1.mat)[which(as.logical(x))])
  new.Z1.des = as.data.frame(t(sapply(strsplit(new.Z1.des, "\\."), function(x) x)))

  design.df = data.frame(Run = as.factor(Z.des), Tag = factor(1:nPlot))

  design.df = cbind(design.df,
              Ani = t(new.Z1.des),
              Trt = phase1DesignEX1[match(as.character(t(new.Z1.des )),
                            as.character(phase1DesignEX1$Ani)),]$Trt)

################################################################################
  #Theortical ANOVA table


  ss = summary.aov.twoPhase(design.df,  blk.str2 = "Run", blk.str1 = "Ani", trt.str = "Tag + Trt")

  rNames = rownames(ss$ANOVA)

  diff(as.numeric(ss$ANOVA[(rev(which(rNames == "   Residual"))[1]+1):length(rNames),"e"]))==0
}

apply(MSoptimal, 2, function(x) testTag(as.numeric(x)))

################################################################################
#Tag Average Efficiency Factor
################################################################################

aveTag = function(bestInd){
  new.Z1.mat =  Z1.mat[bestInd,]
  #new.Z1.mat=  Z1.mat

################################################################################
  #Set the design from the search

  colnames(new.Z1.mat) = sort(levels(interaction(newLETTERS[1:nZ1])))
  new.Z1.des = apply(new.Z1.mat, 1,  function(x)
  colnames(new.Z1.mat)[which(as.logical(x))])
  new.Z1.des = as.data.frame(t(sapply(strsplit(new.Z1.des, "\\."), function(x) x)))

  design.df = data.frame(Run = as.factor(Z.des), Tag = factor(1:nPlot))

  design.df = cbind(design.df,
              Ani = t(new.Z1.des),
              Trt = phase1DesignEX1[match(as.character(t(new.Z1.des )),
                            as.character(phase1DesignEX1$Ani)),]$Trt)

################################################################################
  #Theortical ANOVA table


  #test tag possible
  ss = summary.aov.twoPhase(design.df,  blk.str2 = "Run", blk.str1 = "Ani", trt.str = "Tag + Trt")

  rNames = rownames(ss$EF)
  fracToNum(ss$EF[grep("Tag", rNames),"eff.Tag"])
}

apply(MSoptimal, 2, function(x) aveTag(as.numeric(x)))


################################################################################
#Test treatment test
################################################################################

testTrt = function(bestInd){
  new.Z1.mat =  Z1.mat[bestInd,]
  #new.Z1.mat=  Z1.mat

################################################################################
  #Set the design from the search

  colnames(new.Z1.mat) = sort(levels(interaction(newLETTERS[1:nZ1])))
  new.Z1.des = apply(new.Z1.mat, 1,  function(x)
  colnames(new.Z1.mat)[which(as.logical(x))])
  new.Z1.des = as.data.frame(t(sapply(strsplit(new.Z1.des, "\\."), function(x) x)))

  design.df = data.frame(Run = as.factor(Z.des), Tag = factor(1:nPlot))

  design.df = cbind(design.df,
              Ani = t(new.Z1.des),
              Trt = phase1DesignEX1[match(as.character(t(new.Z1.des )),
                            as.character(phase1DesignEX1$Ani)),]$Trt)

################################################################################
  #Theortical ANOVA table


  #test tag possible
  ss = summary.aov.twoPhase(design.df,  blk.str2 = "Run", blk.str1 = "Ani", trt.str = "Tag + Trt")

  rNames = rownames(ss$ANOVA)
  btAni = ss$ANOVA[(rev(which(rNames == "   Between Ani"))[1]+1):(rev(which(rNames == "   Residual"))[1]-1), c("e", "Ani")]

  rNames = rownames(btAni)
  EF = try(btAni[c(grep("Trt", rNames), grep("Residual", rNames)),], TRUE)

  if(is.null(nrow(EF))) {
    return("Confounded with Tag!")
  } else{
    EF = apply(EF, 1, function(x) paste(x, sep = "", collapse = "."))
    return(EF[1]== EF[2])
  }
}

apply(MSoptimal, 2, function(x) testTrt(as.numeric(x)))


################################################################################
#Test treatment average efficiency factors
################################################################################

aveTrt = function(bestInd){
  new.Z1.mat =  Z1.mat[bestInd,]
  #new.Z1.mat=  Z1.mat

################################################################################
  #Set the design from the search

  colnames(new.Z1.mat) = sort(levels(interaction(newLETTERS[1:nZ1])))
  new.Z1.des = apply(new.Z1.mat, 1,  function(x)
  colnames(new.Z1.mat)[which(as.logical(x))])
  new.Z1.des = as.data.frame(t(sapply(strsplit(new.Z1.des, "\\."), function(x) x)))

  design.df = data.frame(Run = as.factor(Z.des), Tag = factor(1:nPlot))

  design.df = cbind(design.df,
              Ani = t(new.Z1.des),
              Trt = phase1DesignEX1[match(as.character(t(new.Z1.des )),
                            as.character(phase1DesignEX1$Ani)),]$Trt)

################################################################################
  #Theortical ANOVA table
  aveEff = as.character(rep("0", 3))

  #test tag possible
  ss = summary.aov.twoPhase(design.df,  blk.str2 = "Run", blk.str1 = "Ani", trt.str = "Tag + Trt")

  rNames = rownames(ss$EF)

  EF = ss$EF[,"eff.Trt"]
  
  EF = EF[-which(EF == "")]
  
  if(length(EF) == 3){
    aveEff = EF

  } else if( length(EF) ==2){
     if(grepl("Tag", names(EF)[1])){
        aveEff[2] = EF[1]
     }
     if(grepl("Trt", names(EF)[1])){
        aveEff[1] = EF[1]
     }
     if(grepl("Tag", names(EF)[2])){
        aveEff[2] = EF[2]
     }
     if(grepl("Trt", names(EF)[2])){
        aveEff[3] = EF[2]
     }

  } else if(grepl("Tag", names(EF))){
     aveEff[2] = EF
  
  } else if(grepl("Trt", names(EF))){
     aveEff[3] = EF
  }
  
  return(fracToNum(aveEff))
}

apply(MSoptimal, 2, function(x) aveTrt(as.numeric(x)))


 MSoptimal = read.csv(file.choose())
MSoptimal = sapply(strsplit(unique(apply(MSoptimal[,-c(1,2)], 2, function(x) paste(x, sep ="", collapse = "."))), "\\."), function(x) x)

apply(MSoptimal, 2, function(x) obj.fun.m(as.numeric(c(x,1))))
apply(MSoptimal, 2, function(x) obj.fun.s(as.numeric(c(x,1))))
apply(MSoptimal, 2, function(x) test(X.trt = Z1.mat[as.numeric(x),], Pb, Rep=Z1.rep)$ave.eff)
apply(MSoptimal, 2, function(x) test(X.trt = Z1.mat[as.numeric(x),], mI(n) - Pb, Rep=Z1.rep)$ave.eff)
apply(MSoptimal, 2, function(x) testTag(as.numeric(x)))
apply(MSoptimal, 2, function(x) aveTag(as.numeric(x)))
apply(MSoptimal, 2, function(x) testTrt(as.numeric(x)))
apply(MSoptimal, 2, function(x) aveTrt(as.numeric(x)))
                
write.csv(rbind(MSoptimal,
                apply(MSoptimal, 2, function(x) obj.fun.m(as.numeric(c(x,1)))),
                apply(MSoptimal, 2, function(x) obj.fun.s(as.numeric(c(x,1)))),
                apply(MSoptimal, 2, function(x) test(X.trt = Z1.mat[as.numeric(x),], Pb, Rep=Z1.rep)$ave.eff),
                apply(MSoptimal, 2, function(x) test(X.trt = Z1.mat[as.numeric(x),], mI(n) - Pb, Rep=Z1.rep)$ave.eff),
                apply(MSoptimal, 2, function(x) testTag(as.numeric(x))),
                apply(MSoptimal, 2, function(x) aveTag(as.numeric(x))),
                apply(MSoptimal, 2, function(x) testTrt(as.numeric(x))),
                apply(MSoptimal, 2, function(x) aveTrt(as.numeric(x)))), "S Opitmal result.csv")

obj.fun.m = function(ind) {

    newZ1 = Z1.mat[ind[-length(ind)], ]
    info.mat = t(newZ1) %*% (mI(n) - Pb) %*% newZ1

    aveEff = tr(info.mat)

    return(aveEff)
}

apply(MSoptimal, 2, function(x) obj.fun.m(as.numeric(c(x,1))))


