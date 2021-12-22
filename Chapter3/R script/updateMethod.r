#19/04/2012 11:15:34 a.m.


sourceDir <- function(path, trace = TRUE, ...) {
  library(MASS)
  library(lattice)
  for (nm in list.files(path, pattern = "\\.[RrSsQq]$")) {
    if(trace) cat(nm,":")
      source(file.path(path, nm), ...)
      if(trace) cat("\n")
  }
}


sourceDir(path = "C:/Users/kcha193/My Dropbox/R functions for two-phase experiments")



phase1DesignEX1 <- local({
  Ani = as.factor(LETTERS[1:6])
  Trt = as.factor(c("Con", "Dis", "Trt")[rep(1:3, each = 2)])
  data.frame(Ani,Trt)
})
phase1DesignEX1

phase2designEX1 <- local({
  Run = as.factor(rep(1:4, each = 6))
  Ani = as.factor(LETTERS[c(1,2,3,4,5,6,
                            2,3,4,5,6,1,
                            4,5,6,1,2,3,
                            5,6,1,2,3,4)])
  Tag = as.factor(c(114:119)[rep(1:6, 4)])
  Trt = phase1DesignEX1$Trt[match(Ani, phase1DesignEX1$Ani)]
  data.frame(Run, Ani, Tag, Trt)
})
phase2designEX1

summary.aov.twoPhase(phase2designEX1,  blk.str2 = "Run" , blk.str1 = "Ani",
trt.str = "Trt + Tag", response = rnorm(24))


d = phase2designEX1

d$Ani[b] = phase2designEX1$Ani[a]
d$Ani[a] = phase2designEX1$Ani[b]

summary.aov.twoPhase(d,  blk.str2 = "Run" , blk.str1 = "Ani",
trt.str = "Trt + Tag", response = rnorm(24))

Pa = Pb1$Within$Ani +  Pb1$Within$Residual

T =  makeBlkDesMatrix(phase2designEX1, "Ani")$A

Ab = t(T) %*% Pa  %*% T
eigen(Ab)$va/4


T1 = T

a = 1
b = 24

x1 = T1[b,] = T[a,]
x2 = T1[a,] = T[b,]

c1 = Pa[a,a]
c2 = Pa[b,b]
  
d1 = Pa[a, -c(a,b)]
d2 = Pa[b, -c(a,b)]

T2 = T[-c(a,b),]

B = (x1 - x2) %*% t(d1 - d2) %*% T2

Aa = t(T1) %*% Pa  %*% T1
eigen(Aa)$va/4

eigen(Ab - Aa)


identical(Ab - Aa, as.numeric(c1 - c2) * (x1 %*% t(x1) - x2 %*% t(x2)) + B + t(B))

eigen(as.numeric(c1 - c2) * (x1 %*% t(x1) - x2 %*% t(x2)))

eigen(B)
eigen(t(B))


eigen(as.numeric(c1 - c2) * (x1 %*% t(x1) - x2 %*% t(x2)) + B + t(B))

t = t(d1 - d2) %*% T2

a1 = t[1] - t[2]

a2 = c1- c2 + t[1] + t[2]

e = sqrt(a1^2 + a2^2 + as.numeric(2* t %*% t(t)))

c(a1 - e, a1+e)




#Compute the eigenvalues



update.Infor.mat.eigen =
function(int, trt.mat, blk.proj){

  library(MASS)

  Ab = t(trt.mat) %*% blk.proj %*% trt.mat
  ginv.Ab = ginv(Ab)
  Rep = nrow(trt.mat)/ncol(trt.mat)
  
 # if(a[1] == a[2]) return()

  #if((a[1] %% 1) !=0 || (a[2] %% 1) !=0) return()

  a = int[1]
  b = int[2]
  
  T1 = trt.mat

  x1 = T1[b,] = trt.mat[a,]
  x2 = T1[a,] = trt.mat[b,]

  if(isTRUE(all.equal(x1,x2))) return(NULL)
  
  c1 = blk.proj[a,a]
  c2 = blk.proj[b,b]

  d1 = blk.proj[a, -c(a,b)]
  d2 = blk.proj[b, -c(a,b)]

  T2 = trt.mat[-c(a,b),]

  t = t(d1 - d2) %*% T2

  #return( list( n.trt.mat = T1,
  #A = (t(trt.mat) %*% blk.proj %*% trt.mat) -
  #  ((c1-c2) * (x1%*% t(x1) -x2%*%t(x2)) + B + t(B))))

  a = which(as.logical(trt.mat[a,]))
  b = which(as.logical(trt.mat[b,]))

  t1 = t[,a]
  t2 = t[,b]

  a1 = t1- t2
  a2 = c1 - c2 +t1 +t2

  e = sqrt(a1^2 + a2^2 + 2* sum(t[,-c(a,b)]^2))

  e1 = a1 - e
  e2 = a1 + e

  e.value = c(e1, e2)
  
  if(length(which(e.value<=0))==0){
    e.value = e.value
  } else{
    e.value = e.value[-which(e.value<=0)]
  }
  if(length(e.value) == 0) return( NULL)
  
  e.vector = numeric()
 
  #Compute the eigenvectors
  for( i in 1:length(e.value)){
    q1 = (a1*e.value[i] + sum(t[,-c(a,b)]^2))/( e.value[i] - a2)
    q2 = q1 - e.value[i]

    p = t
    p[a] = q1
    p[b] = q2

    p = p/(max(p) * sqrt(2))

     e.vector = rbind(e.vector, p)
  }
  
  e.vector = t(e.vector)
  
  if(length(e.value) == 1){
    psi = (1/e.value) - t(e.vector) %*% ginv.Ab %*%  e.vector
  }else{
    psi = diag(1/e.value) - t(e.vector) %*% ginv.Ab %*%  e.vector
  }
  ginv.Aa = ginv.Ab + ginv.Ab %*%  e.vector %*% ginv(psi)  %*%  t(e.vector) %*% ginv.Ab

  e.va = eigen(ginv(ginv.Aa))$va
  Aa.ave.eff = (length(e.va[-which(e.va<1e-7)])/sum(1/e.va[-which(e.va<1e-7)]))/Rep

  e.va = eigen(Ab)$va
  Ab.ave.eff = (length(e.va[-which(e.va<1e-7)])/sum(1/e.va[-which(e.va<1e-7)]))/Rep


  if(Ab.ave.eff>  Aa.ave.eff){
    return (NULL)
  }else{
    return(list(ginvAb = ginv.Ab, AbAveEff = Ab.ave.eff,
                newTrtMat = T1,
                ginvAa = ginv.Aa, AaAveEff = Aa.ave.eff))
              
  }
}



Ab - Aa

a = 1
b = 7

nTrt = 8
Rep = 3

trt.des = as.factor(rep(1:nTrt, Rep))
X.trt = matrix(0, ncol = nTrt, nrow= nTrt*Rep)
X.trt[cbind(1:( nTrt*Rep), trt.des)] = 1

trt.des = sample(trt.des)
X.trt = matrix(0, ncol = nTrt, nrow= nTrt*Rep)
X.trt[cbind(1:( nTrt*Rep), trt.des)] = 1

phase2designEX1$temp = trt.des

 summary.aov.twoPhase(phase2designEX1,  blk.str2 = "Run" , blk.str1 = "Ani",
trt.str = "temp", response = rnorm(24))


t(X.trt) %*% blk.proj %*% X.trt

states = lapply(1:nrow(X.trt), function(x) cbind(x, (x+1):nrow(X.trt)))
states = states[-length(states)]

bestSwap = numeric(2)
best.ave.eff = 0
best.trt.des = as.factor(numeric(length(trt.des)))

for( i in 1:length(states)){
   S = states[[i]]
  for( j in 1:nrow(S)){
    A = try(update.Infor.mat(int = S[j,], trt.mat = X.trt, blk.proj = blk.proj), TRUE)
    
    if(class(A) =="try-error" || is.na(A$Ab.ave.eff) || is.na(A$Aa.ave.eff)) next
    
    cat("Old Average Efficiency factors =", A$Ab.ave.eff, "\n")
    
    if(best.ave.eff < A$Aa.ave.eff){
      X.trt = A$new.trt.mat
      best.ave.eff =A$Aa.ave.eff
      bestSwap = S[j,]
      cat("New Average Efficiency factors =", best.ave.eff, "\n")
    }
  }
}

phase2designEX1$temp = as.factor(apply(X.trt, 1, function(x) which(as.logical(x))))

 summary.aov.twoPhase(phase2designEX1,  blk.str2 = "Run" , blk.str1 = "Ani",
trt.str = "temp", response = rnorm(24))


#new way
trt.des = as.factor(rep(1:nTrt, each = Rep))
X.trt = matrix(0, ncol = nTrt, nrow= nTrt*Rep)
X.trt[cbind(1:( nTrt*Rep), trt.des)] = 1
best.ave.eff = 0
best.trt.des = as.factor(numeric(length(trt.des)))
#  pb <- txtProgressBar(min = 0, max = 100000, style = 3)

for( i in 1:100000){
   #setTxtProgressBar(pb, i)
  infor.mat = t(X.trt) %*% blk.proj %*% X.trt

  con.eff = eigen(infor.mat)$va
  con.eff = con.eff[which(con.eff>1e-7)]
  new.ave.eff= (length(con.eff)/sum(con.eff))/Rep
  
  if(is.na(new.ave.eff)) next
  
  if(new.ave.eff > best.ave.eff){
    best.trt.des = trt.des
    best.ave.eff = new.ave.eff
    cat("Average Efficiency factors =", best.ave.eff, "\n")
  }
  
  trt.des = sample(trt.des)
  X.trt = matrix(0, ncol = nTrt, nrow= nTrt*Rep)
  X.trt[cbind(1:( nTrt*Rep), trt.des)] = 1
}
#   close(pb)


A = try(update.Infor.mat(int = c(1,2), trt.mat = X.trt, blk.proj = blk.proj), TRUE)

eigen(t(X.trt) %*% blk.proj %*% X.trt)$va/Rep

phase2designEX1$temp = as.factor(apply(X.trt, 1, function(x) which(as.logical(x))))

 summary.aov.twoPhase(phase2designEX1,  blk.str2 = "Run" , blk.str1 = "Ani",
trt.str = "temp", response = rnorm(24))

 summary.aov.onePhase(phase2designEX1,  blk.str = "Run" ,,
trt.str = "temp", response = rnorm(24))




blk.proj = Pb1$Within$Residual

nTrt = 8
Rep = 24/nTrt


X.trt = matrix(0, ncol = nTrt, nrow= nTrt*Rep)

old.temp = list(X.trt = X.trt)
temp = list()

addPt =
function(x, index){
  x[index[1], index[2]] = 1
  return(x)
}

extractTMatrix =
function(y){

  trace.screen = sapply(y,  function(x) tr(t(x) %*% blk.proj %*% x))
  y[which(trace.screen == max(trace.screen))]

}

extractTMatrix1 =
function(y){

  trace.screen = sapply(y,  function(x) tr((t(x) %*% blk.proj %*% x)^2))
  y[which(trace.screen == min(trace.screen))]

}


for( i in 1:nrow(X.trt)){
  print(i)
  for(j in 1:ncol(X.trt)){
      temp.X.trt = old.temp$X.trt

      if(!is.list(temp.X.trt)){
        temp$X.trt[[j]] = addPt(temp.X.trt,  index = c(i,j))
      }else{
        temp$X.trt[[j]] =lapply(temp.X.trt, addPt, index = c(i,j))
        
      }

  }

    if(i != 1){

      newtemp = lapply(temp$X.trt, function(y) extractTMatrix(y))
      newtemp = lapply(newtemp, function(y) extractTMatrix1(y))

      temp = list()
      counter = 1
      for( k in 1:length(newtemp)){
       for(l in 1:length(newtemp[[k]])){

        temp$X.trt[[counter]] = newtemp[[k]][[l]]
        counter = counter + 1
       }
     }
    }

  temp$trace.screen =
   unlist(lapply(temp$X.trt, function(x) tr(t(x) %*% blk.proj %*% x)))

  if(length(temp$trace.screen) >10000){
    index = sample(10000)
    temp$trace.screen = temp$trace.screen[index[1:5000]]
    temp$X.trt = temp$X.trt[index[1:5000]]
  }
  
  if(length(temp$trace.screen == max(temp$trace.screen))>1){
    temp$X.trt = temp$X.trt[which(temp$trace.screen == max(temp$trace.screen))]
    temp$trace.screen = temp$trace.screen[which(temp$trace.screen == max(temp$trace.screen))]
  }
  old.temp = temp
  temp = list()

}


lapply(old.temp$X.trt,  colSums)


X.trt = old.temp$X.trt[[1]]

phase2designEX1$temp = as.factor(apply(X.trt, 1, function(x) which(as.logical(x))))

 summary.aov.twoPhase(phase2designEX1,  blk.str2 = "Run" , blk.str1 = "Ani",
trt.str = "temp", response = rnorm(24))


ave.eff=
function(x, blk.proj, Rep){
  e.va =  eigen(t(x) %*% blk.proj %*% x)$va

  (length(e.va[-which(e.va<1e-7)])/sum(1/e.va[-which(e.va<1e-7)]))/Rep
}

old.temp$trace.screen =
   sapply(old.temp$X.trt, function(x) ave.eff(x, blk.proj = blk.proj, Rep = Rep))

if(length(old.temp$trace.screen == max(old.temp$trace.screen))>1){
    old.temp$X.trt =  old.temp$X.trt[which(old.temp$trace.screen == max(old.temp$trace.screen))]
    old.temp$trace.screen = old.temp$trace.screen[which(old.temp$trace.screen == max(old.temp$trace.screen))]
  }

X.trt = old.temp$X.trt[[1]]

phase2designEX1$temp = as.factor(apply(X.trt, 1, function(x) which(as.logical(x))))

 summary.aov.twoPhase(phase2designEX1,  blk.str2 = "Run" , blk.str1 = "Ani",
trt.str = "temp", response = rnorm(24))

