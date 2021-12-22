ave.eff=
function(x, blk.proj, Rep){
  e.va =  eigen(t(x) %*% blk.proj %*% x)$va
  if(all(e.va<1e-7)){
    (length(e.va)/sum(1/e.va))/Rep
  } else{
    (length(e.va[-which(e.va<1e-7)])/sum(1/e.va[-which(e.va<1e-7)]))/Rep
  }
}

bal.eff=
function(x, blk.proj, Rep){
  e.va =  eigen(t(x) %*% blk.proj %*% x)$va

  sd(e.va[-which(e.va<1e-7)])
}

addPt =
function(x, index){
  temp = rep(0, ncol(x))
  temp[index] = 1
  return(rbind(x, temp))
}

extractMS1=
function(y, index){

  trace.screen = sapply(y,  function(x) tr(t(x) %*% blk.proj[1:index, 1:index] %*% x))
  y[which(trace.screen == max(trace.screen))]

}

extractMS2 =
function(y, index){

  trace.screen = sapply(y,  function(x) tr((t(x) %*% blk.proj[1:index, 1:index] %*% x)^2))
  y[which(trace.screen == min(trace.screen))]

}

c.optimised=
function(X.trt, blk.proj, samp = 10000, tol = .50){
  cat("Finding the optimal design by constructing the treatment design matrix from scratch:\n")

Rep = nrow(X.trt)/ncol(X.trt)
nTrt = ncol(X.trt)
n = nrow(X.trt)
X.trt = matrix(0, ncol = nTrt, nrow= 0)



old.temp = list(X.trt = X.trt)
temp = list()


s.pt = 1

  pb <- txtProgressBar(min = 0, max = n, style = 3)
for( i in s.pt:n){
  setTxtProgressBar(pb, i)

  for(j in 1:ncol(X.trt)){
      temp.X.trt = old.temp$X.trt

      if(!is.list(temp.X.trt)){
        temp$X.trt[[j]] = addPt(temp.X.trt,  index = j)
      }else{
        temp$X.trt[[j]] =lapply(temp.X.trt, addPt, index = j)

      }

  }

    if(i != s.pt){

      newtemp = lapply(temp$X.trt, function(y) extractMS1(y, index = i))
      newtemp = lapply(newtemp, function(y) extractMS2(y, index = i))

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
   unlist(lapply(temp$X.trt, function(x) tr(t(x) %*% blk.proj[1:i, 1:i] %*% x)))

  #print(length(temp$trace.screen))

  if(length(temp$trace.screen) > samp){

    index = sample(samp)
    temp$trace.screen = temp$trace.screen[index[1:(samp*tol)]]
    temp$X.trt = temp$X.trt[index[1:(samp*tol)]]
  }

  old.temp = temp
  temp = list()

}
close(pb)

#M- criteria
old.temp$trace.screen =
   sapply(old.temp$X.trt, function(x) ave.eff(x, blk.proj = blk.proj, Rep = Rep))

    old.temp$X.trt =  old.temp$X.trt[which(old.temp$trace.screen == max(old.temp$trace.screen))]
    old.temp$trace.screen = old.temp$trace.screen[which(old.temp$trace.screen == max(old.temp$trace.screen))]

#S- criteria
old.temp$trace.screen =
   sapply(old.temp$X.trt, function(x) bal.eff(x, blk.proj = blk.proj, Rep = Rep))

    old.temp$X.trt =  old.temp$X.trt[which(old.temp$trace.screen == min(old.temp$trace.screen))]
    old.temp$trace.screen = old.temp$trace.screen[which(old.temp$trace.screen == min(old.temp$trace.screen))]



  return(old.temp$X.trt)
}

