
#######################################################################################

#Function that swap any two random row indexes
swap<-function(ind){

   swap <- sample(ind, 2)   # index of 1 to N2 grid cells

   ind[swap[1]] =  ind[swap[2]]
   ind[swap[2]] =  ind[swap[1]]

   cat("crit = ",obj.fun(ind),"\n")

   return(ind)
 }


#Function to calculate the average efficiency factors
obj.fun=function(sq){
  newX = X.trt[sq,]

  info.mat =  t(newX) %*% blk.proj %*% newX

  e.va = eigen(info.mat)$va
  con.eff = e.va[-which(e.va<1e-7)]/Rep

  return(1/mean(1/con.eff))
}


set.seed(123)
init =  (1:nrow(X.trt))
#shuffle the row indexes and use it as the intial treatment design matrix
  
#  run SANN
res <- optim(init, obj.fun, swap, method="SANN",
  control =
    list(maxit=10000, temp = 1000, tmax = 10, trace = TRUE, REPORT = 100, fnscale = -1))

test(X.trt = X.trt[res$par,], blk.proj, Rep=Rep)

(y2 = y2/(1+0.25))

test(X.trt = X.trt[res$par,], blk.proj, Rep=Rep)

#run SANN with the new temperture generated from the paper
res1 <- optim(init, obj.fun, swap, method="SANN",
  control =
    list(maxit=20000, temp =  y2, tmax = 10, trace = TRUE, REPORT = 100, fnscale = -1))

test(X.trt = X.trt[res1$par,], blk.proj, Rep=Rep)

temp =  0.01080776
 tmax = 100000
 
temp / log((((1:10000)-1) %/% tmax)*tmax+ exp(1))

