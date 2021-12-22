setwd("C:/Documents and Settings/krug001/My Documents/My Research/Proteomics/Paper 1/Decomposition example")

projMat <- function(X) X %*% solve(t(X) %*% X) %*% t(X)
K <- function(n) matrix(1/n, nrow=n, ncol=n)
tr <- function(X) sum(diag(X))
require(MASS)

# Suppose we have 16 subjects randomly assigned to one of four treatment groups (A,B,C,D)
# Let's set up the data for this CRD first

nSub <- 16
nTrt <- 4
nRep <- 4
nPep <- 3
v <- nTrt*nPep
nRow <- nSub*nPep

design <- local({

   Subject <- factor(rep(1:nSub, each=nPep))
   Treatment <- factor(rep(LETTERS[1:4], each=nRep*nPep))
   Sequence <- factor(rep(paste("Q", 1:3, sep=""), times=nSub))
   
   data.frame(Subject, Treatment, Sequence)

})
design 

# Let's now set up the trt and blk design matrices to compute the EMS

Z <- local({

   Zd <- diag(nRow)
   dimnames(Zd) <- list(1:nRow,  paste("Q",1:nRow,sep=""))

   Zs <- matrix(0, nrow=nRow, ncol=nSub)
   Zs[cbind(1:nRow, with(design, match(Subject, 1:nSub)))] <- 1
   dimnames(Zs) <- list(1:nRow,  paste("S",1:nSub,sep=""))

   list(D=Zd, S=Zs)
})
Z

N <- local({

   Xt <- matrix(0, nrow=nRow, ncol=nTrt)
   Xt[cbind(1:nRow, with(design, match(Treatment, LETTERS[1:4])))] <- 1
   dimnames(Xt) <- list(1:nRow,  LETTERS[1:4])

   Xs <- matrix(0, nrow=nSub*nPep, ncol=nPep)
   Xs[cbind(1:nRow, with(design, match(Sequence, paste("Q", 1:3, sep=""))))] <- 1
   dimnames(Xs) <- list(1:nRow,  paste("Q", 1:3, sep=""))

   T.S <- as.matrix(data.frame(apply(data.frame(Xt), 2, function(x)x*data.frame(Xs))))

   #list(Trt=Xt, Seq=Xs, Trt.Seq=T.S)
})
N

X <- local({

   t <- diag(nTrt) %x% rep(1,nPep)
   s <- rep(1,nTrt) %x% diag(nPep)
   dimnames(t) <- list(1:v,  LETTERS[1:4])
   dimnames(s) <- list(1:v,  paste("Q", 1:3, sep=""))
   ts <- as.matrix(data.frame(apply(data.frame(t), 2, function(x)x*data.frame(s))))

   list(t=t, s=s, ts=ts)
})
X   

P <- list(bS=projMat(Z$S)-K(nSub*nPep), wS = diag(nSub*nPep)-projMat(Z$S))
V <- lapply(Z, function(x) x%*%t(x))
T <- list(t=projMat(X$t)-K(nTrt*nPep), s=projMat(X$s)-K(nTrt*nPep), 
          ts=projMat(X$ts)-projMat(X$t)-projMat(X$s)+K(nTrt*nPep))
for(i in 1:length(T)) dimnames(T[[i]])<- list(dimnames(X$ts)[[2]],dimnames(X$ts)[[2]])

# Project onto between and within subjects vector subspace

### Between Subjects
  # Total SS

  # Random effects
  lapply(lapply(V, function(x) P$bS %*% x), tr)

  # Fixed effects
  lapply(T, function(x) round(tr(x %*% P$bS %*% x), 3))

InfMat <- function(C,N,T){

   ei <- eigen(t(T) %*% t(N) %*% C %*% N %*% T)
   nn <- length(ei$values)
   L <- matrix(0, nrow=nn, ncol=nn)
   for(i in 1:(nn)) {
      if(ei$values[i]<1e-6)next
      L <- L + ei$values[i]*ei$vectors[,i]%*%t(ei$vectors[,i])
   }
   return(L)
}

invInfMat <- function(C,N,T){

   ei <- eigen(t(T) %*% t(N) %*% C %*% N %*% T)
   nn <- length(ei$values)
   L <- matrix(0, nrow=nn, ncol=nn)
   for(i in 1:(nn)) {
      if(ei$values[i]<1e-6)next
      L <- L + (1/ei$values[i])*ei$vectors[,i]%*%t(ei$vectors[,i])
   }
   return(L)
}

Ab   <- lapply(T, function(x) InfMat(C=P$bS, N=N, T=x))
Ainv <- lapply(T, function(x) invInfMat(C=P$bS, N=N, T=x))
MbS <- lapply(T, function(x) x %*% t(N) %*% P$bS %*% Z$S)
lapply(lapply(MbS, function(m) t(m) %*% Ainv[[1]] %*% m), function(y) round(tr(y),4))
lapply(lapply(MbS, function(m) t(m) %*% Ainv[[2]] %*% m), function(y) round(tr(y),4))
lapply(lapply(MbS, function(m) t(m) %*% Ainv[[3]] %*% m), function(y) round(tr(y),4))

MbD <- lapply(T, function(x) x %*% t(N) %*% P$bS %*% Z$D)
lapply(lapply(MbD, function(m) t(m) %*% Ainv[[1]] %*% m), function(y) round(tr(y),4))
lapply(lapply(MbD, function(m) t(m) %*% Ainv[[2]] %*% m), function(y) round(tr(y),4))
lapply(lapply(MbD, function(m) t(m) %*% Ainv[[3]] %*% m), function(y) round(tr(y),4))

Aw   <- lapply(T, function(x) InfMat(C=P$wS, N=N, T=x))
Ainv <- lapply(T, function(x) invInfMat(C=P$wS, N=N, T=x))
MwS <- lapply(T, function(x) x %*% t(N) %*% P$wS %*% Z$S)
lapply(lapply(MwS, function(m) t(m) %*% Ainv[[1]] %*% m), function(y) round(tr(y),4))
lapply(lapply(MwS, function(m) t(m) %*% Ainv[[2]] %*% m), function(y) round(tr(y),4))
lapply(lapply(MwS, function(m) t(m) %*% Ainv[[3]] %*% m), function(y) round(tr(y),4))

MwD <- lapply(T, function(x) x %*% t(N) %*% P$wS %*% Z$D)
lapply(lapply(MwD, function(m) t(m) %*% Ainv[[1]] %*% m), function(y) round(tr(y),4))
lapply(lapply(MwD, function(m) t(m) %*% Ainv[[2]] %*% m), function(y) round(tr(y),4))
lapply(lapply(MwD, function(m) t(m) %*% Ainv[[3]] %*% m), function(y) round(tr(y),4))


# Adding run effects

nRun <- 4
nSub <- 16
nTrt <- 4
nRep <- 4
nPep <- 3
v <- nTrt*nPep
nRow <- nSub*nPep

design <- local({

   Run <- factor(rep(1:nRun, each=nTrt*nPep))
   Subject <- factor(rep(1:nSub, each=nPep))
   Treatment <- factor(rep(LETTERS[1:4], each=nPep, times=nRun))
   Sequence <- factor(rep(paste("Q", 1:3, sep=""), times=nSub))
   
   data.frame(Run, Subject, Treatment, Sequence)

})
design 

# Let's now set up the trt and blk design matrices to compute the EMS

Z <- local({

   Zd <- diag(nRow)
   dimnames(Zd) <- list(1:nRow,  paste("Q",1:nRow,sep=""))

   Zs <- matrix(0, nrow=nRow, ncol=nSub)
   Zs[cbind(1:nRow, with(design, match(Subject, 1:nSub)))] <- 1
   dimnames(Zs) <- list(1:nRow,  paste("S",1:nSub,sep=""))

   Zr <- matrix(0, nrow=nRow, ncol=nRun)
   Zr[cbind(1:nRow, with(design, match(Run, 1:nRun)))] <- 1
   dimnames(Zr) <- list(1:nRow,  paste("R",1:nRun,sep=""))

   list(D=Zd, S=Zs, R=Zr)
})
Z

N <- local({

   Xt <- matrix(0, nrow=nRow, ncol=nTrt)
   Xt[cbind(1:nRow, with(design, match(Treatment, LETTERS[1:4])))] <- 1
   dimnames(Xt) <- list(1:nRow,  LETTERS[1:4])

   Xs <- matrix(0, nrow=nSub*nPep, ncol=nPep)
   Xs[cbind(1:nRow, with(design, match(Sequence, paste("Q", 1:3, sep=""))))] <- 1
   dimnames(Xs) <- list(1:nRow,  paste("Q", 1:3, sep=""))

   T.S <- as.matrix(data.frame(apply(data.frame(Xt), 2, function(x)x*data.frame(Xs))))

   #list(Trt=Xt, Seq=Xs, Trt.Seq=T.S)
})
N

X <- local({

   t <- diag(nTrt) %x% rep(1,nPep)
   s <- rep(1,nTrt) %x% diag(nPep)
   dimnames(t) <- list(1:v,  LETTERS[1:4])
   dimnames(s) <- list(1:v,  paste("Q", 1:3, sep=""))
   ts <- as.matrix(data.frame(apply(data.frame(t), 2, function(x)x*data.frame(s))))

   list(t=t, s=s, ts=ts)
})
X   

# Let's now set up the trt and blk design matrices to compute the EMS

P <- list(bR=projMat(Z$R)-K(nSub*nPep), bS=projMat(Z$S)-projMat(Z$R), wS = diag(nSub*nPep)-projMat(Z$S))
V <- lapply(Z, function(x) x%*%t(x))

invL <- function(C,N,T){

   ei <- eigen(T %*% t(N) %*% C %*% N %*% T)
   nn <- length(ei$values)
   L <- matrix(0, nrow=nn, ncol=nn)
   for(i in 1:(nn)) {
      if(ei$values[i]<1e-6)next
      L <- L + (1/ei$values[i])*ei$vectors[,i]%*%t(ei$vectors[,i])
   }
   return(L)
}
### Between Runs Expected Trt SS

# Fixed effects
AbR   <- lapply(T, function(x) InfMat(C=P$bR, N=N, T=x))

# Random effects

lapply(lapply(Z, function(x) P$bR %*% (x %*% t(x))), tr)
lapply(lapply(Z, function(x) P$bS %*% (x %*% t(x))), tr)
lapply(lapply(Z, function(x) P$wS %*% (x %*% t(x))), tr)

Ainv <- lapply(T, function(x) invL(C=P$bR, N=N, T=x))  # Between Runs info matrix for Treatments, Sequences, Trt x Seq
# Treatments
MbRt <- lapply(Z, function(x) T$t %*% t(N) %*% P$bR %*% x)
MbRs <- lapply(Z, function(x) T$s %*% t(N) %*% P$bR %*% x)
MbRts <- lapply(Z, function(x) T$ts %*% t(N) %*% P$bR %*% x)
lapply(lapply(MbRt, function(m) t(m) %*% Ainv[[1]] %*% m), function(y) round(tr(y),4)) # Treatments
lapply(lapply(MbRs, function(m) t(m) %*% Ainv[[2]] %*% m), function(y) round(tr(y),4)) # Sequences
lapply(lapply(MbRts, function(m) t(m) %*% Ainv[[3]] %*% m), function(y) round(tr(y),4)) # Trt x Seq

### Between Subjects within Runs Expected Trt SS

# Fixed effects
AbS   <- lapply(T, function(x) InfMat(C=P$bS, N=N, T=x))
# Random effects
Ainv <- lapply(T, function(x) invL(C=P$bS, N=N, T=x))  # Between Runs info matrix for Treatments, Sequences, Trt x Seq
# Treatments
MbSt <- lapply(Z, function(x) T$t %*% t(N) %*% P$bS %*% x)
MbSs <- lapply(Z, function(x) T$s %*% t(N) %*% P$bS %*% x)
MbSts <- lapply(Z, function(x) T$ts %*% t(N) %*% P$bS %*% x)
lapply(lapply(MbSt, function(m) t(m) %*% Ainv[[1]] %*% m), function(y) round(tr(y),4)) # Treatments
lapply(lapply(MbSs, function(m) t(m) %*% Ainv[[2]] %*% m), function(y) round(tr(y),4)) # Sequences
lapply(lapply(MbSts, function(m) t(m) %*% Ainv[[3]] %*% m), function(y) round(tr(y),4)) # Trt x Seq

### Between Peptides within Subjects Expected Trt SS

# Fixed effects
AwS   <- lapply(T, function(x) InfMat(C=P$wS, N=N, T=x))
# Random effects
Ainv <- lapply(T, function(x) invL(C=P$wS, N=N, T=x))  # Between Runs info matrix for Treatments, Sequences, Trt x Seq
# Treatments
MwSt <- lapply(Z, function(x) T$t %*% t(N) %*% P$wS %*% x)
MwSs <- lapply(Z, function(x) T$s %*% t(N) %*% P$wS %*% x)
MwSts <- lapply(Z, function(x) T$ts %*% t(N) %*% P$wS %*% x)
lapply(lapply(MwSt, function(m) t(m) %*% Ainv[[1]] %*% m), function(y) round(tr(y),4)) # Treatments
lapply(lapply(MwSs, function(m) t(m) %*% Ainv[[2]] %*% m), function(y) round(tr(y),4)) # Sequences
lapply(lapply(MwSts, function(m) t(m) %*% Ainv[[3]] %*% m), function(y) round(tr(y),4)) # Trt x Seq

#### Start old stuff

InfMat <- function(C,N){

   ei <- eigen(t(N) %*% C %*% N)
   nn <- length(ei$values)
   L <- matrix(0, nrow=nn, ncol=nn)
   for(i in 1:(nn)) {
      if(ei$values[i]<1e-6)next
      L <- L + ei$values[i]*ei$vectors[,i]%*%t(ei$vectors[,i])
   }
   return(L)
}

invInfMat <- function(C,N){

   ei <- eigen(t(N) %*% C %*% N)
   nn <- length(ei$values)
   L <- matrix(0, nrow=nn, ncol=nn)
   for(i in 1:(nn)) {
      if(ei$values[i]<1e-6)next
      L <- L + (1/ei$values[i])*ei$vectors[,i]%*%t(ei$vectors[,i])
   }
   return(L)
}

Ab   <- InfMat(C=P$bS, N=N)
Ainv <- invInfMat(C=P$bS, N=N)
MbS <- lapply(T, function(x) x %*% t(N) %*% P$bS %*% Z$S)
lapply(lapply(MbS, function(m) t(m) %*% Ainv %*% m), function(y) round(tr(y),4))

MbD <- lapply(T, function(x) x %*% t(N) %*% P$bS %*% Z$D)
lapply(lapply(MbD, function(m) t(m) %*% Ainv %*% m), function(y) round(tr(y),4))

Aw   <- InfMat(C=P$wS, N=N)
Ainv <- invInfMat(C=P$wS, N=N)
MwS <- lapply(T, function(x) x %*% t(N) %*% P$wS %*% Z$S)
lapply(lapply(MwS, function(m) t(m) %*% Ainv %*% m), function(y) round(tr(y),4))

MwD <- lapply(T, function(x) x %*% t(N) %*% P$wS %*% Z$D)
lapply(lapply(MwD, function(m) t(m) %*% Ainv %*% m), function(y) round(tr(y),4))

#### End old stuff

#-- Within

#- Total Fixed
   Aw <- t(X$Trt.Seq) %*% P$wS %*% X$Trt.Seq
   Aw
#- Total Random
   lapply(lapply(V, function(x) P$wS %*% x), tr)



Px <- list(Trt=projMat(X$Trt), Seq=projMat(X$Seq))
Px$Trt.Seq <- projMat(X$Trt.Seq)-Px$Trt-Px$Seq+K(nSub*nPep)

# Project onto between and within subjects vector subspace

#-- Between subjects
 
# Coeffients of random VCs

sum(diag(P$bS %*% V$D))
sum(diag(P$bS %*% V$S))

# Coefficients of fixed VC

sum(diag(t(X$Trt) %*% P$bS %*% X$Trt))/(nTrt-1)
sum(diag(t(X$Seq) %*% P$bS %*% X$Seq))/(nPep-1)
sum(diag(t(X$Trt.Seq) %*% P$bS %*% X$Trt.Seq))/((nTrt-1)*(nPep-1))

#--- Project onto between treatments vector subspace

# Coeffients of random VCs

sum(diag(Px$Trt %*% P$bS %*% V$D)) # EQUALS sum(diag(Px$Trt %*% P$bS))
sum(diag(Px$Trt %*% P$bS %*% V$S))

# Coefficients of fixed VC

sum(diag(t(X$Trt) %*% P$bS %*% Px$Trt %*% P$bS %*% X$Trt))/(nTrt-1)
sum(diag(t(X$Seq) %*% P$bS %*% Px$Trt %*% P$bS %*% X$Seq))/(nPep-1)
sum(diag(t(X$Trt.Seq) %*% P$bS %*% Px$Trt %*% P$bS %*% X$Trt.Seq))/((nTrt-1)*(nPep-1))

#--- Project onto within treatments vector subspace

# Coeffients of random VCs

iPx <- diag(nrow(Px$Trt)) - projMat(X$Trt)
sum(diag(iPx %*% P$bS %*% V$D)) # EQUALS sum(diag(Px$Trt %*% P$bS))
sum(diag(iPx %*% P$bS %*% V$S))

# Coefficients of fixed VC

sum(diag(t(X$Trt) %*% P$bS %*% iPx %*% P$bS %*% X$Trt))/(nTrt-1)
sum(diag(t(X$Seq) %*% P$bS %*% iPx %*% P$bS %*% X$Seq))/(nPep-1)
sum(diag(t(X$Trt.Seq) %*% P$bS %*% iPx %*% P$bS %*% X$Trt.Seq))/((nTrt-1)*(nPep-1))

#-- Within subjects
 
# Coeffients of random VCs

sum(diag(P$wS %*% V$D))
sum(diag(P$wS %*% V$S))

# Coefficients of fixed VC

sum(diag(t(X$Trt) %*% P$wS %*% X$Trt))/(nTrt-1)
sum(diag(t(X$Seq) %*% P$wS %*% X$Seq))/(nPep-1)
sum(diag(t(X$Trt.Seq) %*% P$wS %*% X$Trt.Seq))/((nTrt-1)*(nPep-1))

#--- Project onto between sequences vector subspace

# Coeffients of random VCs

sum(diag(Px$Seq %*% P$wS %*% V$D)) # EQUALS sum(diag(Px$Seq %*% P$wS))
sum(diag(Px$Seq %*% P$wS %*% V$S))

# Coefficients of fixed VC

sum(diag(t(X$Trt) %*% P$wS %*% Px$Seq %*% P$wS %*% X$Trt))/(nTrt-1)
sum(diag(t(X$Seq) %*% P$wS %*% Px$Seq %*% P$wS %*% X$Seq))/(nPep-1)
sum(diag(t(X$Trt.Seq) %*% P$wS %*% Px$Seq %*% P$wS %*% X$Trt.Seq))/((nTrt-1)*(nPep-1))

#--- Project onto between treatments.sequences vector subspace

# Coeffients of random VCs

sum(diag(Px$Trt.Seq %*% P$wS %*% V$D)) # EQUALS sum(diag(Px$Trt.Seq %*% P$wS))
sum(diag(Px$Trt.Seq %*% P$wS %*% V$S))

# Coefficients of fixed VC

sum(diag(t(X$Trt) %*% P$wS %*% Px$Trt.Seq %*% P$wS %*% X$Trt))/(nTrt-1)
sum(diag(t(X$Seq) %*% P$wS %*% Px$Trt.Seq %*% P$wS %*% X$Seq))/(nPep-1)
sum(diag(t(X$Trt.Seq) %*% P$wS %*% Px$Trt.Seq %*% P$wS %*% X$Trt.Seq))/((nTrt-1)*(nPep-1))

# Within Subjects residual

P$R <- P$wS - Px$Seq - Px$Trt.Seq + K(nSub*nPep)
sum(diag(P$R))