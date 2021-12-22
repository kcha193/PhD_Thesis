setwd("C:/Documents and Settings/krug001/My Documents/My Research/Proteomics/Paper 1/Decomposition example")

projMat <- function(X) X %*% solve(t(X) %*% X) %*% t(X)
K <- function(n) matrix(1/n, nrow=n, ncol=n)

inv <- function(X)
  {ev <- eigen(t(X)%*%X)
   wh <-(ev$values>1.e-6)
   eval <- ev$values[wh]
   ev$vectors[,wh]%*%diag(1/eval,nrow=length(eval))%*%
         t(ev$vectors[,wh])}

# Suppose we have 16 subjects randomly assigned to one of four treatment groups (A,B,C,D)
# Let's set up the data for this CRD first

nSub <- 16
nTrt <- 4
nRep <- 4
nSeq <- 3

design <- local({

   Subject <- factor(rep(1:(nRep*nTrt), each=nSeq))
   Treatment <- factor(rep(LETTERS[1:4], each=nRep*nSeq))
   Sequence <- factor(rep(paste("Q", 1:3, sep=""), times=nRep*nTrt))
   
   data.frame(Subject, Treatment, Sequence)

})
design 

# Let's now set up the trt design contrast matrices to compute the EMS

X <- local({

   nRow <- nSub*nSeq
   Xt <- matrix(-1/nTrt, nrow=nRow, ncol=nTrt)
   ind <- with(design, match(Treatment, LETTERS[1:4]))
   Xt[cbind(1:nRow, ind)] <- Xt[cbind(1:nRow, ind)]+1
   dimnames(Xt) <- list(1:nRow,  LETTERS[1:4])

   Xs <- matrix(-1/nSeq, nrow=nSub*nSeq, ncol=nSeq)
   ind <- with(design, match(Sequence, paste("Q", 1:3, sep="")))
   Xs[cbind(1:nRow, ind)] <- Xs[cbind(1:nRow, ind)]+1
   dimnames(Xs) <- list(1:nRow,  paste("Q", 1:3, sep=""))

   T.S <- as.matrix(data.frame(apply(data.frame(Xt), 2, function(x)x*data.frame(Xs))))

   list(T=Xt, S=Xs, T.S=T.S)
})

Z <- local({

   nRow <- nSub*nSeq
   Zs <- matrix(0, nrow=nRow, ncol=nSub)
   Zs[cbind(1:nRow, with(design, match(Subject, 1:nSub)))] <- 1
   dimnames(Zs) <- list(1:nRow,  paste("S",1:nSub,sep=""))

   Zd <- diag(nRow)
   dimnames(Zd) <- list(1:nRow,  paste("Q",1:nRow,sep=""))

   list(S=Zs, D=Zd)
})

P <- list(bS=projMat(Z$S)-K(nSub*nSeq), wS = diag(nSub*nSeq)-projMat(Z$S))
V <- lapply(Z, function(x) x%*%t(x))
Q <- lapply(X, function(x) x %*% inv(x) %*% t(x))

Pbs <- lapply(X, function(x) P$bS %*% x)
Mbs <- lapply(Pbs, function(PX) PX %*% inv(PX) %*% t(PX))
lapply(Mbs, function(x) sum(diag(x)))

Pws <- lapply(X, function(x) P$wS %*% x)
Mws <- lapply(Pws, function(PX) PX %*% inv(PX) %*% t(PX))
lapply(Mws, function(x) sum(diag(x)))


Px <- list(Trt=projMat(X$Trt), Seq=projMat(X$Seq))
Px$Trt.Seq <- projMat(X$Trt.Seq)-Px$Trt-Px$Seq+K(nSub*nSeq)

# Project onto between and within subjects vector subspace

#-- Between subjects
 
# Coeffients of random VCs

sum(diag(P$bS %*% V$D))
sum(diag(P$bS %*% V$S))

# Coefficients of fixed VC

sum(diag(t(X$T) %*% P$bS %*% X$T))/(nTrt-1)
sum(diag(t(X$S) %*% P$bS %*% X$S))/(nSeq-1)
sum(diag(t(X$T.S) %*% P$bS %*% X$T.S))/((nTrt-1)*(nSeq-1))

#--- Project onto between treatments vector subspace

# Coeffients of random VCs

sum(diag(t(X$T) %*% P$bS %*% V$D)) # EQUALS sum(diag(Px$Trt %*% P$bS))
sum(diag(t(X$T) %*% P$bS %*% V$S))

# Coefficients of fixed VC

sum(diag(t(X$T) %*% P$bS %*% Q$T %*% P$bS %*% X$T))/(nTrt-1)
sum(diag(t(X$Seq) %*% P$bS %*% Px$Trt %*% P$bS %*% X$Seq))/(nSeq-1)
sum(diag(t(X$Trt.Seq) %*% P$bS %*% Px$Trt %*% P$bS %*% X$Trt.Seq))/((nTrt-1)*(nSeq-1))

#--- Project onto within treatments vector subspace

# Coeffients of random VCs

iPx <- diag(nrow(Px$Trt)) - projMat(X$Trt)
sum(diag(iPx %*% P$bS %*% V$D)) # EQUALS sum(diag(Px$Trt %*% P$bS))
sum(diag(iPx %*% P$bS %*% V$S))

# Coefficients of fixed VC

sum(diag(t(X$Trt) %*% P$bS %*% iPx %*% P$bS %*% X$Trt))/(nTrt-1)
sum(diag(t(X$Seq) %*% P$bS %*% iPx %*% P$bS %*% X$Seq))/(nSeq-1)
sum(diag(t(X$Trt.Seq) %*% P$bS %*% iPx %*% P$bS %*% X$Trt.Seq))/((nTrt-1)*(nSeq-1))

#-- Within subjects
 
# Coeffients of random VCs

sum(diag(P$wS %*% V$D))
sum(diag(P$wS %*% V$S))

# Coefficients of fixed VC

sum(diag(t(X$T) %*% P$wS %*% X$T))/(nTrt-1)
sum(diag(t(X$S) %*% P$wS %*% X$S))/(nSeq-1)
sum(diag(t(X$T.S) %*% P$wS %*% X$T.S))/((nTrt-1)*(nSeq-1))

#--- Project onto between sequences vector subspace

# Coeffients of random VCs

sum(diag(Px$Seq %*% P$wS %*% V$D)) # EQUALS sum(diag(Px$Seq %*% P$wS))
sum(diag(Px$Seq %*% P$wS %*% V$S))

# Coefficients of fixed VC

sum(diag(t(X$Trt) %*% P$wS %*% Px$Seq %*% P$wS %*% X$Trt))/(nTrt-1)
sum(diag(t(X$Seq) %*% P$wS %*% Px$Seq %*% P$wS %*% X$Seq))/(nSeq-1)
sum(diag(t(X$Trt.Seq) %*% P$wS %*% Px$Seq %*% P$wS %*% X$Trt.Seq))/((nTrt-1)*(nSeq-1))

#--- Project onto between treatments.sequences vector subspace

# Coeffients of random VCs

sum(diag(Px$Trt.Seq %*% P$wS %*% V$D)) # EQUALS sum(diag(Px$Trt.Seq %*% P$wS))
sum(diag(Px$Trt.Seq %*% P$wS %*% V$S))

# Coefficients of fixed VC

sum(diag(t(X$Trt) %*% P$wS %*% Px$Trt.Seq %*% P$wS %*% X$Trt))/(nTrt-1)
sum(diag(t(X$Seq) %*% P$wS %*% Px$Trt.Seq %*% P$wS %*% X$Seq))/(nSeq-1)
sum(diag(t(X$Trt.Seq) %*% P$wS %*% Px$Trt.Seq %*% P$wS %*% X$Trt.Seq))/((nTrt-1)*(nSeq-1))

# Within Subjects residual

P$R <- P$wS - Px$Seq - Px$Trt.Seq + K(nSub*nSeq)
sum(diag(P$R))