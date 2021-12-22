setwd("C:/Documents and Settings/krug001/My Documents/My Research/Proteomics/Paper 1/Decomposition example")

#  Incidence matrix

N <- matrix(0, nrow = 6, ncol = 3)
N[1,1] <- N[2,2] <- N[3,2]  <- N[4,3] <- N[5,3] <- N[6,1] <- 1

N <- rbind(N,N)

#   Set up Ci matrices

J6 <- matrix(1, nrow = 6, ncol = 6)
J2 <- matrix(1, nrow = 2, ncol = 2)
I6 <- diag(6)
I2 <- diag(2)


C0 <- kronecker(J6,J2)/12
C1 <- (I6 - J6/6) %x% J2/2
C2 <- I6 %x% (I2 - J2/2)

#  Working matrices

NC0N <- t(N) %*% C0 %*% N
NC1N <- t(N) %*% C1 %*% N
NC2N <- t(N) %*% C2 %*% N

C0N <- t(N) %*% C0
C1N <- t(N) %*% C1
C2N <- t(N) %*% C2


H0 <- eigen(C0)$vectors
H1 <- eigen(C1)$vectors
H2 <- eigen(C2)$vectors

round(t(H0) %*% C0 %*% H0,1)
round(t(H1) %*% C1 %*% H1,1)
round(t(H2) %*% C2 %*% H2,1)

# shows r_0 = 1
#  r_1 = 5
#  r_2 = 6

# Set up the matrices H.01, H.11, H.21 as at foot of p166

#  H.01  vector of 1's
#  H.11  5 contrasts between block totals
#  H.21  6 within block contrasts

H.01 <- t(H0)[1,]
H.11 <- t(H1)[1:5,]
H.21 <- t(H2)[1:6,]

HC0N <- H.01 %*% C0 %*% N
HC1N <- H.11 %*% C1 %*% N
HC2N <- H.21 %*% C2 %*% N

NC0HN <- t(HC0N) %*% HC0N
NC1HN <- t(HC1N) %*% HC1N
NC2HN <- t(HC2N) %*% HC2N

#   These NCiHN matrices are the same as NCiN, confirm numerically algebra top of page 167

#   now have to do eigen decompositions of the NCN matrices

C0.e <- eigen(NC0N)$vectors
Q.01 <- C0.e[,1] %*% t(C0.e[,1])

#   Q.01 is J_3/3, and lambda.01 is 4

C1.e <- eigen(NC1N)$vectors
Q.11 <- C1.e[,1:2] %*% t(C1.e[,1:2])

#   Q.11 is I_3 - J_3/3 with 2 eigenvalues of 1

C2.e <- eigen(NC2N)$vectors
Q.22 <- C2.e[,1:2] %*% t(C2.e[,1:2])

#   Q.11 is I_3 - J_3/3 with 2 eigenvalues of 3

# For C0 matrix in expression 2.4, estiamte is grand mean

Q.01 %*% C0N /4

# For C1 matrix in expression 2.4, estimates given by contrasts between 4 block totals
# containing the treatment and 2 block totals not containing it

Q.11 %*% C1N

# For C2 matrix in expression 2.4, estimates given by contrasts within blocks containing the treatment

Q.22 %*% C2N

# This completes Section 2

# Start Section 3

# Some useful functions

projMat <- function(X) X %*% solve(t(X) %*% X) %*% t(X)
J <- function(n) matrix(1, nrow=n, ncol=n)
K <- function(n) J(n)/n
tr <- function(X) sum(diag(X))

nTrt <- 3
T0 <- K(nTrt)
T1 <- diag(nTrt) - T0
T <- list(T0,T1)
Z <- list(b=diag(6)%x%rep(1,2),w=diag(12))

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

Ab   <- lapply(T, function(x) InfMat(C=C1, N=N, T=x))
Ab.Inv <- lapply(T, function(x) invInfMat(C=C1, N=N, T=x))
Mb <- lapply(T, function(x) x %*% t(N) %*% C1 %*% Z$b)
lapply(lapply(Mb, function(m) t(m) %*% Ab.Inv[[1]] %*% m), function(y) round(tr(y),4))
lapply(lapply(Mb, function(m) t(m) %*% Ab.Inv[[2]] %*% m), function(y) round(tr(y),4))

Mw <- lapply(T, function(x) x %*% t(N) %*% C1 %*% Z$w)
lapply(lapply(Mw, function(m) t(m) %*% Ab.Inv[[1]] %*% m), function(y) round(tr(y),4))
lapply(lapply(Mw, function(m) t(m) %*% Ab.Inv[[2]] %*% m), function(y) round(tr(y),4))

Aw   <- lapply(T, function(x) InfMat(C=C2, N=N, T=x))
Aw.Inv <- lapply(T, function(x) invInfMat(C=C2, N=N, T=x))
Mb <- lapply(T, function(x) x %*% t(N) %*% C2 %*% Z$b)
lapply(lapply(Mb, function(m) t(m) %*% Ab.Inv[[1]] %*% m), function(y) round(tr(y),4))
lapply(lapply(Mb, function(m) t(m) %*% Ab.Inv[[2]] %*% m), function(y) round(tr(y),4))

Mw <- lapply(T, function(x) x %*% t(N) %*% C2 %*% Z$w)
lapply(lapply(Mw, function(m) t(m) %*% Ab.Inv[[1]] %*% m), function(y) round(tr(y),4))
lapply(lapply(Mw, function(m) t(m) %*% Ab.Inv[[2]] %*% m), function(y) round(tr(y),4))
