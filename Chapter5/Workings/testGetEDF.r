
################################################################################
#\\sigma(run)^2/\\sigma^2 = 0"
#

EDF = numeric(3)
r.EDF = numeric(3)
VC = numeric(3)

x = c(0.00001, 0.0001, 0.001, 0.01, 0.1, 1, 10, 100, 1000)
n = length(x)

total <- n
# create progress bar
pb <- txtProgressBar(min = 0, max = total, style = 3)

for( i in 1:n){
  setTxtProgressBar(pb, i)
  for( j in 1:1000){
  gamma.run = 0
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

  #G.mat = getGMat(MS.likelihood = MS.likelihood, trt.str = "Trt + Tag",
  #blk.str1 = "Ani", blk.str2 = "Run", design = design)

  temp = getVcEDF(MS.likelihood = MS.likelihood, G.mat = G.mat, real.VC = real.VC)

  EDF = rbind(EDF, temp$S$E)
  r.EDF = rbind(r.EDF, temp$S$"r.EDF")
  VC =  rbind(VC, t(temp$V))
  }
}

EDF = EDF[-1,]
r.EDF = r.EDF[-1,]

apply(EDF[1:1000,], 2, mean)
apply(EDF[1001:2000,], 2, mean)
apply(EDF[2001:3000,], 2, mean)
apply(EDF[3001:4000,], 2, mean)
apply(EDF[4001:5000,], 2, mean)
apply(EDF[5001:6000,], 2, mean)
apply(EDF[6001:7000,], 2, mean)
apply(EDF[7001:8000,], 2, mean)
apply(EDF[8001:9000,], 2, mean)

apply(r.EDF[1:1000,], 2, mean)
apply(r.EDF[1001:2000,], 2, mean)
apply(r.EDF[2001:3000,], 2, mean)
apply(r.EDF[3001:4000,], 2, mean)
apply(r.EDF[4001:5000,], 2, mean)
apply(r.EDF[5001:6000,], 2, mean)
apply(r.EDF[6001:7000,], 2, mean)
apply(r.EDF[7001:8000,], 2, mean)
apply(r.EDF[8001:9000,], 2, mean)

#Concentrate the EDF of Between Animal
plot(1:length(x), EDF[-1,2],  type = "l", xaxt = "n", ylab = "EDF",
xlab = "\\simga_(animal)^2/\\sigma^2",
main = "EDF versus \\simga_(animal)^2/\\sigma^2 (\\simga_(run)^2/\\sigma^2 is fixed)" )
axis(1, 1:length(x), as.character(x))


################################################################################
#\\sigma(animal)^2/\\sigma^2 = 0.25"
#

EDF = numeric(3)
r.EDF = numeric(3)
VC = numeric(3)

x = c(0.00001, 0.0001, 0.001, 0.01, 0.1, 1, 10, 100, 1000)
#x = seq(0.001, 100, length = 1000)
n = length(x)

total <- n
# create progress bar
pb <- txtProgressBar(min = 0, max = total, style = 3)

for( i in 1:n){
  setTxtProgressBar(pb, i)

  gamma.run = x[i]
  gamma.ani =   0.0000001

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

  #G.mat = getGMat(MS.likelihood = MS.likelihood, trt.str = "Trt + Tag",
  #blk.str1 = "Ani", blk.str2 = "Run", design = design)

  temp = getVcEDF(MS.likelihood = MS.likelihood, G.mat = G.mat, real.VC = real.VC)

  EDF = rbind(EDF, temp$S$E)
  r.EDF = rbind(r.EDF, temp$S$"r.EDF")
  VC =  rbind(VC, t(temp$V))

}

#Concentrate the EDF of Between Animal
plot(1:length(x), EDF[-1,2],  type = "l", xaxt = "n", ylab = "EDF",
xlab = "\\simga_(run)^2/\\sigma^2",
main = "EDF versus \\simga_(run)^2/\\sigma^2 (\\sigma(animal)^2/\\sigma^2 is fixed)", col = "blue" )
axis(1, 1:length(x), as.character(x))

################################################################################
#\\sigma(run)^2/\\sigma^2 = 1"
#

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

  gamma.run = 0
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

  #G.mat = getGMat(MS.likelihood = MS.likelihood, trt.str = "Trt + Tag",
  #blk.str1 = "Ani", blk.str2 = "Run", design = design)

  temp = getVcEDF(MS.likelihood = MS.likelihood, G.mat = G.mat, real.VC = real.VC)

  EDF = rbind(EDF, temp$S$E)
  r.EDF = rbind(r.EDF, temp$S$"r.EDF")
  VC =  rbind(VC, t(temp$V))

}

#Concentrate the EDF of Between Animal
lines(1:length(x), r.EDF[-1,2],  type = "l", xaxt = "n", ylab = "EDF",
xlab = "\\simga_(animal)^2/\\sigma^2",
main = "EDF when \\sigma(run)^2/\\sigma^2 = 0", col = "red" )

axis(1, 1:length(x), as.character(x))


################################################################################
#\\sigma(run)^2/\\sigma^2 = 4"
#

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

  gamma.run = 4
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

  #G.mat = getGMat(MS.likelihood = MS.likelihood, trt.str = "Trt + Tag",
  #blk.str1 = "Ani", blk.str2 = "Run", design = design)

  temp = getVcEDF(MS.likelihood = MS.likelihood, G.mat = G.mat, real.VC = real.VC)

  EDF = rbind(EDF, temp$S$E)
  r.EDF = rbind(r.EDF, temp$S$"r.EDF")
  VC =  rbind(VC, t(temp$V))

}

#Concentrate the EDF of Between Animal
plot(1:length(x), r.EDF[-1,2],  type = "l", xaxt = "n", ylab = "EDF",
xlab = "\\simga_(animal)^2/\\sigma^2",
main = "EDF when \\sigma(run)^2/\\sigma^2 = 0" )

axis(1, 1:length(x), as.character(x))


################################################################################
#\\sigma(run)^2/\\sigma^2 = 100"
#

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

  gamma.run = 100
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

  #G.mat = getGMat(MS.likelihood = MS.likelihood, trt.str = "Trt + Tag",
  #blk.str1 = "Ani", blk.str2 = "Run", design = design)

  temp = getVcEDF(MS.likelihood = MS.likelihood, G.mat = G.mat, real.VC = real.VC)

  EDF = rbind(EDF, temp$S$E)
  r.EDF = rbind(r.EDF, temp$S$"r.EDF")
  VC =  rbind(VC, t(temp$V))

}

#Concentrate the EDF of Between Animal
plot(1:length(x), r.EDF[-1,2],  type = "l", xaxt = "n", ylab = "EDF",
xlab = "\\simga_(animal)^2/\\sigma^2",
main = "EDF when \\sigma(run)^2/\\sigma^2 = 0" )

axis(1, 1:length(x), as.character(x))

################################################################################
#\\sigma(run)^2/\\sigma^2 = 100"
#

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

  gamma.run = 100
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

  #G.mat = getGMat(MS.likelihood = MS.likelihood, trt.str = "Trt + Tag",
  #blk.str1 = "Ani", blk.str2 = "Run", design = design)

  temp = getVcEDF(MS.likelihood = MS.likelihood, G.mat = G.mat, real.VC = real.VC)

  EDF = rbind(EDF, temp$S$E)
  r.EDF = rbind(r.EDF, temp$S$"r.EDF")
  VC =  rbind(VC, t(temp$V))

}

#Concentrate the EDF of Between Animal
plot(1:length(x), r.EDF[-1,2],  type = "l", xaxt = "n", ylab = "EDF",
xlab = "\\simga_(animal)^2/\\sigma^2",
main = "EDF when \\sigma(run)^2/\\sigma^2 = 0" )

axis(1, 1:length(x), as.character(x))
