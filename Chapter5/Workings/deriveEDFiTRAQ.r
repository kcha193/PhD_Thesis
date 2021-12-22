#12/03/2012 4:25:34 p.m.

library(infoDecompuTE)
library(asreml)

################################################################################
#Phase 1 design - Completely randomised design
#Assigning four animals (A, B, C, D) to each of two treatment groups (a, b)
#First experiments contain the treatment effects and biological effects
#(variation between different animals).


phase1DesignEX1 <- local({
  Ani = as.factor(LETTERS[1:4])
  Trt = as.factor(c("Con", "Dis")[rep(1:2, each = 2)])
  data.frame(Ani,Trt)
})
phase1DesignEX1

getVCs.onePhase(phase2designEX4, blk.str = "Ani", trt.str = "Trt")

################################################################################
#Phase 2 design - MudPIT-iTRAQ experiments - 4 techincal repicates
phase2designEX4 <- local({
  Run = as.factor(rep(1:4, each = 4))
  Ani = as.factor(LETTERS[c(1,2,3,4,
                            2,3,4,1,
                            3,4,1,2,
                            4,1,2,3)])
  Tag = as.factor(c(114,115,116,117)[rep(1:4, 4)])
  Trt = phase1DesignEX1$Trt[match(Ani, phase1DesignEX1$Ani)]
  data.frame(Run, Ani, Tag, Trt)
})
phase2designEX4





run.eff = rnorm(4, mean = 0, sd = 10)
ani.eff = rnorm(4, mean = 0, sd = 5)
trt.eff = c(1, 2)
tag.eff = c(0,0,0,0)
res.eff = rnorm(16, mean = 0, sd = 0.5)

y = with(phase2designEX4, run.eff[Run] + ani.eff[Ani] + tag.eff[Tag] + trt.eff[Trt]) + res.eff


asreml.fit = asreml(y ~ Trt + Tag, random = ~ Run + Ani, data = phase2designEX4)

summary(asreml.fit)
svc.asreml(asreml.fit)

wald(asreml.fit)

summary(aov(y ~ Tag + Trt + Error(Run + Ani), data = phase2designEX4))


getVCs.twoPhase(phase2designEX4, blk.str1 = "Ani", blk.str2 = "Run",
trt.str = "Trt + Tag")

getVCs.onePhase(phase2designEX4, blk.str = "Run + Ani",
trt.str = "Trt + Tag")


getVCs.twoPhase(phase2designEX4[-c(1:4),], blk.str1 = "Ani", blk.str2 = "Run",
trt.str = "Trt + Tag")


getVCs.twoPhase(phase2designEX4[-c(1:8),], blk.str1 = "Ani", blk.str2 = "Run",
trt.str = "Trt + Tag")


################################################################################
#Phase 1 design - Completely randomised design
#Assigning six animals (A, B, C, D, E, F) to each of two treatment groups (a, b)
#First experiments contain the treatment effects and biological effects
#(variation between different animals).


phase1DesignEX3 <- local({
  Ani = as.factor(LETTERS[1:6])
  Trt = as.factor(c("Con", "Dis")[rep(1:2, each = 3)])
  data.frame(Ani,Trt)
})
phase1DesignEX3

getVCs.onePhase(phase1DesignEX3, blk.str = "Ani", trt.str = "Trt")

################################################################################
#Phase 2 design - MudPIT-iTRAQ experiments - 4 techincal repicates
phase2designEX6 <- local({
  Run = as.factor(rep(1:3, each = 4))
  Ani = as.factor(LETTERS[c(1,2,4,5,
                            5,6,2,3,
                            3,1,6,4)])
  Tag = as.factor(c(114,115,116,117)[rep(1:4, 3)])
  Trt = phase1DesignEX3$Trt[match(Ani, phase1DesignEX3$Ani)]
  data.frame(Run, Ani, Tag, Trt)
})
phase2designEX6

getVCs.twoPhase(phase2designEX6, blk.str1 = "Ani", blk.str2 = "Run",
trt.str = "Trt + Tag")

################################################################################
#Phase 1 design - Completely randomised design
#Assigning eight animals (A, B, C, D, E, F, G, H) to each of two treatment groups (a, b)
#First experiments contain the treatment effects and biological effects (variation between different animals).

phase1DesignEX2 <- local({
  Ani = as.factor(LETTERS[1:8])
  Trt = as.factor(c("Con", "Dis")[c(rep(1:2, each = 2), rep(2:1, each = 2))])
  data.frame(Ani, Trt)
})
phase1DesignEX2

getVCs.onePhase(phase1DesignEX2, blk.str = "Ani", trt.str = "Trt")



################################################################################
#Phase 2 design - MudPIT-iTRAQ experiments -2 techincal repicates
phase2designEX5 <- local({
  Run = as.factor(rep(1:4, each = 4))
  Ani = as.factor(LETTERS[c(1,2,3,4,
                            3,4,1,2,
                            5,6,7,8,
                            7,8,5,6)])
  Tag = as.factor(c(114,115,116,117)[rep(1:4, 4)])
  Trt = phase1DesignEX2$Trt[match(Ani, phase1DesignEX2$Ani)]
  Rep = as.factor(1:16)
  data.frame(Run, Ani, Tag, Trt, Rep)
})
phase2designEX5


getVCs.twoPhase(phase2designEX5, blk.str1 = "Ani", blk.str2 = "Run",
trt.str = "Trt + Tag")

getVCs.twoPhase(phase2designEX5[-c(1),], blk.str1 = "Ani", blk.str2 = "Run",
trt.str = "Trt + Tag")


y = rnorm(16)



summary(aov(y ~ Tag + Trt + Error(Run) + Error(Ani), data = phase2designEX5))

