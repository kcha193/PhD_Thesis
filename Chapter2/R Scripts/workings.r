#23/02/2011 11:09:32 a.m.
#Testing
sourceDir <- function(path, trace = TRUE, ...) {
  library(MASS)
  library(lattice)
  for (nm in list.files(path, pattern = "\\.[RrSsQq]$")) {
    if(trace) cat(nm,":")
      source(file.path(path, nm), ...)
      if(trace) cat("\n")
  }
}



sourceDir(path = "C:/Users/kcha193/My Dropbox/Packaging tools/packaging/infoDecompuTE/R")

                      Packaging tools\packaging\infoDecompuTE\R
design2 <- local({
  Run = as.factor(rep(1:4, each = 4))
  Ani = as.factor(LETTERS[c(1,2,3,4,
                            5,6,7,8,
                            3,4,1,2,
                            7,8,5,6)])
  Sam = as.factor(rep(1:2, each = 8))
  Tag = as.factor(c(114,115,116,117)[rep(1:4, 4)])
  Trt = as.factor(letters[c(1,2,1,2,
                            2,1,2,1,
                            1,2,1,2,
                            2,1,2,1)])
  data.frame(Run, Ani, Sam, Tag, Trt)
})


summary.aov.twoPhase(design2, blk.str1 = "Ani/Sam", blk.str2 = "Run", 
 trt.str = "Trt + Tag", trt.contr = list(Trt = NA, Tag = list(Tag1 = Tag1, Tag2 = Tag2, Tag3 = Tag3)))   


summary.aov.twoPhase(design2[-c(13:15),], blk.str1 = "Ani/Sam", blk.str2 = "Run", 
 trt.str = "Trt + Tag", trt.contr = list(Trt = NA, Tag = list(Tag1 = Tag1, Tag2 = Tag2, Tag3 = Tag3)))   

eigen(t(T[[1]]) %*% t(N) %*% C %*% N %*% T[[1]])

eigen(rr %*% t(T[[1]]) %*% t(N) %*% C %*% N %*% T[[1]] %*% rr)

design.df = design2[-13,] 

 summary.aov.twoPhase(design2, blk.str1 = "Ani/Sam", blk.str2 = "Run", 
 trt.str = "Trt")

 summary.aov.twoPhase(design2[-c(13:15),], blk.str1 = "Ani/Sam", blk.str2 = "Run", 
 trt.str = "Trt")









design = data.frame(Blk = factor(rep(1:4, time = 3)),
                    Ani = factor(c(1,2,3,4,
                                   2,3,4,1,
                                   3,4,1,2)),
                    Trt = factor(LETTERS[c(1,2,3,1,
                                           2,3,1,2,
                                           3,1,2,3)]))

y = rnorm(12)

summary(aov(y~ Trt + Ani + Error(Blk), design))

summary(aov(y~ Trt +  Error(Blk + Ani), design))


  getVCs.onePhase(design, random.terms = "Blk/Ani", fixed.term = "Trt")

  summary.aov.twoPhase(design2,  blk.str1 = "Ani",  blk.str2 ="Run",  trt.str  = "Trt", response = rnorm(16))

  summary.aov.onePhase(design2,  blk.str = "Ani", trt.str  = "Trt + Tag", response = rnorm(16))


summary.aov.twoPhase(design2, blk.str1 = "Ani/Sam", blk.str2 = "Run", 
  trt.str = "Trt + Tag", trt.contr = list(Trt = NA, Tag = list(Tag1 = Tag1, Tag2 = Tag2, Tag3 = Tag3)), 
  contr.matrix = TRUE)   

  summary.aov.twoPhase(design2, blk.str1 = "Ani/Sam", blk.str2 = "Run", 
  trt.str = "Trt + Tag", trt.contr = list(Trt = NA, Tag = list(Tag1 = Tag1, Tag2 = Tag2, Tag3 = Tag3)), 
  contr.matrix = FALSE, table.legend = TRUE)   




#Testing Richard's and Kathy's design
#Design 1
                                             

design.df = data.frame(Set = factor(rep(1:7, each = 6)),
                    Blk = factor(rep(1:7, each = 6)),
                    Array = factor(rep(1:3, time = 7, each = 2)),
                    Dye = factor(c("Gre", "Red")),
                    Treat = factor(LETTERS[c(1,2,2,4,4,1,
                                             2,3,3,5,5,2,
                                             3,4,4,6,6,3,
                                             4,5,5,7,7,4,
                                             5,6,6,1,1,5,
                                             6,7,7,2,2,6,
                                             7,1,1,3,3,7)]),
                    Plant = factor(letters[c( c(1,2,2,3,3,1),
                                              c(1,2,2,3,3,1)+3,
                                              c(1,2,2,3,3,1)+3*2,
                                              c(1,2,2,3,3,1)+3*3,
                                              c(1,2,2,3,3,1)+3*4,
                                              c(1,2,2,3,3,1)+3*5,
                                              c(1,2,2,3,3,1)+3*6)]))

getVCs.twoPhase(design.df , random.terms1 = "Set/Array", random.terms2 ="Plant", 
fixed.term = "Treat + Dye", var.comp = c("Plant", "Set:Array", "Blk"))


VC$random
VC$fixed

getVCs.twoPhase(design.df, random.terms1 = "Plant", random.terms2 ="Set/Array", 
fixed.term = "Treat + Dye", var.comp = c("Plant", "Set:Array", "Blk"))

VC$random
VC$fixed

#Design 2: multiple dye-swap design
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
                       
design.df 

#Two Phase                       
getVCs.twoPhase(design.df , random.terms2 = "Set/Array", random.terms1 = "Plant", 
fixed.terms = "Treat + Dye", var.comp = c("Plant", "Set:Array"))


design.df = data.frame(Set = factor(rep(1:4, each = 4)),
                       Array = factor(rep(1:2, time = 4, each = 2)),
                       Dye = factor(c("Gre", "Red")),
                       Treat = factor(LETTERS[c(1,2,2,1,
                                                1,2,2,1,
                                                1,2,2,1,
                                                1,2,2,1)]),
                       Plant = factor(LETTERS[c(1,2,2,1,
                                                1,2,2,1,
                                                1,2,2,1,
                                                1,2,2,1)]))                       
                       

#One Phase                       
getVCs.onePhase(design.df , random.terms = "Set/(Array*Plant)", 
fixed.terms = "Treat + Dye", var.comp = c("Plant", "Set:Array"))

#Two Phase                       
getVCs.twoPhase(design.df , random.terms2 = "Set/Array", random.terms1 = "Plant", 
fixed.terms = "Treat + Dye", var.comp = c("Plant", "Set:Array"))



#Design 3: alternating loop design
design.df = data.frame(Set = factor(rep(1:2, each = 8)),
                       Array =interaction(factor(rep(1:2, each = 8)), 
                            factor(rep(1:4, time = 2, each = 2))),
                       Dye = factor(c("Gre", "Red")),
                       Treat = factor(LETTERS[c(1,2,1,2,
                                                1,2,1,2,
                                                2,1,2,1,
                                                2,1,2,1)]),
                       Plant = interaction(c(1,1,2,2,3,3,4,4,
                                             1,2,2,3,3,4,4,1),
                                 factor(LETTERS[c(1,2,1,2,
                                                1,2,1,2,
                                                2,1,2,1,
                                                2,1,2,1)])))

design.df

#Fitting Array first
getVCs.twoPhase(design.df , random.terms2 = "Set", 
random.terms1 = "Array + Plant", 
fixed.terms = "Treat + Dye", var.comp = c("Plant", "Set:Array"))

#Fitting Plant first
getVCs.twoPhase(design.df , random.terms2 = "Set", 
random.terms1 = "Plant + Array", 
fixed.terms = "Treat + Dye", var.comp = c("Plant", "Set:Array"))

#Design 3: alternating loop design (Break down the Array and Plant to 3 orthogonal contrasts using the 2-K contrasts)
design.df = data.frame(Set = factor(rep(1:2, each = 8)),
                       Array =interaction(factor(rep(1:2, each = 8)), 
                            factor(rep(1:4, time = 2, each = 2))),
                       Dye = factor(c("Gre", "Red")),
                       Treat = factor(LETTERS[c(1,2,1,2,
                                                1,2,1,2,
                                                2,1,2,1,
                                                2,1,2,1)]),
                       Plant = interaction(c(1,1,2,2,3,3,4,4,
                                             1,2,2,3,3,4,4,1),
                                 factor(LETTERS[c(1,2,1,2,
                                                1,2,1,2,
                                                2,1,2,1,
                                                2,1,2,1)])))

design.df$Array1 = factor(rep(c(2,2,1,1,-2,-2,-1,-1),2))
design.df$Array2 = factor(rep(c(2,2,-2,-2,1,1,-1,-1),2))
design.df$Array3 = factor(rep(c(2,2,1,1,-2,-2,-1,-1),2) *rep(c(2,2,-2,-2,1,1,-1,-1),2)) 

design.df$Plant1 = factor(c(2,2,1,1,-2,-2,-1,-1,2,1,1,-2,-2,-1,-1,2))
design.df$Plant2 = factor(c(2,2,-2,-2,1,1,-1,-1,2,-2,-2,1,1,-1,-1,2))
design.df$Plant3 = factor(c(2,2,1,1,-2,-2,-1,-1,2,1,1,-2,-2,-1,-1,2) * c(2,2,-2,-2,1,1,-1,-1,2,-2,-2,1,1,-1,-1,2))

#Fitting Array first
getVCs.twoPhase(design.df , random.terms2 = "Set", 
random.terms1 = "Array1 + Array2 + Array3 + Plant1 + Plant2 + Plant3", 
fixed.terms = "Treat + Dye", var.comp = c("Plant", "Set:Array"))

#Fitting Plant first
getVCs.twoPhase(design.df , random.terms2 = "Set", 
random.terms1 = "Plant1 + Plant2 + Plant3 + Array1 + Array2 + Array3", 
fixed.terms = "Treat + Dye", var.comp = c("Plant", "Set:Array"))















#Test2
design = data.frame(
Sets = factor(rep(1:6,each=10)),
ArrwS = factor(rep(1:5,times=12)),
Col = factor(rep(c("R","G"),each=5,times=6)),
CO2 = factor(rep(c("B","Y","R"),each=20)),
ParP = factor(rep(c(1,2,3,4,2,4,1,3,4,1,2,3),each=5)),
HST = factor(rep(c(1:5,2,3,4,5,1,1:5,3,4,5,1,2),times=3)))

 design[,5] =as.factor(paste("P", design[,5], sep = ""))
 design[,6] =as.factor(paste("H", design[,6], sep = ""))


VC = getVCs.twoPhase.complex(design, random.terms1 = "ParP", random.terms2 ="Sets/ArrwS", fixed.term = "Col + CO2 + HST", var.comp = c("ParP", "Sets:ArrwS"))



VC = getVCs.twoPhase.complex(design, random.terms1 = "ParP", random.terms2 ="Sets/ArrwS + Sets:Col", fixed.term = "Col + CO2 + HST", var.comp = c("ParP", "Sets:ArrwS"))

VC$random
VC$effFactor
VC$fixed
                                                                      


lCO2 <- (as.numeric(CO2)-2)   # Linear effect 
qCO2 <- (3*lCO2^2-2)         # Quadratic effect    

    
CO2m <- cbind(lCO2,qCO2)


   
aHST <- as.numeric(HST)    
HST12 <- ((aHST==2)-(aHST==1)) 
HST45 <- ((aHST==5)-(aHST==4)) 
HST3 <- (-1+5*(aHST==3))    
HSTo <- ((aHST==1|aHST==2)- (aHST==4|aHST==5))     
HSTm <- cbind(HST12,HST45,HST3,HSTo)
HCm <- cbind(lCO2*HSTm,qCO2*HSTm)

trt.contr = list(Col =ifelse(as.character(Col) =="R", 1, -1), CO2 = CO2m, HST = HSTm)  






T = as.list(rep(1,3))
names(T)  <- colnames(effectsMatrix)

T[[1]] = cCol %x%  mK(ncol(cHCm)) 
T[[2]] = mK(ncol(cCol)) %x%cHST %x% mK(3) 
T[[3]] = mK(ncol(cCol))  %x%cHCm 



design = data.frame(
Sets = factor(rep(1:6,each=10)),
ArrwS = factor(rep(1:5,times=12)),
Col = factor(rep(c("R","G"),each=5,times=6)),
CO2 = factor(rep(c("B","Y","R"),each=20)),
ParP = factor(rep(c(1,2,3,4,2,4,1,3,4,1,2,3),each=5)),
HST = factor(rep(c(1:5,2,3,4,5,1,1:5,3,4,5,1,2),times=3)))




HSTms = apply(HSTm, 1, sum)

HCms = apply(HCm, 1, sum)


anova(lm(ydes ~ Sets/ArrwS+Col*Sets+ParP+HSTm+HCm))

summary(aov(ydes~ Col + HCm + HSTm +  Error( Sets/ArrwS + ParP)))

summary(aov(ydes~ Col + HCms + HSTms +  Error( Sets/ArrwS + ParP)))

design.df = design
random.terms1 = "ParP"
random.terms2 ="Sets/ArrwS + Col*Sets"
fixed.terms = "Col + HST*CO2" 
var.comp = c("ParP", "Sets:ArrwS")
trt.contr = trt.contr

VC = getVCs.twoPhase.complex(design, , , , , )

VC$random
VC$effFactor
VC$fixed

                   DF ParP Sets:ArrwS
Between Sets                         
   Between ParP    3  5    2         
   Residual                          
      CO2          2  0    2         
Between Sets:ArrwS                   
   Residual                          
      HST          4  0    2         
      HST:CO2      8  0    2         
      Residual     12 0    2         
Between Col                          
   Between ParP                      
      Col          1  5/3  0         
Between Sets:Col                     
   Between ParP    3  85/9 0         
   Residual        2  0    0         
Between Within                       
   Residual                          
      HST          4  0    0         
      HST:CO2      8  0    0         
      Residual     12 0    0         

                            Col HST  CO2 HST:CO2 eff.Col eff.HST eff.CO2 eff.HST:CO2
Between Sets ParP                                                                   
Between Sets Residual                20                          1                  
Between Sets:ArrwS Residual     9/2      3/2             3/8             3/8        
Between Col ParP            30                   1                                  
Between Sets:Col ParP                                                               
Between Sets:Col Residual                                                           
Between Within Residual         15/2     5/2             5/8             5/8        


#Brien and Payne's design

design.df = design
random.terms1 ="(Row*(Squ/Col))/Hal"
random.terms2 ="((Oc/In/St)*Ju)/Pos"
fixed.terms = "Tre*Met"
var.comp = NA
trt.contr = NA


VC = summary.aov.twoPhase(design, blk.str1 = "(Row*(Squ/Col))/Hal", 
blk.str2 = "((Oc/In/St)*Ju)/Pos", trt.str = "Tre*Met", 
table.legend = TRUE)

                           DF  a  b  c  d  e   f   g h i  j  k  l  m  n  
Between Oc                                                               
   Between Squ             1   12 24 96 72 288 0   1 4 16 48 0  24 96 288
Between Oc:In              4   0  0  0  0  0   0   1 4 16 0  0  24 96 0  
Between Oc:In:St                                                         
   Between Squ:Col                                                       
      Tre                  3   4  8  0  24 0   0   1 4 0  0  0  24 0  0  
      Residual             3   4  8  0  24 0   0   1 4 0  0  0  24 0  0  
   Residual                12  0  0  0  0  0   0   1 4 0  0  0  24 0  0  
Between Ju                 5   0  0  0  0  0   0   1 4 16 48 96 0  0  0  
Between Oc:Ju              5   0  0  0  0  0   0   1 4 16 48 0  0  0  0  
Between Oc:In:Ju                                                         
   Between Row             2   12 24 96 0  0   192 1 4 16 0  0  0  0  0  
   Between Row:Squ         2   12 24 96 0  0   0   1 4 16 0  0  0  0  0  
   Residual                16  0  0  0  0  0   0   1 4 16 0  0  0  0  0  
Between Oc:In:St:Ju                                                      
   Between Squ:Col                                                       
      Tre                  3   8  16 0  48 0   0   1 4 0  0  0  0  0  0  
      Residual             3   8  16 0  48 0   0   1 4 0  0  0  0  0  0  
   Between Row:Squ:Col                                                   
      Tre                  3   12 24 0  0  0   0   1 4 0  0  0  0  0  0  
      Residual             9   12 24 0  0  0   0   1 4 0  0  0  0  0  0  
   Residual                72  0  0  0  0  0   0   1 4 0  0  0  0  0  0  
Between Oc:In:St:Ju:Pos                                                  
   Between Row:Squ:Col:Hal                                               
      Met                  1   12 0  0  0  0   0   1 0 0  0  0  0  0  0  
      Tre:Met              3   12 0  0  0  0   0   1 0 0  0  0  0  0  0  
      Residual             20  12 0  0  0  0   0   1 0 0  0  0  0  0  0  
   Residual                408 0  0  0  0  0   0   1 0 0  0  0  0  0  0  

> Legend
[1] "a=Row:Squ:Col:Hal, b=Row:Squ:Col, c=Row:Squ, d=Squ:Col, e=Squ, f=Row, g=Oc:In:St:Ju:Pos, h=Oc:In:St:Ju, i=Oc:In:Ju, j=Oc:Ju, k=Ju, l=Oc:In:St, m=Oc:In, n=Oc"

                                        Tre  Met Tre:Met eff.Tre eff.Met eff.Tre:Met
Between Oc Squ                                                                      
Between Oc:In Residual                                                              
Between Oc:In:St Squ:Col                16/3             1/27                       
Between Oc:In:St Residual                                                           
Between Ju Residual                                                                 
Between Oc:Ju Residual                                                              
Between Oc:In:Ju Row                                                                
Between Oc:In:Ju Row:Squ                                                            
Between Oc:In:Ju Residual                                                           
Between Oc:In:St:Ju Squ:Col             32/3             2/27                       
Between Oc:In:St:Ju Row:Squ:Col         128              8/9                        
Between Oc:In:St:Ju Residual                                                        
Between Oc:In:St:Ju:Pos Row:Squ:Col:Hal      288 72              1       1          
Between Oc:In:St:Ju:Pos Residual      

#iTRAQ design
design.df = design
random.terms1 ="Cage/Animal/Sample/Subsample"
random.terms2 ="Run + Tag"
fixed.terms = "Treat*Position"
var.comp = NA
trt.contr = NA
table.legend = FALSE


VCs = getVCs.twoPhase(design.df, random.terms1 = random.terms1, random.terms2 = random.terms2, fixed.terms = fixed.terms, var.comp = NA, trt.contr = NA, table.legend = FALSE)

 VCs$random
                                        DF Cage:Animal:Sample:Subsample Cage:Animal:Sample Cage:Animal Cage Run
Between Run                                                                                                    
   Between Cage:Animal                  1  1                            38/21              24/7        0    8  
   Between Cage:Animal:Sample           2  1                            5/3                0           0    8  
   Between Cage:Animal:Sample:Subsample 3  1                            0                  0           0    8  
Between Within                                                                                                 
   Between Cage                                                                                                
      Tag                               7  1                            13/7               25/7        7    0  
   Between Cage:Animal                                                                                         
      Tag                               7  1                            38/21              24/7        0    0  
   Between Cage:Animal:Sample                                                                                  
      Tag                               6  1                            22/15              0           0    0  
      Position                          1  1                            2                  0           0    0  
      Treat:Position                    1  1                            2                  0           0    0  
      Residual                          6  1                            79/45              0           0    0  
   Between Cage:Animal:Sample:Subsample                                                                        
      Tag                               4  1                            0                  0           0    0  
      Residual                          17 1                            0                  0           0    0  

                                         Tag Treat Position Treat:Position eff.Tag eff.Treat eff.Position eff.Treat:Position
Between.Run.Cage:Animal                                                                                                     
Between.Run.Cage:Animal:Sample                                                                                              
Between.Run.Cage:Animal:Sample:Subsample                                                                                    
Within.Cage                              1/7 28    4/7      2/7            1/49    1         1/49         1/49              
Within.Cage:Animal                       8/7       16/21    8/21           8/49              4/147        4/147             
Within.Cage:Animal:Sample                8/3       80/3     40/3           8/21              20/21        20/21             
Within.Cage:Animal:Sample:Subsample      6                                 6/7                                       
