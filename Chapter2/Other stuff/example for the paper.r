library(devtools)

install_github("kcha193/infoDecompuTE", "kcha193")

#Example for the paper
library(infoDecompuTE)


package = function(packName){
  if(!require(packName, character.only = TRUE)) 
    install.packages(packName)
  
  require(packName, character.only = TRUE)
}


sourceDir <- function(path, trace = TRUE, ...) {
  package("MASS")
  package("formatR")
  package("ggplot2")
  package("reshape")
  for (nm in list.files(path, pattern = "\\.[RrSsQq]$")) {
    if(trace) cat(nm,":")
    source(file.path(path, nm), ...)
    if(trace) cat("\n")
  }
}

sourceDir(path = "C:/Users/kcha193/Dropbox/Packaging tools/packaging/infoDecompuTE/R")


#Phase 1 experiment
design1 <- local({ 
  Ani = as.factor(LETTERS[c(1,2,3,4,
                            5,6,7,8)])
   Trt = as.factor(c("healthy", "diseased")[c(1,2,1,2,
                            2,1,2,1)])
  data.frame(Ani, Trt)
})
design1

summaryAovOnePhase(design1, blk.str = "Ani", trt.str = "Trt")     

#Phase 2 experiment  
design2 <- local({ 
  Run = as.factor(rep(1:4, each = 4))
  Ani = as.factor(LETTERS[c(1,2,3,4,
                            5,6,7,8,
                            3,4,1,2,
                            7,8,5,6)])
  Sam = as.factor(as.numeric(duplicated(Ani)) + 1)
  Tag = as.factor(c(114,115,116,117)[rep(1:4, 4)])
  Trt = as.factor(c("healthy", "diseased")[c(1,2,1,2,
                            2,1,2,1,
                            1,2,1,2,
                            2,1,2,1)])
  data.frame(Run, Ani, Sam, Tag, Trt)
})
design2
                                  
summaryAovTwoPhase(design2, blk.str1 = "Ani", blk.str2 = "Run", 
trt.str = "Tag + Trt")                                            
   
#Add the sample into the Phase 1 block structure                                           
summaryAovTwoPhase(design2, blk.str1 = "Ani/Sam", blk.str2 = "Run", 
trt.str = "Tag + Trt")                                            

attr(terms(~Ani / Sam), "term.labels")

attr(terms(~Ani * Sam), "term.labels")


#Assuming there is crossing between the animals and samples 
summaryAovTwoPhase(design2, blk.str1 = "Ani*Sam", blk.str2 = "Run", 
trt.str = "Tag + Trt")                                            

 
#Set Artificial strata 
design2$AniSet = as.factor(as.numeric(design2$Run)%%2 + 1 )
design2

summaryAovTwoPhase(design2, blk.str1 =  "Ani/Sam", blk.str2 = "AniSet/Run", 
trt.str = "Tag + Trt", var.comp = c("Ani(Sam)", "Ani", "AniSet(Run)"))                                    

#Fit contrast Matrix                                   
TagA = rep(c(1,1,-1,-1),time = 4)                
TagB = rep(c(1,-1,1,-1),time = 4)                
TagC = TagA * TagB
Tag = cbind(TagA, TagB, TagC)
Tag

Trt = as.numeric(design2$Trt)-1.5
Trt

summaryAovTwoPhase(design2, blk.str1 = "Ani/Sam", blk.str2 = "Run", 
trt.str  = "Tag + Trt", trt.contr = list(Tag = Tag, Trt = Trt) )                                


#Fit contrast Vector    
Tag = list(TagA = TagA, TagB = TagB, TagC = TagC)
Tag

summaryAovTwoPhase(design2, blk.str1 = "Ani/Sam", blk.str2 = "Run", 
trt.str  = "Tag + Trt", trt.contr = list(Tag = Tag, Trt = Trt),
table.legend = TRUE)

Tag = list(B = TagB, AC = cbind(TagA, TagC))
Tag

summaryAovTwoPhase(design2, blk.str1 = "Ani/Sam", blk.str2 = "Run", 
trt.str  = "Tag + Trt", trt.contr = list(Tag = Tag, Trt = Trt),
table.legend = TRUE)

summaryAovTwoPhase(design2, blk.str1 = "Ani/Sam", blk.str2 = "Run", 
trt.str  = "Tag + Trt", trt.contr = list(Tag = Tag),
table.legend = TRUE)


summaryAovTwoPhase(design2, blk.str1 = "Ani/Sam", blk.str2 = "Run", 
trt.str  = "Tag * Trt", trt.contr = list(Tag = Tag),
table.legend = TRUE)

#Compute MS 
set.seed(527)
summaryAovTwoPhase(design2, blk.str1 = "Ani/Sam", blk.str2 = "Run", 
trt.str = "Tag + Trt", response = rnorm(16))$ANOVA                                            

#Generate Latex scripts
summaryAovTwoPhase(design2, blk.str1 = "Ani/Sam", blk.str2 = "Run", 
trt.str = "Tag + Trt", latex = TRUE, fixed.names = c("\\gamma", "\\tau"))  

#Generate Latex scripts with MS
set.seed(527)
summaryAovTwoPhase(design2, blk.str1 = "Ani/Sam", blk.str2 = "Run", 
trt.str = "Tag + Trt", response = rnorm(16), latex = TRUE, 
fixed.names = c("\\gamma", "\\tau") )                                             

                                         
