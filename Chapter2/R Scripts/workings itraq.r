
design = read.csv("unRandDesign.csv", header = TRUE, colClasses = "factor")



getVCs.twoPhase(design, random.terms1 = "Cage/Animal/Sample/Subsample", 
                random.terms2 = "Run", fixed.terms = "Tag +Tag + Treat*Position", table.legend = TRUE )
                
design.df = design[,-1]

design.df = design
                
TagA = rep(c(1,1,1,1,-1,-1,-1,-1),each = 8)                
TagB = rep(c(1,1,-1,-1,1,1,-1,-1),each = 8)                
TagC = rep(c(1,-1,1,-1,1,-1,1,-1),each = 8)                
                
TagAB = TagA * TagB          
TagBC = TagB * TagC      
TagAC = TagA * TagC       
                
TagABC = TagAB * TagC                
                
TreatM =  as.numeric(design.df$Treat)               
PosM =  as.numeric(design.df$Position)               

TreatM = ifelse(TreatM == 2, -1, 1)
PosM = ifelse(PosM == 2, -1, 1)      

InterM = TreatM * PosM                

Tag = list(TagA = TagA, TagB = TagB, TagC = TagC, 
          TagAB= TagAB, TagBC = TagBC, TagAC = TagAC, TagABC = TagABC)



getVCs.twoPhase(design, random.terms1 = "Cage/Animal/Sample/Subsample", 
                random.terms2 = "Run", 
                fixed.terms = "Tag + Treat*Position", 
                trt.contr = list(Tag,  TreatM, PosM, InterM))

            
design.df = design[-c(which(design$Run=="1")),]

                
TagA = rep(c(1,1,1,1,-1,-1,-1,-1),each = 7)                
TagB = rep(c(1,1,-1,-1,1,1,-1,-1),each = 7)                
TagC = rep(c(1,-1,1,-1,1,-1,1,-1),each = 7)                
                
TagAB = TagA * TagB          
TagBC = TagB * TagC      
TagAC = TagA * TagC       
                
TagABC = TagAB * TagC                

Tag = list(TagA = TagA, TagB = TagB, TagC = TagC, 
          TagAB= TagAB, TagBC = TagBC, TagAC = TagAC, TagABC = TagABC)

TreatM =  as.numeric(design.df$Treat)               
PosM =  as.numeric(design.df$Position)               

TreatM = ifelse(TreatM == 2, -1, 1)
PosM = ifelse(PosM == 2, -1, 1)      

InterM = TreatM * PosM                

               
getVCs.twoPhase(design.df, random.terms1 = "Cage/Animal/Sample/Subsample", 
                random.terms2 = "Run", 
                fixed.terms = "Tag  +  Treat*Position", 
                trt.contr = list(Tag, TreatM, PosM, InterM), table.legend = TRUE )

 
getVCs.twoPhase(design.df, random.terms1 = "Cage/Animal/Sample/Subsample", 
                random.terms2 = "Run", 
                fixed.terms = "Treat*Position + Tag", 
                trt.contr = list(TreatM, PosM, InterM, Tag) )

getVCs.onePhase(design.df, random.terms = "Cage/Animal/Sample/Subsample", 
                fixed.terms = "Treat*Position + Tag", 
                trt.contr = list(TreatM, PosM, InterM, Tag))
                
                
                
                