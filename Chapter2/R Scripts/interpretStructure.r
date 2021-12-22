
   trtStruct <- "(Set*(Run/Animal))* (Genotype /Age)"
   trtStruct <- gsub(" ", "", trtStruct)
   leftBracket  <- gregexpr("\\(*", trtStruct)
   rightBracket <- gregexpr("\\)*", trtStruct, perl=TRUE)

trtStruct <- "(S*A)/(D*F)/(G*R)"
   structFormula <- as.formula(paste("y ~", trtStruct, sep="", collapse=""))
   attr(terms(structFormula, keep.order=TRUE), "term.labels")
 


   trtStructure <- "(Treatment /   Sequence) *Tag"
   trtStructure <- gsub(" ", "", trtStructure)
   strsplit(trtStructure, "\\*|/")
