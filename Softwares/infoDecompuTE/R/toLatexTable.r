##' Convert the R output to Latex Table
##' 
##' Print the Latex scripts on the screen for the user to output the table from
##' the Latex output.
##' 
##' Once the Latex script is generated, it requires the user to install and
##' load two Latex packages: \code{booktabs} and \code{bm} to compile the Latex
##' script.
##' 
##' @param ANOVA a matrix containing the coefficients of the variance
##' components in EMS of ANOVA table generated by
##' \code{\link{getCoefVC.onePhase}} or \code{\link{getCoefVC.twoPhase}}.
##' @param EF a matrix containing the coefficient of the fixed effects
##' components and the treatment average efficiency factors generated by
##' \code{\link{getFixedEF.onePhase}} or \code{\link{getFixedEF.onePhase}}
##' function.
##' @param fixed.names a vector of character allows the users to modify symbols
##' for the fixed effects.
##' @author Kevin Chang
##' @examples
##' 
##' design1 <- local({ 
##'   Ani = as.factor(LETTERS[c(1,2,3,4,
##'                             5,6,7,8)])
##'   Trt = as.factor(letters[c(1,1,1,1,
##'                             2,2,2,2)])
##'   data.frame(Ani, Trt)
##' })
##' 
##' blk.str <- "Ani"
##'     
##' rT <- terms(as.formula(paste("~", blk.str, sep = "")), keep.order = TRUE) 
##' blkTerm <- attr(rT,"term.labels")
##'      
##' Z <- makeBlkDesMat(design1, blkTerm)
##' 
##' 
##' trt.str = "Trt"              
##' fT <- terms(as.formula(paste("~", trt.str, sep = "")), keep.order = TRUE)  
##' 
##' trtTerm <- attr(fT, "term.labels")
##' effectsMatrix <- attr(fT, "factor")        
##' 
##' T <- makeContrMat(design1, trtTerm, effectsMatrix, contr.vec = NA)
##' 
##' N <- makeOverDesMat(design1, trtTerm)
##' 
##' Replist = getTrtRep(design1, trtTerm)   
##'  
##' Rep <- Replist$Rep
##' trt.Sca <- Replist$Sca
##'     
##' effFactors = lapply(makeOrthProjectors(Z), function(z) 
##'       getEffFactor(z, T, N, Rep, trt.Sca))
##' 
##' effFactors <- effFactors[sort(1:length(effFactors), decreasing=TRUE)]
##' 
##' v.mat <- getVMat.onePhase(Z.Phase1 = Z, design.df = design.df, var.comp = NA)
##'     
##' ANOVA <- getCoefVC.onePhase(Pb = effFactors, design.df = design1, v.mat = v.mat, 
##'     response = NA, table.legend = FALSE, decimal = FALSE, digits = 2)
##' 		
##' EF <- getFixedEF.onePhase(effFactors = effFactors, trt.Sca = trt.Sca,  T = T, 
##'   Rep = Rep, 
##' 	table.legend = FALSE, decimal = FALSE, digits = 2, list.sep = FALSE)
##' 
##' toLatexTable(ANOVA = ANOVA, EF = EF, fixed.names = c("\\tau"))
##' 
##' @export toLatexTable
toLatexTable <- function(ANOVA, EF, fixed.names) {
    
    
    matchRowNames <- rownames(ANOVA)
    random.ColNames <- colnames(ANOVA)
    fixed.ColNames <- colnames(EF)
    
    #browser()
    
    if (is.na(fixed.names[1])) 
        fixed.names <- c("\\tau", "\\gamma", "\\rho", "\\phi", "\\delta", "\\omega")
    
    fixed.ColNames <- gsub("\\.", "", fixed.ColNames)
    
    match.fixed.names <- unique(unlist(strsplit(fixed.ColNames[1:(length(fixed.ColNames)/2)], 
        "[[:punct:]]")))
    fixed.names <- fixed.names[1:length(match.fixed.names)]
    match.fixed.names <- c(match.fixed.names, "*", "(", ")")
    fixed.names <- c(fixed.names, "*", "(", ")")
    
    fixed.names <- sapply(strsplit(gsub("([[:punct:]])", "\\.\\1\\.", fixed.ColNames[1:(length(fixed.ColNames)/2)]), 
        "\\."), function(x) paste(fixed.names[match(x, match.fixed.names)], collapse = ""))
    
    fixed.names <- gsub("\\*", "", fixed.names)
    
    # browser()
    
    if (length(grep("MS", random.ColNames)) == 1) {
        ColNamesLen <- length(random.ColNames) - 1
    } else {
        ColNamesLen <- length(random.ColNames)
    }
    
    
    if (any(random.ColNames == "e")) {
        startColNames <- 3
    } else {
        startColNames <- 2
    }
    
    
    # check for interaction effects
    random.names <- sapply(strsplit(gsub("([[:punct:]])", "\\.\\1\\.",
                                         random.ColNames[startColNames:ColNamesLen]), 
        "\\."), function(x) paste(substr(x, 1, 1), collapse = ""))
    
    random.names <- tolower(random.names)
    
    random.names <- gsub("\\*", "", random.names)
    
    #random.names <- gsub("\\:", "", random.names)
    
    if (ColNamesLen > 2 && startColNames == 3) {
        random.names <- c("\\sigma^2", paste("\\sigma_{", random.names, "}^2", sep = ""))
    } else if (startColNames == 2) {
        random.names <- c(paste("\\sigma_{", random.names, "}^2", sep = ""))
    } else {
        random.names <- c("\\sigma^2")
    }
    
    ANOVA <- ifelse(ANOVA == "0", "", ANOVA)
    finalTable <- cbind(ANOVA, EF[match(matchRowNames, rownames(EF)), ])
    tempEF <- EF
    
    # avoid repeat in the fixed componenets <- check!!!!!
    for (i in 1:nrow(ANOVA)) {
        index <- which(matchRowNames[i] == rownames(tempEF))[1]
        
        if (is.na(index)) 
            next
        
        finalTable[i, (ncol(ANOVA) + 1):ncol(finalTable)] <- tempEF[index, ]
        if (nrow(tempEF) == 2) 
            next
        
        tempEF <- tempEF[-index, ]
        # rownames(EF)[grep(matchRowNames[i], rownames(EF))[1]] <- ''
    }
    
    finalTable <- ifelse(is.na(finalTable), "", finalTable)
    
    output <- "\\begin{table}[ht]\n\\centering\n \\caption{Theoretical ANOVA table}\n"
    
    if (length(grep("MS", random.ColNames)) == 1) {
        output <- c(output, paste("\\begin{tabular}[t]{lrll", paste(rep("l", length(fixed.names)), 
            collapse = ""), "} \n", sep = ""))
    } else {
        output <- c(output, paste("\\begin{tabular}[t]{lrl", paste(rep("l", length(fixed.names)), 
            collapse = ""), "} \n", sep = ""))
    }
    
    output <- c(output, "\\toprule \n")
    
    
    firstRow <- paste("\\multicolumn{1}{l}{\\textbf{Source of Variation}} & \\multicolumn{1}{l}{\\textbf{DF}} & \\multicolumn{1}{l}{\\textbf{EMS}}&", 
        sep = "")
    
    if (length(grep("MS", random.ColNames)) == 1) 
        firstRow <- rbind(firstRow, "\\multicolumn{1}{l}{\\textbf{MS}} &")
    
    firstRow <- rbind(firstRow, paste(paste("\\multicolumn{1}{l}{$\\bm{E_{", fixed.names, 
        sep = "", collapse = "}}$}&"), "}}$}\\\\ \n", sep = ""))
    
    # browser()
    
    output <- c(output, firstRow)
    output <- c(output, "\\midrule \n")
    for (i in 1:length(matchRowNames)) {
        
        # row names
        SV <- matchRowNames[i]
        SV <- gsub("   ", "\\\\quad ", SV)
        
        # DF
        DF <- paste("$", finalTable[i, 1], "$", sep = "")
        DF <- ifelse(DF == "$$", "", DF)
        
        # Random VC
        coef.VC <- finalTable[i, 2:(1 + length(random.names))]
        random.VC <- ""
        for (j in 1:length(coef.VC)) {
            if (coef.VC[j] == "") 
                next
            if (coef.VC[j] == 1) 
                coef.VC[j] <- ""
            
            random.VC <- paste(random.VC, coef.VC[j], random.names[j], sep = "")
            
            random.VC <- paste(random.VC, "+", sep = "")
        }
        
        # Fixed VC
        if (length(grep("MS", random.ColNames)) == 1) {
            coef.VC <- finalTable[i, (3 + length(random.names)):(2 + length(random.names) + 
                length(fixed.names))]
        } else {
            coef.VC <- finalTable[i, (2 + length(random.names)):(1 + length(random.names) + 
                length(fixed.names))]
        }
        
        fixed.VC <- ""
        for (j in 1:length(coef.VC)) {
            if (coef.VC[j] == "") 
                next
            
            if (coef.VC[j] == 1) 
                coef.VC[j] <- ""
            
			if( grepl(",", coef.VC[j])){
				fixed.VC <- paste(fixed.VC, "(", coef.VC[j], ")\\theta_{", fixed.names[j], "}", 
					sep = "")
			}else {			
				fixed.VC <- paste(fixed.VC, coef.VC[j], "\\theta_{", fixed.names[j], "}", 
					sep = "")
            }
			
            fixed.VC <- paste(fixed.VC, "+", sep = "")
        }
        
        total.VC <- ifelse(random.VC == "", "", paste("$", random.VC, "+", fixed.VC, 
            "$", sep = ""))
        total.VC <- ifelse(fixed.VC == "", paste("$", random.VC, "$", sep = ""), total.VC)
        total.VC <- gsub("\\+\\+", "\\+", total.VC)
        total.VC <- gsub("\\+\\$", "\\$", total.VC)
        total.VC <- gsub("\\$\\$", "", total.VC)
        
        if (length(grep("MS", random.ColNames)) == 1) {
            eff <- paste("$", finalTable[i, (3 + length(random.names) + length(fixed.names)):ncol(finalTable)], 
                "$", sep = "", collapse = " & ")
        } else {
            eff <- paste("$", finalTable[i, (2 + length(random.names) + length(fixed.names)):ncol(finalTable)], 
                "$", sep = "", collapse = " & ")
        }
        
        
        eff <- gsub("\\$\\$", "", eff)
        
        # browser()
        finalTable[i, 2:(1 + length(random.names))]
        
        if (length(grep("MS", random.ColNames)) == 1) {
            currentRow <- paste(SV, " & ", DF, " & ", total.VC, " & ", finalTable[i, 
                2 + length(random.names)], " &", eff, "\\\\", sep = "")
        } else {
            currentRow <- paste(SV, " & ", DF, " & ", total.VC, " &", eff, "\\\\", sep = "")
        }
        
        if ((grepl("Between", matchRowNames[i + 1]) || grepl("Within", matchRowNames[i + 
            1])) && !grepl("^Within", matchRowNames[i]) && !grepl("^Between", matchRowNames[i])) {
            currentRow <- paste(currentRow, "\\hline \n")
        } else {
            currentRow <- paste(currentRow, "\n")
        }
        
        output <- c(output, currentRow)
    }
    
    output <- c(output, "\\bottomrule \n \\end{tabular} \n \\label{tab:} \n\\end{table} \n")
    output <- c(output, "")
    
    return( cat(paste0(output, collapse = "")))
} 