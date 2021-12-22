sourceDir <- function(path, trace = TRUE, ...) {
    library(MASS)
    library(inline)
    library(compiler)
    library(formatR)
    library(Rcpp)
    library(RcppArmadillo)
    library(ggplot2)
    library(gridExtra)
    library(reshape)
    for (nm in list.files(path, pattern = "\\.[RrSsQq]$")) {
        if (trace)
            cat(nm, ":")
        source(file.path(path, nm), ...)
        if (trace)
            cat("\n")
    }
}

sourceDir(path = "C:/Users/kcha193/My Dropbox/Packaging tools/packaging/infoDecompuTE/R")
sourceDir(path = "C:/Users/kcha193/My Dropbox/Packaging tools/packaging/optimTE/R")


sourceDir(path = "C:/Users/Kevin/Dropbox/Packaging tools/packaging/infoDecompuTE/R")
sourceDir(path = "C:/Users/Kevin/Dropbox/Packaging tools/packaging/optimTE/R")


design.df = optCRD(nTrt = 6, bRep = 3, tRep = 2, nPlot = 4, iter = 1000)

phase2Design.df = design.df[, 1:2]
n = nrow(phase2Design.df)

blk.str1 = "Ani"
blk.str2 = "Run"
trt.str1 = "Trt"
trt.str2 = "Tag"

# Extract the fixed and random terms
rT2 = terms(as.formula(paste("~", blk.str2, sep = "")), keep.order = TRUE)  #random terms phase 2

rT2 = terms(as.formula(paste("~", paste(blk.str2, trt.str2, sep = " + "), sep = "")), keep.order = TRUE)  #random terms phase 2


blkTerm2 = attr(rT2, "term.labels")

Z2 = makeBlkDesMat(design.df, rev(blkTerm2))

# write('2. Defining the block structures of second Phase.', '')
Pb <- makeOrthProjectors(Z2)

blk.proj = Pb[[1]]  #OrthProjectors with mimial confounding


phase1Design.df = data.frame(Ani = factor(LETTERS[1:18]), Trt = factor(letters[1:6]))

summaryAovOnePhase(phase1Design.df, blk.str = blk.str1, trt.str = trt.str1)


tRep = nrow(phase2Design.df)/nrow(phase1Design.df)

if (!is.wholenumber(tRep)) {
    stop("The samples from Phase 1 experiment can not be fitted to the Phase 2 experiment!")
}


############################################################################# 
#Block treatament factors of interst from the Phase 1 experiment

fillIn = function(phase2Factor, phase1Design.df, design.df, count, ani.limit, trt.limit) {

    if (length(count) == 1) {
        return(names(count)[1])
    }


    for (i in 1:length(count)) {
        ani.check = rep(FALSE, length(ani.limit))
        trt.check = rep(FALSE, length(trt.limit))

        for (j in 1:length(ani.limit)) {
            ani.check[j] = sum(design.df[which(as.character(design.df[, colnames(design.df)[j]]) == as.character(phase2Factor[j])), colnames(phase1Design.df)[1]] ==
                names(count)[i]) < ani.limit[j]

            trt = as.character(phase1Design.df[, 2][which(phase1Design.df[, 1] == names(count)[i])])

            trt.check[j] = sum(design.df[which(as.character(design.df[, colnames(design.df)[j]]) == as.character(phase2Factor[j])), colnames(phase1Design.df)[2]] ==
                trt) < trt.limit[j]

        }

        if (all(ani.check) && all(trt.check)) {
            return(names(count)[i])
        } else if(as.character(phase2Factor[1]) ==  as.character(rev(design.df[,1])[1])){
             return(names(count)[i])        
        } 

    }

    return(names(which(count == max(count)))[1])
}

rT1 = terms(as.formula(paste("~", blk.str1, sep = "")), keep.order = TRUE)  #random terms phase 1
fT1 = terms(as.formula(paste("~", trt.str1, sep = "")), keep.order = TRUE)  #fixed terms  phase 1
blkTerm1 = attr(rT1, "term.labels")


ani.limit = (n/sapply(phase2Design.df, nlevels))/nrow(phase1Design.df)
trt.limit = (n/sapply(phase2Design.df, nlevels))/nlevels(phase1Design.df[, attr(fT1, "term.labels")])

ani.char = as.character(phase1Design.df[, 1])
count = rep(tRep, nrow(phase1Design.df))
names(count) = ani.char

design.df = cbind(phase2Design.df, phase1Design.df)

design.df[, colnames(phase1Design.df)] = "-1"

for (i in 1:nrow(design.df)) {

    design.df[i, colnames(phase1Design.df)[1]] = 
        fillIn(design.df[i, 1:ncol(phase2Design.df)], phase1Design.df, design.df, count, ani.limit,
        trt.limit)

    design.df[i, colnames(phase1Design.df)[2]] = 
        as.character(phase1Design.df[, 2][which(phase1Design.df[, 1] == design.df[i, colnames(phase1Design.df)[1]])])
        
    count[design.df[i, colnames(phase1Design.df)[1]]] = count[design.df[i, colnames(phase1Design.df)[1]]] - 1


    if (any(count == 0))
        count = count[-which(count == 0)]
}


################################################################################



Z1.mat = makeBlkDesMat(design.df, rev(blkTerm1))[[blkTerm1]]




trtTerm1 = attr(fT1, "term.labels")



# Parameter's of treatment
trt.rep = nrow(design.df)/6

X.trt = makeBlkDesMat(design.df, rev(trtTerm1))[[trtTerm1]]


new.Z1.mat = Z1.mat[optThreeStage(init = 1:nrow(design.df), iter = 1000, obj.fun = obj.fun.CRD, 
    test.obj.fun = test.obj.fun.CRD, swap.stage1.new = swap.stage1.new,
    swap.stage2.new = swap.stage2.new, swap.stage3.new = swap.stage3.new, 
    Z1.mat = Z1.mat, X.trt = X.trt, nPlot = 4, Z1.rep = 2, trt.rep = trt.rep,
    blk.proj = blk.proj, nTrt = 6, C.cage = NA, cage.Rep = NA, C.ani = NA, 
    ani.Rep = NA, newResDF = NA, resDF = FALSE, upperValue = 1), ]


design.df$Ani = as.factor(apply(new.Z1.mat, 1, function(x) LETTERS[which(x==1)]))
design.df$Trt = phase1Design.df[match( design.df$Ani, phase1Design.df[, 1]),2]


summaryAovTwoPhase(design.df, blk.str2 = blk.str2, blk.str1 = blk.str1, trt.str = "Tag + Trt")

