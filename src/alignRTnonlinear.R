library(fANCOVA)
library(stringr)
library(ggplot2)
library(data.table)
## Script to align the retention times of a (fractionated) library with loess

args = commandArgs(trailingOnly=TRUE)

if (length(args)<4) {
  stop("At least four arguments must be supplied: evidence.txt CiRT.txt pdfout.pdf output_name ", call.=FALSE)
} else{
  evidence <- args[1]
  irt <- args[2]
  pdfout <- args[3]
  outfile <- gsub("*\\.tsv|*\\.txt|$", ".tsv", args[4])
}
### Function definitions

alignFractions <- function(lib, pirt,optimizeSpan=c("once", "always"), span=NULL, plot=T, pdfout="alignment.pdf",...){
  optimizeSpan <- match.arg(optimizeSpan)
  fractions <- unique(lib$`Raw file`)
  libAnn <- data.table(NULL)

  if(plot){
    pdf(pdfout)
  }
  for(fraction in fractions){
    fr <- alignFraction(lib[`Raw file` == fraction], pirt, plot=plot, ...)

    if(nrow(libAnn) > 0){
      libAnn <- rbind(libAnn, fr$lib)
    }else{
      libAnn <- fr$lib
    }
    if(optimizeSpan == "once"){
      span <- fr$param$span
    }
  }
  if(plot){
    dev.off()
  }
  return(libAnn)
}

alignFraction <- function(lib_fr, pirt, holdout_size=0.1,
                          outlier_tol=7, plot=T, span=NULL,
                          optimization_criterion=c("gcv", "aicc"),
                          PDF=NULL){
  ocrit <- match.arg(optimization_criterion)
  ## Get evidence into shape
  lib_fr$`Modified sequence`<-str_replace_all(lib_fr$`Modified sequence`,"_","")
  lib_fr$`Modified sequence`<-str_replace_all(lib_fr$`Modified sequence`,"(ox)","Oxidation")
  lib_fr$`Modified sequence`<-str_replace_all(lib_fr$`Modified sequence`,"(ph)","Phospho")
  lib_fr$`Modified sequence`<-str_replace_all(lib_fr$`Modified sequence`,"C","C\\(Carbamidomethyl\\)")
  lib_fr$`Modified sequence`<-str_replace_all(lib_fr$`Modified sequence`,"(ac)","Acetylation")
  lib_frBest <- lib_fr[lib_fr[,.I[which.min(PEP)], by=.(`Modified sequence`, Charge)]$V1]
  ## Get irts into shape
  pirt <- copy(pirt)
  setnames(pirt, "ModifiedPeptideSequence", "sequence")
  pirt$sequence <- str_replace_all(pirt$sequence, "\\(UniMod:4\\)", "\\(Carbamidomethyl\\)")
  pirt$sequence <- str_replace_all(pirt$sequence, "\\(UniMod:1\\)", "Acetylation")
  pirt$sequence <- str_replace_all(pirt$sequence, "\\(UniMod:35\\)", "Oxidation")
  setnames(pirt, "sequence", "Modified sequence")
  setnames(pirt, "PrecursorCharge", "Charge")
  pirt <- unique(pirt[, .(`Modified sequence`, Charge, NormalizedRetentionTime)])

  alignTable <- merge(pirt, lib_frBest, by = c("Modified sequence", "Charge"))
  alignTable <- alignTable[!(is.na(`Retention time`) | is.na(NormalizedRetentionTime))]
  alignTable <- alignTable[order(`Retention time`)]
  reg <- lm(NormalizedRetentionTime ~ `Retention time` , data = alignTable)
  outliers <- which(abs(resid(reg)) >= outlier_tol)
  alignTable$outlier <- FALSE
  alignTable[outliers, outlier := TRUE]
  message(paste("Removed", length(outliers), "Outliers by linear regression."))
  holdoutRows <- sample(1:nrow(alignTable[-nrow(alignTable)][-1][!(outlier)]),
                        floor(nrow(alignTable[!(outlier)])*holdout_size), replace = F)
  alignTable[, holdout := FALSE]
  alignTable[holdoutRows, holdout := TRUE]
  alignTableHoldout <- alignTable[holdoutRows]
  alignTableTrain <- alignTable[-holdoutRows][!(outlier)]
  ## Regress linear and fit loess on training set
  reg <- lm(NormalizedRetentionTime ~ `Retention time` , data = alignTableTrain)
  ## surface direct lets you do extrapolation in cases where the extreme values are outliers
  ## The plot should always be inspected for weird behaviour at the edges
  loessMod <- loess.as(alignTableTrain$`Retention time`,
                       alignTableTrain$NormalizedRetentionTime, user.span = span,
                       criterion = ocrit, plot=F,
                       control=loess.control(surface = "direct"))
  message(paste("Span value: ", loessMod$pars$span))
  ## Evaluate fit on holdout set
  loessR <- try(alignTableHoldout$NormalizedRetentionTime - predict(loessMod, alignTableHoldout$`Retention time`), silent=T)
  loessR <- sum(loessR^2)
  lmR <- try(alignTableHoldout$NormalizedRetentionTime -
             (reg$coefficients[1] + reg$coefficients[2] * alignTableHoldout$`Retention time`), silent=T)
  lmR <- sum(lmR^2)
  meanR <- alignTableHoldout$NormalizedRetentionTime - mean(alignTableHoldout$NormalizedRetentionTime)
  meanR <- sum(meanR^2)
  loessRsq <- 1- loessR/meanR
  lmRsq <- 1- lmR/meanR
  message(paste0("Linear fit Rsq: ", lmRsq, ", Squared sum of residuals: ", lmR))
  message(paste0("Loess fit Rsq: ", loessRsq, ", Squared sum of residuals: ", loessR))

  normRT <- predict(loessMod, lib_fr$`Retention time`)
  ## prediction_errors <- is.na(normRT)
  ## if(any(prediction_errors)){
  ##   normRT[prediction_errors] <- reg$coefficients[1] +
  ##     reg$coefficients[2] * alignTable$`Retention time.y`[prediction_errors]
  ## }
  file <- unique(lib_fr$`Raw file`)
  p <- ggplot(alignTable, aes(x=`Retention time`, y=NormalizedRetentionTime)) +
    geom_point(aes(color=outlier), size=0.8) +
    scale_colour_manual(values = c("black", "red")) +
    geom_abline(intercept = reg$coefficients[1], slope = reg$coefficients[2]) +
    geom_line(data = lib_fr, aes(x=`Retention time`, y=normRT), color="blue", size=1) +
    ggtitle(paste0("File: ", file, "\n",
                   "Linear fit Rsq: ", lmRsq, ", Squared sum of residuals: ", lmR, "\n",
                   "Loess fit Rsq: ", loessRsq, ", Squared sum of residuals: ", loessR))

  residTable <- data.table(RT=alignTableTrain$`Retention time`,
                           linearAlignment=resid(reg),
                           loessAlignment=resid(loessMod))
  p2 <- ggplot(residTable, aes(x=RT)) +
    geom_point(aes(y=linearAlignment), color="red", size=0.8) +
    geom_point(aes(y=loessAlignment), color="black", size=0.8) +
    geom_hline(yintercept=0) +
    ylab("Fit residuals (s)")
  if(plot){
    if(!is.null(PDF)){
      pdf(PDF)
    }
    plot(p)
    plot(p2)
    if(!is.null(PDF)){
      dev.off()
    }
  }
  lib_fr$NormalizedRetentionTime <- normRT
  return(list(pars=loessMod$pars, lib=lib_fr))
}


## Main


fracLib <- fread(evidence)
cirts <- fread(irt)
# sort the data by PEP
fracLib <- fracLib[fracLib[,.I[which.min(PEP)], by= c("Charge", "Modified sequence", "Raw file")]$V1]
# Record the protein names and the id numbers
ProteinNames <- subset(fracLib, select=c(`V1`, `Protein names`))
# Clean up the protein names and remove all the \t
ProteinNames$`Protein names` <- gsub("\t", "", ProteinNames$`Protein names`)

# Remove the protein names column to avoid extra tabs
fracLib <- subset(fracLib, select= -`Protein names`)

# Testing run with one sample run
files <- unique(fracLib$`Raw file`)
lib <- fracLib[`Raw file` == files[10]]
fr <- alignFraction(lib, pirt = cirts, optimization_criterion = "gcv", PDF=NULL)

# Alignment for the entire file and all runs
libAnn <- alignFractions(fracLib, pirt = cirts, optimizeSpan = "always", pdfout = pdfout)

# Merge the protein names back with the aligned library by their id number
libAnn <- merge(libAnn, ProteinNames, by="V1")

write.table(libAnn, outfile, row.names = F, quote = T, sep = "\t")
