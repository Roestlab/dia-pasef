###### Finding IM coverage of the old library and then the new library ###################
library(stringr)
library(ggplot2)
library(data.table)


#Load data
file <- list.files('.', pattern="scored.tsv", full.names=T)
lib_file <- list.files('.', pattern = "hela(.*).tsv", full.names=T)

lib <- fread(lib_file)
lib <- lib[Decoy == 0]
lib$ModifiedPeptideSequence <- gsub("^\\.","",lib$ModifiedPeptideSequence)
stat <- lapply(file, fread)

stat_rep <- lapply(file, function(x){
  name <- strsplit(basename(x),"_")[[1]][16]
})

names(stat) <- stat_rep

#subset data
peptide_stat <- lapply(1:3, function(x){
  rep <- subset(stat[[x]], select=c("transition_group_id","run_id" , "RT","id", "Sequence", "FullPeptideName" ,
                            "Charge", "m/z","Intensity", "decoy", "ProteinName","q_value","pep","peak_group_rank" ))
  rep <- rep[Intensity != 0]
  rep <- rep[decoy == 0]
  rep <- rep[peak_group_rank == 1]
  rep <- rep[q_value <= 0.01]
  rep[, rep := names(stat)[x]]
  rep
})

fullstat <- do.call('rbind', peptide_stat)

# create windows
win_start <- seq(400,1200-12.5, 12.5)
win_end <- seq(400+12.5, 1200, 12.5)

win_mz <- data.frame(start=win_start,
                     end = win_end,
                     num_win = seq(1, 64))




for(i in 1:length(fullstat$`m/z`)){
  fullstat$win_mz[i] <- which( win_mz$start < fullstat$`m/z`[i] & win_mz$end > fullstat$`m/z`[i])
} # takes too long, change it!!


fullstat$win_mz <- lapply(fullstat$`m/z`, function(x){
  s <- which( win_mz$start < x & win_mz$end > x)
})


statplot <- ggplot(fullstat)
statplot + geom_point(aes(x=`m/z`, y=RT))

fullstat$lib_im <- lapply(1:nrow(fullstat), function(a){
  im <- unique(lib[which(lib$ModifiedPeptideSequence == fullstat[a]$FullPeptideName & 
                           round(lib$PrecursorMz,3) == round(fullstat[a]$`m/z`,3) &lib$PrecursorCharge == fullstat[a]$Charge)]$PrecursorIonMobility)
  im
})

fullstat$lib_im <- numeric(length(fullstat$transition_group_id))

for (i in 1:nrow(fullstat)) {
  fullstat$lib_im[i] <- unique(lib[which(lib$ModifiedPeptideSequence == fullstat[i]$FullPeptideName & 
                                        round(lib$PrecursorMz,3) == round(fullstat[i]$`m/z`,3))]$PrecursorIonMobility)
    
}

#define functions to plot distributions of IM values

imcutoffs <- function(window_number, data){
  # Returns the im window cutoffs at 95 quantile for each mz window
  df <- data[win_mz == window_number]
  x <- quantile(df$lib_im, probs = c(.025,.975))
  return(x)
}

window_cutoffs <- NULL
lib.short <- lib[!duplicated(lib$TransitionGroupId)]
for(i in win_mz$num_win){
  window_cutoffs <- rbind(window_cutoffs, quantile(fullstat[win_mz== i]$lib_im, probs=c(0.025, 0.975)))
}


# Calculate for each range, how many peptides are identified for the runs, how many are identified for lib

coverage_ratio <- data.frame(im_25 = window_cutoffs[,1],
                             im_975= window_cutoffs[,2],
                             experiment = numeric(64),
                             library = numeric(64))



lib.short$win_mz <- lapply(lib.short$PrecursorMz, function(x){
  s <- which( win_mz$start < x & win_mz$end > x)
})
fullstat.short <- fullstat[!duplicated(fullstat$transition_group_id)]

for (k in 1:64) {
  coverage_ratio$experiment[k] <- sum(fullstat.short[win_mz == k]$lib_im > coverage_ratio$im_25[k] & fullstat.short[win_mz == k]$lib_im < coverage_ratio$im_975[k])
  coverage_ratio$library[k] <- sum(lib.short[win_mz == k]$PrecursorIonMobility > coverage_ratio$im_25[k] & lib.short[win_mz == k]$PrecursorIonMobility < coverage_ratio$im_975[k])
}

coverage_ratio$window <- win_mz$num_win

coverage_ratio$percent <- coverage_ratio$experiment/coverage_ratio$library

ggplot(coverage_ratio, aes(x = window, y = percent)) +
  geom_bar(stat='identity', fill='cornflowerblue', alpha = 0.8)


library(tidyr)
df <- data.frame(experiment = coverage_ratio$experiment,
                 library = coverage_ratio$library,
                 window = coverage_ratio$window)
rownames(df) <- 1:64

df.long <- reshape2::melt(df, id.vars="window")

ggplot(df.long, aes(x = window, y=value, fill=variable)) +
  geom_bar(stat='identity', position='dodge' ,alpha=0.8)


#coverage ratios:

sum(coverage_ratio$experiment)/sum(coverage_ratio$library)

# extrapolated assumption of optimal id from coverage ratio:
sum(coverage_ratio$experiment)/sum(coverage_ratio$library) * 91461



