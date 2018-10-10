##### Compare Library RT #####

## Purpose: Compare and find optimal parameters for OpenSwath in Loess Library Generation

##### Libraries Required ########
library(data.table)
library(stringr)
library(ggplot2)
setwd("D:/TIMS-DIAPASEF/Library_Improvement/20180815_LibraryImprovement_CompareRT_Span")
###### Paths ####################
#project_root <- "/project/6011811/frankmak/tims/Mar2018_helaBenchmark/"
rt_data_path <- 'RTcompare'
span_data_path <- 'spancompare'

######### RT comparison ##########

rtfiles <- list.files(rt_data_path, pattern="full_stat.csv", full.names=T)

stat <- lapply(rtfiles, fread)

rt_runs <- unlist(lapply(basename(rtfiles), function(x){
  name <- strsplit(basename(x),'_')[[1]][19]
  name
}))

names(stat) <- rt_runs

stat <- do.call("rbind", lapply(names(stat), function(x){
  res <- stat[[x]]
  res$rt <- x
  res
} ))

View(stat)
ggplot(stat) + geom_bar(aes(x=))


ggplot(stat[fdr <= 0.1], aes(x=fdr, y=tp, color=rt)) +
  geom_line() +
  theme_bw() +
  geom_vline(xintercept=0.01, linetype="dotted")

# ZOOMED:	
ggplot(stat[fdr <= 0.1], aes(x=fdr, y=tp, color=rt)) +
  geom_line() + coord_cartesian(xlim = c(0, 0.02), ylim = c(20000, 23500)) + 
  theme_bw() +
  geom_vline(xintercept=0.01, linetype="dotted")


ggplot(stat[fdr <= 0.01]) + geom_bar(aes(x = rt))
ggplot(stat[fdr <= 0.1], aes(x=tp, fill=rt)) +
  geom_density(alpha = 0.3)


rtlin <- list.files(rt_data_path, pattern = "_scored.tsv", full.names = T)
rtstat <- fread(rtlin[1])

#rtlin2 <- fread(paste0(project_root, "data/openswath/Library_improvement/decoy_comparison/witholdCiRT/20180323_AnBr_SA_diaPASEF_200ng_HeLa_Rost_Method_2_b_TIMSon_02_A2_01_2166_LibShuffle_scored.tsv"))

rtstat <- rtstat[peak_group_rank == 1]



rtstat <- rtstat[q_value <= 0.01]

ggplot(rtstat) + geom_density(aes(x=delta_rt), alpha = 0.4)
hist(rtstat$RT)


######## span comparison #########################
spanfiles <- list.files(span_data_path, pattern="full_stat.csv", full.names=T)

spanstat <- lapply(spanfiles, fread)

span_runs <- unlist(lapply(basename(spanfiles), function(x){
  name <- strsplit(basename(x),'_')[[1]][19]
  name
}))

names(spanstat) <- span_runs

spanstat <- do.call("rbind", lapply(names(spanstat), function(x){
  res <- spanstat[[x]]
  res$span <- x
  res
} ))

View(spanstat)

ggplot(spanstat[fdr <= 0.1 ], aes(x=fdr, y=tp, color=span)) +
  geom_line() +
  theme_bw() +
  geom_vline(xintercept=0.01, linetype="dotted")

#ZOOMED:
ggplot(spanstat[fdr <= 0.1 ], aes(x=fdr, y=tp, color=span)) +
  geom_line() + coord_cartesian(xlim = c(0, 0.02), ylim = c(20000, 23500)) +
  theme_bw() +
  geom_vline(xintercept=0.01, linetype="dotted")

ggplot(spanstat[fdr <= 0.01]) + geom_bar(aes(x = span))
ggplot(spanstat[fdr <= 0.1], aes(x=tp, fill=span)) +
  geom_density(alpha = 0.3)

plot(spanstat$fdr, spanstat$tp)



  
