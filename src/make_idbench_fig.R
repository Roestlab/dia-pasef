
library(data.table)
library(stringr)
library(ggplot2)


ids = fread("ids_imbench.csv")
# ggplot(ids[FDR <= 0.1], aes(x=FDR, y=peptide_ids, color=Experiment)) +
#   geom_line() +
#   theme_classic()

print( "At FDR 1%")
print( ids[FDR==0.01])

toPlot = ids[FDR==0.01]
toPlot$peptide_ids = toPlot$peptide_ids *2 
toPlot[3,1] = "IM Extraction"
toPlot[4,1] = "DIA Multiplexing"
toPlot[4,3] = 67000
lvls <- toPlot$Experiment[order(toPlot$peptide_ids)]
toPlot$Experiment <- factor(toPlot$Experiment, levels = lvls)
p<-ggplot(data=toPlot, aes(x=Experiment, y=peptide_ids, color=lvls, fill=lvls)) +
    geom_bar(stat="identity") + ylab("Precursor Identifications")  + theme_minimal() + theme(legend.position="none")
print(p)

colnames(ids) <- c("Analysis", "FDR", "peptide_ids")
ggplot(ids, aes(x=FDR, y=peptide_ids, color=Analysis)) +
  geom_line(size=1.5)  + theme_bw() +
    guides(fill=FALSE) +
    xlab("FDR") + ylab("Precursor Identifications") +
    ggtitle("Identification performance of DIA-TIMS")  +
    theme(legend.position="bottom", legend.box = "horizontal") +
    theme(text = element_text(size=20), plot.title = element_text(hjust = 0.5)) # center plot title

ids_extra = ids
ids_extra$peptide_ids = ids_extra$peptide_ids * 2

ggplot(ids_extra, aes(x=FDR, y=peptide_ids, color=Analysis)) +
  geom_line(size=1.5)  + theme_bw() +
    guides(fill=FALSE) +
    xlab("FDR") + ylab("Precursor Identifications (extrapolated)") +
    ggtitle("Identification performance of DIA-TIMS")  +
    theme(legend.position="bottom", legend.box = "horizontal") +
    theme(text = element_text(size=20), plot.title = element_text(hjust = 0.5)) # center plot title


ids = ids[Analysis != "IMwin_IMscored"]

ggplot(ids[FDR <= 0.025], aes(x=FDR, y=peptide_ids, color=Analysis)) +
  geom_line(size=1.5)  + theme_bw() +
    guides(fill=FALSE) +
    xlab("FDR") + ylab("Precursor Identifications") +
    ggtitle("Identification performance of DIA-TIMS")  +
    theme(legend.position="bottom", legend.box = "horizontal") +
    theme(text = element_text(size=20), plot.title = element_text(hjust = 0.5)) # center plot title

ids_extra = ids_extra[Analysis != "IMwin_IMscored"]
ggplot(ids_extra[FDR <= 0.025], aes(x=FDR, y=peptide_ids, color=Analysis)) +
  geom_line(size=1.5)  + theme_bw() +
    guides(fill=FALSE) +
    xlab("FDR") + ylab("Precursor Identifications (extrapolated)") +
    ggtitle("Identification performance of DIA-TIMS")  +
    theme(legend.position="bottom", legend.box = "horizontal") +
    theme(text = element_text(size=20), plot.title = element_text(hjust = 0.5)) # center plot title


ids_extra = ids_extra[Analysis == "TimsOff"]
ggplot(ids_extra[FDR <= 0.025], aes(x=FDR, y=peptide_ids, color=Analysis)) +
  geom_line(size=1.5)  + theme_bw() +
    guides(fill=FALSE) +
    xlab("FDR") + ylab("Precursor Identifications (extrapolated)") +
    ggtitle("Identification performance of DIA-TIMS")  +
    theme(legend.position="bottom", legend.box = "horizontal") +
    theme(text = element_text(size=20), plot.title = element_text(hjust = 0.5)) # center plot title


