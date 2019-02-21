
library(data.table)
library(stringr)
library(ggplot2)


ids = fread("ids_imbench.csv")
# ggplot(ids[FDR <= 0.1], aes(x=FDR, y=peptide_ids, color=Experiment)) +
#   geom_line() +
#   theme_classic()

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



