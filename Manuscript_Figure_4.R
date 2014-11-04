## Script creating the figure 2 - second visualization type - for the IMMPUTE paper
#
# Copyright 08/11/2014 Antoine Lizee @ UCSF.
# antoine.lizee@ucsf.edu

options(java.parameters="-Xmx4g")
library(XLConnect)
library(ggplot2)
require(reshape2)
require(plyr)

source("Manuscript_palettes.R")


# Prepare the data --------------------------------------------------------

# Find by hand the correct paths
parentDir <- "~/Dropbox/IMMPUTE/Imputation Scoring"
dir(pattern = "M1AllLoci_MinCut", path = parentDir, full.names = T, recursive = T)
filePaths <- matrix(c("~/Dropbox/IMMPUTE/Imputation Scoring/UWASH/Experiment 2/Score_2F_NoMask/4D_UWASHexpII_PcurveM1AllLoci_MinCut.txt", "UWASH",
                      "~/Dropbox/IMMPUTE/Imputation Scoring/FHCRC/Experiment 2/Score_2F_NoMask/4D_FHCRC_PcurveM1AllLoci_MinCut.txt", "FHCRC",
                      "~/Dropbox/IMMPUTE/Imputation Scoring/CCHMC/SVM AI/Experiment 2/Score_2F_NoMask/4D_CCHMC_SVMAI_PcurveM1AllLoci_MinCut.txt", "CCHMC",
                      "~/Dropbox/IMMPUTE/Imputation Scoring/MCRI/Experiment 2/Score_2F_NoMask/4D_MCRIexp2_PcurveM1AllLoci_MinCut.txt", "MCRI"),
                    ncol = 2, byrow = T)

## Some verifications about coherence of the data we get
res <- read.delim(filePaths[1,1])
cbind(res$SubCount / res$SubCount[1], res$CallRate) ##OK
cbind(res$CorrSub.Count / res$CorrSub.Count[1], res$CorrSub.TotRatio) ##OK
cbind(res$CorrSub.Count / res$SubCount, res$CorrSub.Ratio) ##OK


##Import all the data

data_subject <- NULL

institutes <- c("UWASH", "FHCRC", "CCHMC", "MCRI")
methods <- c("HIBAG", "MAGPrediction", "e-HLA", "HLA*IMP:02")

for (method_ind in 1:length(methods)) {
  res_i <- read.delim(filePaths[method_ind,1]) 
  df_i <- data.frame( Threshold = res_i$Cut,
                      TotAccuracy = res_i$CorrSub.Count / res_i$SubCount[1],
                      SubsetAccuracy = res_i$CorrSub.Count / res_i$SubCount,
                      CR = res_i$SubCount / res_i$SubCount[1] * 100)
  df_i$MaxTot <- pmin(df_i$TotAccuracy[1], df_i$CR/100)
  data_subject <- rbind(data_subject, data.frame(df_i, Method=methods[method_ind]))
}

##Reorder the methods
ll <- levels(data_subject$Method)
data_subject$Method <- factor(data_subject$Method, levels = ll[order(as.character(ll))])

require(plyr)
perf0.5_data <- ddply(data_subject, ~ Method, function(df) df[which.min(abs(df$Threshold-0.5)),c("CR", "SubsetAccuracy", "TotAccuracy")])

giveRowAndCols <- function(df) {
  df$row <- (df$Method %in% c("e-HLA", "HLA*IMP:02")) + 1
  df$col <- (df$Method %in% c("e-HLA", "HIBAG")) + 1
  return(df)
}

data_subject <- giveRowAndCols(data_subject)
perf0.5_data <- giveRowAndCols(perf0.5_data)

# Plots -------------------------------------------------

require(ggplot2)
require(grid)

g1 <- ggplot(data_subject) + theme_bw(base_size=16) +
  #   geom_text(data=locus_data, aes(label="90% Performance", y=y_line_text, x=x_line_text), size=3.5, color="grey60") +
  #   geom_hline(y=90, linetype=2, size=0.8, color="grey85") +
  geom_line(aes(x=CR, y=SubsetAccuracy*100, color=Method), size=0.6) +
  #   geom_text(data=locus_data, aes(label=Locus, y=y_text), x=75, size=4.5 ) +
  geom_point(data=perf0.5_data, aes(x=CR, y=SubsetAccuracy*100, fill=Method), color="black", size=2.5, shape=23) +
  labs(y="Performance (%)", x="Call Rate(%)") +
  #   coord_cartesian(xlim=c(47,103)) +
  scale_x_reverse() +
  palette_perso +
  palette_perso_fill 
#   theme(strip.background = element_blank(), strip.text = element_blank(),
#         panel.grid.major = element_line(colour = "grey85", size = 0.2),
#         axis.title.x = element_text(vjust = 0),
#         axis.title.y = element_text(vjust = 0.4, angle=90),
#         panel.border = element_rect(fill=NA, colour="grey40", size=0.5),
#         #         panel.background = element_rect(fill = "grey98", color = NA),
#         #         panel.margin = unit(0.25, "lines"), # DEFAULT is good
#         #         plot.title = element_blank(), # USELESS if no title specified
#         legend.margin = unit(-1.5, "lines"), 
#         legend.text = element_text(size = rel(0.9)),
#         legend.key = element_rect(colour = NA),
#         legend.key.height = unit(1.0, "lines"),
#         legend.key.width = unit(1.0, "lines"),
#         legend.title = element_blank(), #element_text(size = rel(0.8), face = "bold", hjust = 0),
#         legend.title.align = 0.2,
#         #         legend.background = element_rect(colour = "black"), # DEBUG ONLY
#         plot.margin = unit(c(0.1, 0.5, 0.6, 0.5), "lines")) ## INTERACTS WITH LEGEND MARGIN on the right, for the second element

g2 <- ggplot(data_subject) + theme_bw(base_size=16) +
  #   geom_text(data=locus_data, aes(label="90% Performance", y=y_line_text, x=x_line_text), size=3.5, color="grey60") +
  #   geom_hline(y=90, linetype=2, size=0.8, color="grey85") +
  geom_line(aes(x=CR, y=TotAccuracy *100, color=Method), size=0.6) +
  #   geom_text(data=locus_data, aes(label=Locus, y=y_text), x=75, size=4.5 ) +
  geom_point(data=perf0.5_data, aes(x=CR, y=TotAccuracy*100, fill=Method), color="black", size=2.5, shape=23) +
  labs(y="Total Accuracy (% of the whole dataset)", x="Call Rate(%)") +
  #   coord_cartesian(xlim=c(47,103)) +
  scale_x_reverse() +
  palette_perso +
  palette_perso_fill 

g3 <- ggplot(data_subject) + theme_bw(base_size=16) +
  #   geom_text(data=locus_data, aes(label="90% Performance", y=y_line_text, x=x_line_text), size=3.5, color="grey60") +
  #   geom_hline(y=90, linetype=2, size=0.8, color="grey85") +
  geom_line(aes(x=CR, y=TotAccuracy *100, color=Method), size=0.6) +
  geom_line(aes(x=CR, y=SubsetAccuracy*100, color=Method), linetype = "31", size=0.7) +
  #   geom_text(data=locus_data, aes(label=Locus, y=y_text), x=75, size=4.5 ) +
  geom_point(data=perf0.5_data, aes(x=CR, y=SubsetAccuracy*100, fill=Method), color="black", size=2.5, shape=23) +
  labs(y="Accuracies (%)", x="Call Rate(%)") +
  #   coord_cartesian(xlim=c(47,103)) +
  scale_x_reverse() +
  palette_perso +
  palette_perso_fill 

g4 <- ggplot(data_subject) + theme_bw(base_size=16) +
  #   geom_text(data=locus_data, aes(label="90% Performance", y=y_line_text, x=x_line_text), size=3.5, color="grey60") +
  #   geom_hline(y=90, linetype=2, size=0.8, color="grey85") +
  geom_line(aes(x=TotAccuracy *100, y=SubsetAccuracy *100, color=Method), size=0.6) +
  #   geom_text(data=locus_data, aes(label=Locus, y=y_text), x=75, size=4.5 ) +
  geom_point(data=perf0.5_data, aes(x=TotAccuracy *100, y=SubsetAccuracy *100, fill=Method), color="black", size=2.5, shape=23) +
  labs(y="Performance (or Subset Accuracy) (%)", x="Total accuracy (%)") +
  #   coord_cartesian(xlim=c(47,103)) +
  scale_x_reverse() +
  palette_perso +
  palette_perso_fill 



# Export (mainly from figure 1) -------------------------------------------


theme_perso <-   theme(
  axis.title = element_text(size = rel(0.9)),
  axis.title.x = element_text(vjust = 0),
  axis.title.y = element_text(vjust = 0.4, angle=90)
)

theme_perso2 <- theme (
  panel.grid.major = element_line(colour = "grey85", size = 0.2),
  axis.title = element_text(size = rel(0.9)),
  axis.title.x = element_text(vjust = 0.5),
  axis.title.y = element_text(vjust = 0.6, angle=90),
  panel.border = element_rect(fill=NA, colour="grey40", size = 0.2),
  #         panel.background = element_rect(fill = "grey98", color = NA),
  panel.margin = unit(0.2, "lines"), # 
  #         plot.title = element_blank(), # USELESS if no title specified
  legend.margin = grid::unit(-0.9, "lines"), 
  legend.text = element_text(size = rel(0.8)),
  legend.key = element_rect(colour = NA),
  legend.key.size = unit(0.5, "lines"),
  #   legend.key.height = unit(0.4, "lines"),
  #   legend.key.width = unit(0.4, "lines"),
  legend.title = element_text(size = rel(0.85), face = "bold", hjust = 0, vjust = -0.2),
  legend.title.align = 0,
  #           legend.background = element_rect(colour = "black"), # DEBUG ONLY
  plot.margin = unit(c(0.4, 0., 0.11, 0.08), "lines"), ## INTERACTS WITH LEGEND MARGIN on the right, for the second element
  axis.ticks.margin = unit(c(0.05), "lines"),
  axis.ticks = element_line(size = 0.1),
  axis.ticks.length = unit(0.25, "lines"))

theme_perso <- theme_bw(base_size=7) %+replace% theme_perso2 #AL OVERWRITE for final output

exportFigure1 <- function (name, ggp) {
  #   png(paste0("/media/FD/Dropbox/IMMPUTE/Manuscript/", name,".png"), w=620, h=500, res=100) # Original
  png(paste0("/media/FD/Dropbox/IMMPUTE/Manuscript/", name,".png"),  w=8.9, h=7.2, units="cm", res=300)
  print(ggp + theme_perso)
  dev.off()
  tiff(paste0("/media/FD/Dropbox/IMMPUTE/Manuscript/", name,".tiff"),  w=8.9, h=7.2, units="cm", res=300)
  print(ggp + theme_perso)
  dev.off()
}

exportFigure1(name = "Figure_Subject_1_perf_CR", g1)
exportFigure1(name = "Figure_Subject_2_totPerf_CR", g2)
exportFigure1(name = "Figure_Subject_3_allPerfs_CR", g3)
exportFigure1(name = "Figure_Subject_4_perf_totPerf", g4)


# Confidence metric assessment ------------------------------------------------------------

areas <- ddply(data_subject, ~Method, function(dfi) {
  totArea <- -sum(dfi$MaxTot[-1] * diff(dfi$CR))
  subArea <- -sum(dfi$TotAccuracy[-1] * diff(dfi$CR))
  data.frame(Tot = totArea,
             Sub = subArea,
             Diff = totArea - subArea,
             Ratio = (totArea - subArea)/totArea)
})

areas$x <- c(30, 45, 25, 35)
areas$y <- c(15, 25, 13, 17)

g5 <- ggplot(data_subject) + theme_bw(base_size=16) +
  geom_line(aes(x=CR, y=MaxTot *100), color="black", size=0.6) +
  geom_line(aes(x=CR, y=TotAccuracy *100, color=Method), size=0.6) +
  geom_point(data=perf0.5_data, aes(x=CR, y=TotAccuracy*100, fill=Method), color="black", size=2.5, shape=23) +
  geom_text(data = areas, aes(label = round(Ratio,3), x = x, y=y), size = rel(3), color = "grey30")+
  labs(y="Total Accuracy (% of the whole dataset)", x="Call Rate(%)") +
  scale_x_reverse() +
  palette_perso +
  palette_perso_fill +
  #   facet_grid(row ~ col)
  facet_grid(~Method)

## Custom export inspired from above

theme_perso_below <- theme_perso + theme(
                           panel.grid.major = element_line(colour = "grey85", size = 0.2),
                           axis.title.x = element_text(vjust = 0.5),
                           axis.title.y = element_text(vjust = 0.6, angle=90),
                           panel.border = element_rect(fill=NA, colour="grey40", size = 0.2),
                           #         panel.background = element_rect(fill = "grey98", color = NA),
                           panel.margin = unit(0.2, "lines"), # 
                           #         plot.title = element_blank(), # USELESS if no title specified
                           legend.margin = grid::unit(-1.2, "lines"), 
                           legend.text = element_text(size = rel(0.9)),
                           legend.key = element_rect(colour = NA),
                           #         legend.key.size = unit(0.2, "lines"),
                           legend.key.height = unit(0.4, "lines"),
                           legend.key.width = unit(0.5, "lines"),
                           legend.title = element_text(size = rel(1.0), face = "bold", hjust = 0, vjust = -0.2),
                           legend.title.align = 0,
                           legend.position = "bottom",
                           #         legend.background = element_rect(colour = "black"), # DEBUG ONLY
                           plot.margin = unit(c(0.35, 0.25, 0.25, 0.08), "lines"), ## INTERACTS WITH LEGEND MARGIN on the right, for the second element
                           axis.ticks.margin = unit(c(0.02), "lines"),
                           axis.ticks = element_line(size = 0.1),
                           axis.ticks.length = unit(0.15, "lines"))


name = "Figure_Subject_5_confidenceComparison"
ggp = g5
png(paste0("/media/FD/Dropbox/IMMPUTE/Manuscript/", name,".png"),  w=18.3, h=6.5, units="cm", res=300)
print(ggp + theme_perso_below)
dev.off()
