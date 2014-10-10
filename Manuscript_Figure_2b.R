## Script creating the figure 2 - second visualization type - for the IMMPUTE paper
#
# Copyright 08/11/2014 Antoine Lizee @ UCSF.
# antoine.lizee@ucsf.edu

options(java.parameters="-Xmx4g")
require(XLConnect)

source("Manuscript_palettes.R")

# Prepare the data --------------------------------------------------------

data <- NULL

institutes <- c("UWASH", "FHCRC", "CCHMC", "MCRI")
methods <- c("HIBAG", "MAGPrediction", "e-HLA", "HLA*IMP:02")
loci <- c("A", "B", "C", "DRB1")

for (method_ind in 1:length(methods)) {
  for (locus_ind in 1:length(loci)) {
    df_i <- readWorksheetFromFile("/media/FD/Dropbox/IMMPUTE/Manuscript/Figure 2.xlsx", sheet=method_ind, 
                                  startCol=locus_ind*4-3, endCol=locus_ind*4-1, startRow=3,
                                  header=F, rownames=F)
    data <- rbind(data, data.frame(df_i, Locus=loci[locus_ind], Method=methods[method_ind]))
  }
}

ll <- levels(data$Method)
data$Method <- factor(data$Method, levels = ll[order(as.character(ll))])

colnames(data)[1:3] <- c("Threshold", "CR", "Perf")
data <- data[data$CR!=0,]

giveRowAndCols <- function(df) {
  df$row <- (!df$Locus %in% c("A", "C")) + 1
  df$col <- (!df$Locus %in% c("A", "B")) + 1
  return(df)
}

data <- giveRowAndCols(data)

locus_data <- data.frame(Locus=loci)
locus_data <- giveRowAndCols(locus_data)
locus_data$y_text <- ifelse(locus_data$row-1, 98, 99.5)
locus_data$x_line_text <- ifelse(locus_data$row-1, 85, 65)
locus_data$y_line_text <- ifelse(locus_data$row-1, 91.8, 90.5)
require(plyr)
perf0.5_data <- ddply(data, ~ Locus + Method, function(df) df[which.min(abs(df$Threshold-0.5)),c("CR", "Perf")])
perf0.5_data <- giveRowAndCols(perf0.5_data)


# Changes the data structure ----------------------------------------------

datab <- data
datab$Recall <- datab$CR * datab$Perf
datab <- datab[datab$Locus == "DRB1",]

# Plot performance recall -------------------------------------------------

require(ggplot2)
require(grid)
g1b <- ggplot(datab) + theme_bw(base_size=11) +
#   geom_text(data=locus_data, aes(label="90% Performance", y=y_line_text, x=x_line_text), size=3.5, color="grey60") +
  geom_line(aes(x=Recall, y=Perf*100, color=Method), size=0.7) +
#   facet_grid(row~col, scales="free_y", ) +
#   geom_text(data=locus_data, aes(label=Locus, y=y_text), x=75, size=4.5 ) +
#   geom_point(data=perf0.5_data, aes(x=Recall, y=Perf*100, fill=Method), color="black", size=2.5, shape=23) +
  labs(y="Performance (%)", x="Global Recall (%)") +
#   coord_cartesian(xlim=c(103,47)) +
  scale_x_reverse() +
  palette_perso +
  palette_perso_fill +
  geom_hline(y=90, linetype=2, size=0.8, color="grey85") +
#   theme(strip.background = element_blank(), strip.text = element_blank(),
#         panel.grid.major = element_line(colour = "grey85", size = 0.2),
#         axis.title.x = element_text(vjust = 0.5),
#         axis.title.y = element_text(vjust = 0.6, angle=90),
#         panel.border = element_rect(fill=NA, colour="grey40", size = 0.2),
#         #         panel.background = element_rect(fill = "grey98", color = NA),
#         panel.margin = unit(0.2, "lines"), # 
#         #         plot.title = element_blank(), # USELESS if no title specified
#         legend.margin = grid::unit(-2.0, "lines"), 
#         legend.text = element_text(size = rel(0.7)),
#         legend.key = element_rect(colour = NA),
#         #         legend.key.size = unit(0.2, "lines"),
#         legend.key.height = unit(0.4, "lines"),
#         legend.key.width = unit(0.3, "lines"),
#         legend.title = element_text(size = rel(0.8), face = "bold", hjust = 0, vjust = -0.2),
#         legend.title.align = 0.3,
#         #         legend.background = element_rect(colour = "black"), # DEBUG ONLY
#         plot.margin = unit(c(-0.1, 0.4, 0.11, 0.08), "lines"), ## INTERACTS WITH LEGEND MARGIN on the right, for the second element
#         axis.ticks.margin = unit(c(0.02), "lines"),
#         axis.ticks = element_line(size = 0.1),
#         axis.ticks.length = unit(0.15, "lines"))

png("/media/FD/Dropbox/IMMPUTE/Manuscript/Figure_Perf_Recall.png", w=12, h=8, units="cm", res=300)
print(g1b)
dev.off()
