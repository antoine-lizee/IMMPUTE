## Script creating the figure 2 for the IMMPUTE paper
# This script creates the main validation Figure by plotting a 
# threshold-sliding Accuracy vs call Rate curve.
#
# Copyright 08/11/2014 Antoine Lizee @ UCSF.
# antoine.lizee@ucsf.edu

options(java.parameters="-Xmx4g")
library(XLConnect)
library(ggplot2)
library(grid)
library(plyr)

source("Manuscript_palettes.R")
source("Manuscript_Utilities.R")

inputFolder <- "/media/FD/Dropbox/IMMPUTE/Manuscript/Data Tables and Figures/"
outputFolder <- "/media/FD/Dropbox/IMMPUTE/Manuscript/Data Tables and Figures/All Figures/"

# Prepare the data --------------------------------------------------------

data <- NULL

institutes <- c("UWASH", "FHCRC", "CCHMC", "MCRI")
methods <- c("HIBAG", "MAGPrediction", "e-HLA", "HLA*IMP:02")
loci <- c("A", "B", "C", "DRB1")

for (method_ind in 1:length(methods)) {
  for (locus_ind in 1:length(loci)) {
    df_i <- readWorksheetFromFile(paste0(inputFolder, "Figure 2.xlsx"), sheet=method_ind, 
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

perf0.5_data <- ddply(data, ~ Locus + Method, function(df) df[which.min(abs(df$Threshold-0.5)),c("CR", "Perf")])
perf0.5_data <- giveRowAndCols(perf0.5_data)



# Create the charts -------------------------------------------------------


## Print two columns

g1 <- ggplot(data) + theme_bw(base_size=11) +
  #   geom_abline(a=1, b=0, linetype=2, size=1.2, color="grey85") +
  geom_text(data=locus_data, aes(label="90% Accuracy", y=y_line_text, x=x_line_text), size=3.5, color="grey60") +
  geom_line(aes(x=CR, y=Perf*100, color=Method), size=0.7) +
  facet_grid(row~col, scales="free_y", ) +
  geom_text(data=locus_data, aes(label=paste0("HLA-",Locus), y=y_text), x=-75, size=4.5 ) +
  #   geom_point(data=perf0.5_data, aes(x=CR, y=Perf*100), color="black", size=3, shape=23, fill="yellow") +
  geom_point(data=perf0.5_data, aes(x=CR, y=Perf*100, fill=Method), color="black", size=2.5, shape=23) +
  #   scale_shape_manual(values=c(15,19,18,17)) +
  labs(y="Imputation Accuracy (%)", x="Call Rate (%)") +
  #   scale_y_continuous(expand=c(0,0)) +
  coord_cartesian(xlim=c(47,103)) +
  scale_x_reverse() +
  palette_perso +
  #   scale_color_brewer(type="qual", palette=2) +
  palette_perso_fill +
  #   scale_fill_brewer(type="qual", palette=2, guide="none") +
  geom_hline(y=90, linetype=2, size=0.8, color="grey85") +
  theme(strip.background = element_blank(), strip.text = element_blank(),
        panel.grid.major = element_line(colour = "grey85", size = 0.2),
        axis.title.x = element_text(vjust = 0),
        axis.title.y = element_text(vjust = 0.4, angle=90),
        panel.border = element_rect(fill=NA, colour="grey40", size=0.5),
        #         panel.background = element_rect(fill = "grey98", color = NA),
        #         panel.margin = unit(0.25, "lines"), # DEFAULT is good
        #         plot.title = element_blank(), # USELESS if no title specified
        legend.margin = unit(-1.5, "lines"), 
        legend.text = element_text(size = rel(0.9)),
        legend.key = element_rect(colour = NA),
        legend.key.height = unit(1.0, "lines"),
        legend.key.width = unit(1.0, "lines"),
        legend.title = element_blank(), #element_text(size = rel(0.8), face = "bold", hjust = 0),
        legend.title.align = 0.2,
        #         legend.background = element_rect(colour = "black"), # DEBUG ONLY
        plot.margin = unit(c(0.1, 0.5, 0.6, 0.5), "lines")) ## INTERACTS WITH LEGEND MARGIN on the right, for the second element

printGGplot(g1, paste0(outputFolder, "Figure_AllelesAcc_2cols"), w=18.3, h=12, units="cm", res=300)

####
if (oneColumnRightLegend <- FALSE) {
  g1 <- ggplot(data) + theme_bw(base_size=7) +
    #   geom_abline(a=1, b=0, linetype=2, size=1.2, color="grey85") +
    geom_line(aes(x=CR, y=Perf*100, color=Method), size=0.45) +
    facet_grid(row~col, scales="free_y", ) +
    geom_text(data=locus_data, aes(label=paste0("HLA-",Locus), y=y_text), x=-75, size=2.4 ) +
    geom_text(data=locus_data, aes(label="90% Accuracy", y=y_line_text, x=x_line_text), size=1.7, color="grey55") +
    geom_point(data=perf0.5_data, aes(x=CR, y=Perf*100), color="black", size=1, shape=23, fill="yellow") +
    #   scale_shape_manual(values=c(15,19,18,17)) +
    labs(y="Imputation Accuracy (%)", x="Call Rate (%)") +
    #   scale_y_continuous(expand=c(0,0)) +
    coord_cartesian(xlim=c(47,103)) +
    scale_x_reverse() +
    palette_perso +
    #     scale_color_brewer(type="qual", palette=2) +
    geom_hline(y=90, linetype=2, size=0.4, color="grey70") +
    theme(strip.background = element_blank(), strip.text = element_blank(),
          panel.grid.major = element_line(colour = "grey85", size = 0.2),
          axis.title.x = element_text(vjust = 0.5),
          axis.title.y = element_text(vjust = 0.6, angle=90),
          panel.border = element_rect(fill=NA, colour="grey40", size = 0.2),
          #         panel.background = element_rect(fill = "grey98", color = NA),
          panel.margin = unit(0.2, "lines"), # 
          #         plot.title = element_blank(), # USELESS if no title specified
          legend.margin = grid::unit(-2.0, "lines"), 
          legend.text = element_text(size = rel(0.7)),
          legend.key = element_rect(colour = NA),
          #         legend.key.size = unit(0.2, "lines"),
          legend.key.height = unit(0.4, "lines"),
          legend.key.width = unit(0.3, "lines"),
          legend.title = element_text(size = rel(0.8), face = "bold", hjust = 0, vjust = -0.2),
          legend.title.align = 0.3,
          #         legend.background = element_rect(colour = "black"), # DEBUG ONLY
          plot.margin = unit(c(-0.1, 0.4, 0.11, 0.08), "lines"), ## INTERACTS WITH LEGEND MARGIN on the right, for the second element
          axis.ticks.margin = unit(c(0.02), "lines"),
          axis.ticks = element_line(size = 0.1),
          axis.ticks.length = unit(0.15, "lines"))
  
  printGGplot(g1, paste0(outputFolder, "Figure_AllelesAcc_1col"),  w=8.9, h=6, units="cm", res=300)
}
####

g1 <- ggplot(data) + theme_bw(base_size=7) +
  #   geom_abline(a=1, b=0, linetype=2, size=1.2, color="grey85") +
  geom_line(aes(x=CR, y=Perf*100, color=Method), size=0.65) +
  facet_grid(row~col, scales="free_y", ) +
  geom_text(data=locus_data, aes(label=paste0("HLA-",Locus), y=y_text), x=-75, size=2.4 ) +
  geom_text(data=locus_data, aes(label="90% Accuracy", y=y_line_text, x=x_line_text), size=2, color="grey55") +
  #   geom_point(data=perf0.5_data, aes(x=CR, y=Perf*100), color="black", size=1.5, shape=23, fill="yellow", linewidth = 0.1) + ## YELLOW VERSION
  geom_point(data=perf0.5_data, aes(x=CR, y=Perf*100, fill = Method), color="black", size=1.5, shape=23, alpha = 0.7, linewidth = 0.1) + 
  #   scale_shape_manual(values=c(15,19,18,17)) +
  labs(y="Imputation Accuracy (%)", x="Call Rate (%)", colour = "Method", fill = "Method") +
  #   scale_y_continuous(expand=c(0,0)) +
  scale_x_reverse() +
  coord_cartesian(xlim=c(47,103)) +
  palette_perso +
  palette_perso_fill +
  #     scale_color_brewer(type="qual", palette=2) +
  geom_hline(y=90, linetype=2, size=0.4, color="grey70") +
  theme(strip.background = element_blank(), strip.text = element_blank(),
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
        plot.margin = unit(c(-0.1, -0.1, 0.25, 0.08), "lines"), ## INTERACTS WITH LEGEND MARGIN on the right, for the second element
        axis.ticks.margin = unit(c(0.02), "lines"),
        axis.ticks = element_line(size = 0.1),
        axis.ticks.length = unit(0.15, "lines"))

printGGplot(g1, paste0(outputFolder, "Figure_AllelesAcc_1col_bottom"), w=8.9, h=7.8, units="cm", res=300)


# Try different color schemes ---------------------------------------------

if (colorSchemes.b <- FALSE) {
  
  tryColorScheme <- function(colorScaleNum){
    png(paste0(outputFolder,"colorSchemes/Figure_2_", colorScaleNum, ".png"), w=8.00, h=7.20, units="in", res=300)
    require(ggplot2)
    
    g1 <- ggplot(data) + theme_bw(base_size=16) +
      #   geom_abline(a=1, b=0, linetype=2, size=1.2, color="grey85") +
      geom_line(aes(x=CR, y=Perf*100, color=Method), size=1) +
      facet_grid(row~col, scales="free_y", ) +
      geom_text(data=locus_data, aes(label=paste0("HLA-",Locus), y=y_text), x=-75, size=8) +
      geom_text(data=locus_data, aes(label="90% Accuracy", y=y_line_text, x=x_line_text), size=4, color="grey60") +
      geom_point(data=perf0.5_data, aes(x=CR, y=Perf*100), color="black", size=2, shape=23, fill="yellow") +
      #   scale_shape_manual(values=c(15,19,18,17)) +
      labs(y="Imputation Accuracy (%)", x="Call Rate (%)") +
      #   scale_y_continuous(expand=c(0,0)) +
      coord_cartesian(xlim=c(47,103)) +
      scale_x_reverse() +
      scale_color_brewer(type="qual", palette=colorScaleNum) +
      geom_hline(y=90, linetype=2, size=1, color="grey85") +
      #   theme(axis.ticks.y = element_blank(), axis.text.y=element_blank(), 
      #         axis.title.y = element_blank()) +
      theme(strip.background = element_blank(), strip.text = element_blank(),
            panel.grid.major = element_line(colour = "grey85", size = 0.2),
            axis.title.x = element_text(vjust = 0),
            axis.title.y = element_text(vjust = 0.4, angle=90))
    
    print(g1)
    dev.off()
  }
  
  for (i in 1:8){
    tryColorScheme(i)
  }
  
}