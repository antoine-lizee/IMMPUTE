## Script creating the figures 1 & 3 for the IMMPUTE paper
#
# Copyright 08/11/2014 Antoine Lizee @ UCSF.
# antoine.lizee@ucsf.edu

options(java.parameters="-Xmx4g")
require(XLConnect)
require(grid)
library(ggplot2)
library(reshape2)
source("geom_point_width_modified.R")

source("Manuscript_palettes.R")

# Preapre the data --------------------------------------------------------


methods <- c("HIBAG", "MAGPrediction", "e-HLA", "HLA*IMP:02")
loci <- c("A", "B", "C", "DRB1")

data <- data.matrix(readWorksheetFromFile("/media/FD/Dropbox/IMMPUTE/Manuscript/Figure 1.xlsx", sheet=1, 
                                          startRow=3, endRow=7, startCol=4, endCol=12, #Not used
                                          region="D3:L7", header=T, rownames=1))

data2 <- array(c(data[,1:4 *2-1], data[,1:4 *2]), c(4,4,2))
dimnames(data2) <- list(loci, methods, c("cr100","cr80"))

dataDf <- melt(data2, varnames=c("Locus", "Method", "CallRate"), value.name="Performance")
dataDf2 <- dcast(dataDf, Locus + Method ~ CallRate )



# Create the helper functions ---------------------------------------------


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

theme_perso_full <- theme_bw(base_size=16) %+replace% theme_perso
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

exportFigure1Dummy <- function (name, ggp) {
  png(paste0(name,".png"),  w=8.9, h=7.2, units="cm", res=300)
  print(ggp + theme_perso)
  dev.off()
}

# plotFigure1 <- function (name, dataDfi, export = T) {
#   if (export) png(paste0("/media/FD/Dropbox/IMMPUTE/Manuscript/", name,".png"), w=620, h=500)
#   require(ggplot2)
#   g1 <- ggplot(dataDfi) + theme_perso +
#     geom_abline(a=1, b=0, linetype=2, size=1.2, color="grey85") +
#     geom_point(aes(y=cr80*100, x=cr100*100, shape=Locus, color=Method), size=5) +
#     scale_shape_manual(values=c(15,19,18,17)) +
#     labs(y="Performance with 80% call rate (%)", x="Performance with 100% call rate(%)") +
#     ylim(60, 100) + xlim(60,100) +
#     coord_equal() +
#     scale_color_brewer(type="qual", palette=3)
#   
#   print(g1)
#   if (export) dev.off()
# }



# Create and export the charts --------------------------------------------

g1 <- ggplot(dataDf2) + theme_bw(base_size=16) +
  geom_abline(a=1, b=0, linetype=2, size=0.6, color="grey80") +
#   geom_point(aes(y=Mask*100, x=None*100, shape=Locus, color=Method), size=2.5) +
  geom_point(aes(y=cr80*100, x=cr100*100, shape=Locus, fill=Method), colour = "black", size=2.2, width = 0.4, alpha = 0.7) +
#   scale_shape_manual(values=c(15,19,18,17)) +
  scale_shape_manual(values=c(21,22,23,24)) +
  guides(fill = guide_legend(override.aes = list(shape = 25))) +
  labs(y="Imputation Accuracy with 80% call rate (%)", x="Imputation Accuracy with 100% call rate (%)") +
  ylim(60, 100) + xlim(60,100) +
  coord_equal() +
  palette_perso_fill2 +
  palette_perso2
  #     scale_color_brewer(type="qual", palette=2)

# exportFigure1Dummy("Figure_1", g1)

print(g1)
exportFigure1("Figure_1", g1)


# Try different color schemes ---------------------------------------------

if (colorSchemes.b <- FALSE) {
  
  tryColorScheme <- function(colorScaleNum){
    gi <- ggplot(dataDf2) + theme_bw(base_size=16) +
      geom_abline(a=1, b=0, linetype=2, size=0.6, color="grey80") +
      geom_point(aes(y=cr80*100, x=cr100*100, shape=Locus, fill=Method), colour = "black", size=2.2, width = 0.4, alpha = 0.75) +
      scale_shape_manual(values=c(21,22,23,24)) +
      guides(colour = guide_legend(override.aes = list(shape = 25))) +
      labs(y="Imputation Accuracy with 80% call rate (%)", x="Imputation Accuracy with 100% call rate(%)") +
      ylim(60, 100) + xlim(60,100) +
      coord_equal() +
      palette_perso_fill2 +
      palette_perso2
      #     scale_color_brewer(type="qual", palette=colorScaleNum)
    
    exportFigure1(paste0("colorSchemes/Figure_1_", colorScaleNum), gi)
  }
  
  for (i in 1:8){
    tryColorScheme(i)
  }
  
}


# Figure 3 ----------------------------------------------------------------
# (Very Similar)

if (figure_3_boolean <- T){
  
  data <- data.matrix(readWorksheetFromFile("/media/FD/Dropbox/IMMPUTE/Manuscript/Figure 3.xlsx", sheet=1, 
                                            region="D3:L7", header=T, rownames=1))
  
  data2 <- array(c(data[,1:4 *2-1], data[,1:4 *2]), c(4,4,2))
  dimnames(data2) <- list(loci, methods, c("None","Mask"))
  
  dataDf <- melt(data2, varnames=c("Locus", "Method", "Masking"), value.name="Performance")
  dataDf2 <- dcast(dataDf, Locus + Method ~ Masking )
  
  require(ggplot2)
  g2 <- ggplot(dataDf2) + theme_bw(base_size=16) +
    geom_abline(a=1, b=0, linetype=2, size=0.6, color="grey80") +
    geom_point(aes(y=Mask*100, x=None*100, shape=Locus, fill=Method), colour = "black", size=2.2, width = 0.4, alpha = 0.75) +
    scale_shape_manual(values=c(21,22,23,24)) +
    guides(fill = guide_legend(override.aes = list(shape = 25))) +
    labs(y="Imputation Accuracy - Mask Non-Equivalence (%)", x="Imputation Accuracy - No Mask (%)") +
    ylim(60, 100) + xlim(60,100) +
    coord_equal() +
    palette_perso_fill2 +
    palette_perso2
    #     scale_color_brewer(type="qual", palette=2)
  
  print(g2)
  exportFigure1("Figure_3", g2)
}

