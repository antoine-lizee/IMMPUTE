
source("Manuscript_palettes.R")

# Get the KG data ---------------------------------------------------------

library(plyr)

## Loci positions, from direct look-up in the Ucsc database
GenePosition_hg19 <- list("A" = c(29910247,29913661)/1000, 
                          "B" = c(31321649, 31324989)/1000, 
                          "C" = c(31236526, 31239913)/1000, 
                          "DRB1" = c(32546547, 32557613)/1000)
#                           "DQA1" = c(32605183, 32611429)/1000, 
#                           "DQB1" = c(32627241, 32634466)/1000,
#                           "DPB1" = c(33043703, 33057473)/1000)

## Launch the graph generating loop, per population
genesDf <- do.call(rbind, 
                   lapply(names(GenePosition_hg19), 
                          function(allelei) { dat <- GenePosition_hg19[[allelei]]
                                              within(data.frame(allele = allelei, 
                                                                x = mean(dat), 
                                                                xmin = dat[1], 
                                                                xmax = dat[2], 
                                                                stringsAsFactors=F), {
                                                                  xminD <- xmin - x
                                                                  xmaxD <- xmax - x
                                                                })}))

KGSNPs <- read.table("~/Dropbox/IMMPUTE/Side_Investigations/mhcKG_pre_QC/KG.map")[,-c(1,3)]
colnames(KGSNPs) <- c("snp.id", "snp.position")
KGSNPs$snp.position <- KGSNPs$snp.position/1000

KGSNPs.preQC <- read.table("~/Dropbox/IMMPUTE/Side_Investigations/mhcKG_pre_QC/1000Genome.bim")[,c(2,4)]
colnames(KGSNPs.preQC) <- c("snp.id", "snp.position")
KGSNPs.preQC$snp.position <- KGSNPs.preQC$snp.position/1000

KGSNPs.preMerge <- read.table("~/Dropbox/IMMPUTE/Side_Investigations/mhcKG_pre_QC/KG_MHCSnps.bim")[,c(2,4)]
colnames(KGSNPs.preMerge) <- c("snp.id", "snp.position")
KGSNPs.preMerge$snp.position <- KGSNPs.preMerge$snp.position/1000



# Prepare the data --------------------------------------------------------

isBtwn <- function(vect, range) {
  vect >= range[1] & vect <= range[2]
}

data_allele_0 <- do.call(rbind, 
                         lapply(names(GenePosition_hg19), 
                                function(allele_i) { 
                                  range_i <- GenePosition_hg19[[allele_i]] + c(-500, 500)
                                  within(data.frame(allele = allele_i, 
                                                    KGSNPs[isBtwn(KGSNPs$snp.position, range_i),]),
                                         snp.distance <- snp.position - genesDf[genesDf$allele == allele_i, "x"])
                                }))
data_allele <- ddply(data_allele_0, ~ allele, function(df) {
  df2 <- df[order(abs(df$snp.distance)),]
  df3 <- data.frame(df2, snp.order = 1:nrow(df2), 
                    cumsnp = (1:nrow(df2)),
                    cumsnpRate = (1:nrow(df2))/nrow(df2))
  df2bis <- df3[order(df3$snp.distance),]
  df3 <- within(df2bis,
                cumsnpRel <- abs(1:(nrow(df2bis))-which(df2bis$snp.distance>0)[1]))})

##Do the same for the preQC SNPs
data_allele_0.pQC <- do.call(rbind, 
                             lapply(names(GenePosition_hg19), 
                                    function(allele_i) { 
                                      range_i <- GenePosition_hg19[[allele_i]] + c(-500, 500)
                                      within(data.frame(allele = allele_i, 
                                                        KGSNPs.preQC[isBtwn(KGSNPs.preQC$snp.position, range_i),]),
                                             snp.distance <- snp.position - genesDf[genesDf$allele == allele_i, "x"])
                                    }))
data_allele.pQC <- ddply(data_allele_0.pQC, ~ allele, function(df) {
  df2 <- df[order(abs(df$snp.distance)),]
  df3 <- data.frame(df2, snp.order = 1:nrow(df2), 
                    cumsnp = (1:nrow(df2)),
                    cumsnpRate = (1:nrow(df2))/nrow(df2))
  df2bis <- df3[order(df3$snp.distance),]
  df3 <- within(df2bis,
                cumsnpRel <- abs(1:(nrow(df2bis))-which(df2bis$snp.distance>0)[1]))})

##Do the same for the preMerge SNPs
data_allele_0.pM <- do.call(rbind, 
                            lapply(names(GenePosition_hg19), 
                                   function(allele_i) { 
                                     range_i <- GenePosition_hg19[[allele_i]] + c(-500, 500)
                                     within(data.frame(allele = allele_i, 
                                                       KGSNPs.preMerge[isBtwn(KGSNPs.preMerge$snp.position, range_i),]),
                                            snp.distance <- snp.position - genesDf[genesDf$allele == allele_i, "x"])
                                   }))
data_allele.pM <- ddply(data_allele_0.pM, ~ allele, function(df) {
  df2 <- df[order(abs(df$snp.distance)),]
  df3 <- data.frame(df2, snp.order = 1:nrow(df2), 
                    cumsnp = (1:nrow(df2)),
                    cumsnpRate = (1:nrow(df2))/nrow(df2))
  df2bis <- df3[order(df3$snp.distance),]
  df3 <- within(df2bis,
                cumsnpRel <- abs(1:(nrow(df2bis))-which(df2bis$snp.distance>0)[1]))})


#have a look
cbind("post-QC" = table(data_allele_0$allele), "pre-QC" = table(data_allele_0.pQC$allele), "pre-Merge" = table(data_allele_0.pM$allele))


# Graph it ----------------------------------------------------------------
gp <- list()

gp$pos <- ggplot(data = data_allele_0) + 
  #   geom_density(data = data_allele_0,
  #                aes(x = snp.position), color = NA, fill = 'grey30', alpha=1,  adjust = 0.1) + #  
  geom_rect(data = genesDf, aes(xmin = xmin, xmax = xmax, fill = allele), alpha = 0.3, color = NA, ymin = 0, ymax = Inf-1) + # ymin = -1.2, ymax = -0.2)+ 
  geom_rect(data = genesDf, aes(xmin = xmin, xmax = xmax, fill = allele, color = allele), ymin = -1200, ymax = 200) + # 
  geom_density(aes(x = snp.position, y=..count.., fill = allele, color = allele), size = 0.5, alpha=0.3,  adjust = 0.1) + #  
  labs(title = paste("Comparison of SNP density around the four HLA loci") ,
       x = "SNP position (kbp)", y="SNP density", fill = "locus", color = "locus") +
  scale_y_continuous(breaks = NULL) +
  facet_grid(allele~.) +
  theme_bw() 

gp$dist <- ggplot(data = data_allele_0) + 
  geom_density(data = data_allele_0,
               aes(x = snp.distance, y=..count.., alpha = "post-QC"), color = NA, fill = 'grey50', alpha = 0.8, adjust = 0.05) + #  
  geom_rect(data = genesDf, aes(xmin = xminD, xmax = xmaxD, fill = allele), alpha = 0.3, color = NA, ymin = 0, ymax = Inf-1) + # ymin = -1.2, ymax = -0.2)+ 
  geom_rect(data = genesDf, aes(xmin = xminD, xmax = xmaxD, fill = allele, color = allele), ymin = -0.3, ymax = -0.1) + # 
  geom_density(data = data_allele_0.pQC, aes(x = snp.distance, y=..count.., fill = allele, color = allele, alpha = "pre-QC"), size = 0.5,  adjust = 0.05) + #  
  labs(title = paste("Comparison of SNP density around the four HLA loci, centered on each locus") ,
       x = "SNP distance to the locus (kbp)", y="SNP density (SNP/kb)", fill = "locus", color = "locus") +
  #   scale_y_continuous(breaks = NULL) +
  facet_grid(allele~.) +
  theme_bw() +
  scale_alpha_manual("", limits = c("pre-QC", "post-QC"), values = c(0.3, 1))

gp$dist_cum <- ggplot(data = data_allele) + 
  geom_rect(data = genesDf, aes(xmin = xminD, xmax = xmaxD, fill = allele), alpha = 0.3, color = NA, ymin = 0, ymax = Inf-1) + # ymin = -1.2, ymax = -0.2)+ 
  geom_rect(data = genesDf, aes(xmin = xminD, xmax = xmaxD, fill = allele, color = allele), ymin = -5, ymax = -1) + # 
  geom_line(aes(x = snp.distance, y=cumsnpRel, fill = allele, color = allele), size = 0.5) + #  
  geom_line(data = data_allele.pQC, aes(x = snp.distance, y=cumsnpRel, fill = allele, color = allele), size = 0.5, alpha = 0.4) + #  
  labs(title = paste("Comparison of SNP cumulative number around the four HLA loci, centered on each locus") ,
       x = "SNP distance to the locus (kbp)", y="SNP cumulative number", fill = "locus", color = "locus") +
  #   scale_y_continuous(breaks = NULL) +
  facet_grid(allele~.) +
  theme_bw() 
# gp$dist_cum

gp$dist_cum_abs <- ggplot(data = data_allele) + 
  geom_line(aes(x = abs(snp.distance), y=cumsnp, fill = allele, color = allele), size = 1) + #  
  geom_line(data = data_allele.pQC, aes(x = abs(snp.distance), y=cumsnp, fill = allele, color = allele), size = 0.5, alpha = 0.5) + #  
  labs(title = paste("Comparison of SNP cumulative number around the four HLA loci, depending on the distance") ,
       x = "SNP absolute distance to the locus (kbp)", y="SNP cumulative number", fill = "locus", color = "locus") +
  theme_bw() 
# gp$dist_cum_abs

pdf("/media/FD/Dropbox/IMMPUTE/Side_Investigations/Extra_Fig_SNP_DRB1.pdf", w = 12, h= 8)
dummy <- lapply(gp, print)
print(ggplot(data = data_allele_0) + 
        geom_density(aes(x = snp.distance, y=..count.., alpha = "post-QC"), color = NA, fill = 'grey50', alpha = 0.8, adjust = 0.02) + #  
        geom_rect(data = genesDf, aes(xmin = xminD, xmax = xmaxD, fill = allele), alpha = 0.3, color = NA, ymin = 0, ymax = Inf-1) + # ymin = -1.2, ymax = -0.2)+ 
        geom_rect(data = genesDf, aes(xmin = xminD, xmax = xmaxD, fill = allele, color = allele), ymin = -0.3, ymax = -0.1) + # 
        geom_density(data = data_allele_0.pQC, aes(x = snp.distance, y=..count.., fill = allele, color = allele, alpha = "pre-QC"), size = 0.5,  adjust = 0.02) + #  
        labs(title = paste("Comparison of SNP density around the four HLA loci, centered on each locus") ,
             x = "SNP distance to the locus (kbp)", y="SNP density (SNP/kb)", fill = "locus", color = "locus") +
        facet_grid(allele~.) +
        theme_bw() +
        scale_alpha_manual("", limits = c("pre-QC", "post-QC"), values = c(0.3, 1)) + 
        xlim(-200,200))
print(gp$dist_cum + xlim(-200,200) + ylim(0,700))
print(gp$dist_cum_abs + xlim(0,200) + ylim(0,1300))
dev.off()



# Final tries -------------------------------------------------------------

## Classic density plot / locus

pdf("/media/FD/Dropbox/IMMPUTE/Side_Investigations/Extra_Fig_SNP_DRB1_1s.pdf", w = 12, h= 8)
print(ggplot(data = data_allele_0) + 
        geom_density(data = data_allele_0.pM, aes(x = snp.distance, y=..count../10, alpha = "pre-Merge"), color = "grey30", fill = NA, size = 0.8,  adjust = 0.1) + #  
        geom_density(aes(x = snp.distance, y=..count.., alpha = "post-QC"), color = NA, fill = 'grey50', alpha = 0.8, adjust = 0.05) + #  
        geom_rect(data = genesDf, aes(xmin = xminD, xmax = xmaxD, fill = allele), alpha = 0.3, color = NA, ymin = 0, ymax = Inf-1) + # ymin = -1.2, ymax = -0.2)+ 
        geom_rect(data = genesDf, aes(xmin = xminD, xmax = xmaxD, fill = allele, color = allele), ymin = -0.3, ymax = -0.1) + # 
        geom_density(data = data_allele_0.pQC, aes(x = snp.distance, y=..count.., fill = allele, color = allele, alpha = "pre-QC"), size = 0.5,  adjust = 0.05) + #  
        labs(title = paste("Comparison of SNP density around the four HLA loci, centered on each locus") ,
             x = "SNP distance to the locus (kbp)", y="SNP density (SNP/kb)", fill = "locus", color = "locus") +
        #   scale_y_continuous(breaks = NULL) +
        facet_grid(allele~.) +
        scale_alpha_manual("", limits = c( "pre-Merge", "pre-QC", "post-QC"), values = c(1, 0.3, 1)) +
        scale_color_brewer(type="qual", palette=2) + 
        scale_fill_brewer(type="qual", palette=2) +
        theme_bw()
      #         theme(legend.key=element_rect(fill='pink'))
      #         guides()
)
dev.off()

# badass version with inset

b.bigger <- TRUE
data_locus <- data.frame(allele = unique(data_allele_0$allele), x_text = -485)

theme_perso_grid <- function(base_size, ...) {
  theme_bw(base_size, ...) %+replace% theme(
    strip.background =   element_blank(),
    panel.border =       element_blank(),
    axis.title.x = element_text(vjust = 0),
    axis.title.y = element_text(vjust = 1, angle = 90),
    strip.background = element_blank(), strip.text = element_blank(),
    legend.text = element_text(size = rel(0.9)),
    legend.key = element_rect(colour = NA),
    legend.key.height = unit(0.8, "lines"),
    legend.key.width = unit(0.8, "lines"),
    axis.ticks.y = element_line(size = 0.3),
    axis.ticks.margin = unit(0.1, "lines"),
    title = element_text(vjust = 1.2),
    legend.justification = c("left", "bottom"),
    legend.box = ifelse(b.bigger, "horizontal", "vertical"),
    legend.box.just = "left",
    legend.position = c(1.01, ifelse(b.bigger, 0.62, 0.5)),
    plot.margin = unit(c(1, ifelse(b.bigger, 8, 6), 0.5, 0.5), "lines") # default @ http://docs.ggplot2.org/dev/vignettes/themes.html
    
  )
}

gdens <- ggplot(data = data_allele_0) + 
  geom_density(data = data_allele_0.pM, aes(x = snp.distance, y=..count../10, alpha = "pre-Merge"), color = "grey30", fill = NA, size = 0.5,  adjust = 0.1) + #  
  geom_density(aes(x = snp.distance, y=..count.., alpha = "post-QC"), color = NA, fill = 'grey50', alpha = 0.8, adjust = 0.05) + #  
  geom_rect(data = genesDf, aes(xmin = xminD, xmax = xmaxD, fill = allele), alpha = 0.3, color = NA, ymin = 0, ymax = Inf-1) + # ymin = -1.2, ymax = -0.2)+ 
  geom_rect(data = genesDf, aes(xmin = xminD, xmax = xmaxD, fill = allele, color = allele), ymin = -0.3, ymax = -0.1) + # 
  geom_density(data = data_allele_0.pQC, aes(x = snp.distance, y=..count.., fill = allele, color = allele, alpha = "pre-QC"), size = 0.3,  adjust = 0.05) + #  
  labs(title = paste("Comparison of SNP density around the four HLA loci, centered on each locus") ,
       x = "SNP distance to the locus (kbp)", y="SNP density (SNP/kb)", fill = "Locus", color = "Locus", alpha = "Status \nof the SNPs") +
  #   scale_y_continuous(breaks = NULL) +
  facet_grid(allele~.) +
  scale_alpha_manual(limits = c( "pre-Merge", "pre-QC", "post-QC"), values = c(1, 0.3, 1)) +
  scale_color_brewer(type="qual", palette=2) + 
  scale_fill_brewer(type="qual", palette=2) +
  theme_perso_grid(base_size = 8) + 
  geom_text(data = data_locus, aes(label = allele, x = x_text), y= 11, size = 4, hjust = 0) + 
  scale_x_continuous(expand = c(0.02,0.02)) + 
  scale_y_continuous(breaks = c(0,5,10)) +
  guides(alpha = guide_legend(order = 1), fill = guide_legend(order = 2), color = guide_legend(order = 2))

ginset <- ggplot(data = data_allele) + 
  geom_line(aes(x = abs(snp.distance), y=cumsnp, fill = allele, color = allele), size = 0.8) + #  
  #         geom_line(data = data_allele.pQC, aes(x = abs(snp.distance), y=cumsnp, fill = allele, color = allele), size = 0.5, alpha = 0.5) + #  
  (if(b.bigger) {labs(x = "absolute distance to the locus (kbp)", y="cumulative count of SNPs", fill = "locus", color = "locus")
  } else { labs(x = NULL, y=NULL, fill = NULL, color = NULL)})  +
  theme_bw(base_size = 6)  + xlim(0,200) + ylim(0,1300) +
  scale_color_brewer(guide = FALSE, type="qual", palette=2) + 
  theme(
    axis.text.y = element_text(angle = 90, hjust = 0.5),
    plot.margin = unit(c(0.2,0.2,0.2,0.2), "lines"))

pdf("/media/FD/Dropbox/IMMPUTE/Manuscript/Fig_SNPs.pdf", w=cmToInches(18.3), h=cmToInches(12))
print(gdens)
#A viewport taking up a fraction of the plot area
if (b.bigger) {
  vp <- viewport(width = 0.23, height = 0.33, x = 0.750, y = 0.13, just = c("left", "bottom"))
} else {
  vp <- viewport(width = 0.18, height = 0.27, x = 0.815, y = 0.15, just = c("left", "bottom"))
}
print(ginset,
      vp = vp)
dev.off()



## Simple plain version ofthe inset

pdf("/media/FD/Dropbox/IMMPUTE/Side_Investigations/Extra_Fig_SNP_DRB1_2s.pdf", w = 5, h= 4)
print(ggplot(data = data_allele) + 
        geom_line(aes(x = abs(snp.distance), y=cumsnp, fill = allele, color = allele), size = 1) + #  
        #         geom_line(data = data_allele.pQC, aes(x = abs(snp.distance), y=cumsnp, fill = allele, color = allele), size = 0.5, alpha = 0.5) + #  
        labs(title = paste("Comparison of cumulative count of SNPs around the four HLA loci") ,
             x = "SNP absolute distance to the locus (kbp)", y="cumulative count of SNPs", fill = "locus", color = "locus") +
        theme_bw()  + xlim(0,200) + ylim(0,1300) +
        scale_color_brewer(type="qual", palette=2))
dev.off()
