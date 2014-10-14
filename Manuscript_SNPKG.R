
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

KGSNPs <- read.table("../Data/KG.map")[,-c(1,3)]
colnames(KGSNPs) <- c("snp.id", "snp.position")
KGSNPs$snp.position <- KGSNPs$snp.position/1000

KGSNPs.preQC <- read.table("../Data/1000Genome.bim")[,c(2,4)]
colnames(KGSNPs.preQC) <- c("snp.id", "snp.position")
KGSNPs.preQC$snp.position <- KGSNPs.preQC$snp.position/1000



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

#have a look
cbind("post-QC" = table(data_allele_0$allele), "pre-QC" = table(data_allele_0.pQC$allele))


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
               aes(x = snp.distance, y=..count.., alpha = "post-QC"), color = NA, fill = 'grey50', adjust = 0.05) + #  
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
print(gp$dist + xlim(-200,200))
print(gp$dist_cum + xlim(-200,200) + ylim(0,700))
print(gp$dist_cum_abs + xlim(0,200) + ylim(0,1300))
dev.off()


