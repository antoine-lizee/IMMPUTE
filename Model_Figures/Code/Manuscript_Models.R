
require(reshape2)
require(plyr)
require(aod)
require(ggplot2)
require(grid)

rm(list = ls())

source("Model_Figures/Code/Helpers/Manuscript_palettes.R")
source("Model_Figures/Code/Helpers/Manuscript_Utilities.R")
outputFolder <- "/media/FD/Dropbox/IMMPUTE/Manuscript/Data Tables and Figures/All Figures/"

# Import MATCH data from dropbox -------------------------------------------

methods <- c("BROAD", "CCHMC", "FHCRC", "MCRI", "UWASH")
loci <- c("A", "B", "C", "DRB1")
# Find by hand the correct paths
parentDir <- "~/Dropbox/IMMPUTE/Imputation Scoring"
dir(pattern = "MATCH_out.txt", path = parentDir, full.names = T, recursive = T)
filePaths <- matrix(c("~/Dropbox/IMMPUTE/Imputation Scoring/CCHMC/SVM AI/Experiment 2/Score_2F_NoMask/MATCH_out.txt", "CCHMC",
                      "~/Dropbox/IMMPUTE/Imputation Scoring/FHCRC/Experiment 2/Score_2F_NoMask/MATCH_out.txt", "FHCRC",
                      "~/Dropbox/IMMPUTE/Imputation Scoring/MCRI/Experiment 2/Score_2F_NoMask/MATCH_out.txt", "MCRI",
                      "~/Dropbox/IMMPUTE/Imputation Scoring/MCRI/Experiment 2 update fast/Score_2F_NoMask/MATCH_out.txt", "MCRI2fast",
                      "~/Dropbox/IMMPUTE/Imputation Scoring/MCRI/Experiment 2 update standard/Score_2F_NoMask/MATCH_out.txt", "MCRI2std",
                      "~/Dropbox/IMMPUTE/Imputation Scoring/UWASH/Experiment 2/Score_2F_NoMask/MATCH_out.txt", "UWASH",
                      "~/Dropbox/IMMPUTE/Imputation Scoring/UWASH/Experiment 2 update/Score_2F_NoMask/MATCH_out.txt", "UWASH2"),
                    ncol = 2, byrow = T)
#"~/Dropbox/IMMPUTE/Imputation Scoring/BROAD/Official/4D/PGM conservative Mask_None/MATCH_out.txt", "BROAD",

# Import data
matchData <- NULL
for (i in 1:nrow(filePaths)) {
  res <- read.delim(filePaths[i,1])
  res2 <- with(res, data.frame(sample.id = as.factor(SampleID),
                               A = A1 + A2,
                               B = B1 + B2,
                               C = C1 + C2,
                               DRB1 = DRB1.1 + DRB1.2)
  )
  matchData <- rbind(matchData, data.frame(melt(res2, value.name = "success", id.vars = "sample.id", variable.name = "locus"), method = filePaths[i,2]))
}

##Apply LUT
institutes <- c("UWASH", "UWASH2", "FHCRC", "CCHMC", "MCRI", "MCRI2fast", "MCRI2std")
methods <- c("HIBAG", "HIBAG2", "MAGPrediction", "e-HLA", "HLA*IMP:02Leg", "HLA*IMP:02Fast", "HLA*IMP:02Std")
levels(matchData$method)[levels(matchData$method) %in% institutes] <- methods[match(levels(matchData$method), institutes, nomatch=0)]

# Quick look at the raw data
gp <- ggplot(matchData) + geom_bar(aes(x=success)) + facet_grid(method ~ locus) + theme_bw()


# Preprocessing ------------------------------------------

# Prepare data
matchData2 <- ddply(matchData, ~ locus + method, summarize, success = sum(success, na.rm = T), fail = length(sample.id)*2 - success)
within(matchData2, sumsf <- success + fail) # just to check

# Inverse-link function:
plogit <- function(t) {
  return(1/(1+exp(-t)))
}

# Filter out the lesser versions for every model
matchData3 <- matchData2[ ! (matchData2$method %in% c("HIBAG", "HLA*IMP:02Leg", "HLA*IMP:02Std")), ]
matchData3$method <- factor(matchData3$method) 
levels(matchData3$method)[levels(matchData3$method)=="HLA*IMP:02Fast"] <- "HLA*IMP:02"
levels(matchData3$method)[levels(matchData3$method)=="HIBAGv1.3"] <- "HIBAG"

# Create European data from the table
pct <- c(96.1, 91.7, 89.8, 93.7, 
         89.3, 86.9, 84.0, 83.0, 
         95.6, 95.1, 90.3, 96.6, 
         59.2, 52.9, 53.4, 53.4)
methods <- c( "HIBAG", "MAGPrediction", "e-HLA", "HLA*IMP:02")
matchData3Euros <- within(matchData3, {
  method <- rep(methods, 4)
  success <- round(pct*1.03*2)
  fail <- round((100-pct)*1.03*2)}
)

# Change the representation of the data for tweaking the form of the models
matchData4 <- booleaniseFactors(matchData3, c("method", "locus"))
matchData4 <- swapColumns(matchData4, "DRB", "locus.A")
matchData4 <- swapColumns(matchData4, "e-HLA", "HIBAG")
matchData4Euros <- booleaniseFactors(matchData3Euros, c("method", "locus"))
matchData4Euros <- swapColumns(matchData4Euros, "DRB", "locus.A")
matchData4Euros <- swapColumns(matchData4Euros, "e-HLA", "HLA*IMP")


# Creation of the models --------------------------------------------------

MOD <- glm(data = matchData4, cbind(success, fail) ~ . + 0, family = binomial(logit))
MOD.intercept <- glm(data = matchData4, cbind(success, fail) ~ ., family = binomial(logit))
MOD.grouped <- glm(data = matchData4, cbind(success, fail) ~ locus.B + locus.C + locus.DRB1 + locus.A + method.HIBAG, family = binomial(logit))

# Strictly similar, just a different method for verification
mod <- glm(data = matchData3, cbind(success, fail) ~ method + locus + 0, family = binomial(logit))
mod.intercept <- glm(data = matchData3, cbind(success, fail) ~ method + locus, family = binomial(logit))
mod.grouped <- glm(data = matchData3, cbind(success, fail) ~ locus + (method=="HIBAG"), family = binomial(logit))

# Europeans
MODEURO <- glm(data = matchData4Euros, cbind(success, fail) ~ . + 0, family = binomial(logit))
MODEURO.intercept <- glm(data = matchData4Euros, cbind(success, fail) ~ ., family = binomial(logit))
modEuro.grouped <- glm(data = matchData3Euros, cbind(success, fail) ~ locus + (method=="HIBAG"), family = binomial(logit))
modEuro.intercept <- glm(data = matchData3Euros, cbind(success, fail) ~ method + locus, family = binomial(logit))

MODEURO.grouped <- glm(data = matchData4Euros, cbind(success, fail) ~ locus.B+ locus.DRB1 + (locus.C|locus.A)  + method.HIBAG + (`method.HLA*IMP:02`|`method.e-HLA`|method.MAGPrediction), family = binomial(logit))
names(MODEURO.grouped$coefficients)[c(4,6)] <- c("locus.AorC", "method.NON-HIBAG")

# With duplicates
modAll.intercept <- glm(data = matchData2, cbind(success, fail) ~ method + locus, family = binomial(logit))


# Significance testing ----------------------------------------------------
# Wald Tests mainly

## Base model, full populations
# Methods 
wald.test(b = coef(mod.intercept), Sigma = vcov(mod.intercept), Terms = 2:4, verbose = T)
# Loci
wald.test(b = coef(mod.intercept), Sigma = vcov(mod.intercept), Terms = 5:7)

# Loci, separately, hence in comparison to locus A
ltot1 <- list(list(hyp = "B = A", term = 5),
              list(hyp = "C = A", term = 6),
              list(hyp = "DRB1 = A", term = 7))
for (l in ltot1) {
  cat("## Testing hypothesis: \n ### ", l[["hyp"]], "### \n")
  w <- wald.test(b = coef(mod.intercept), Sigma = vcov(mod.intercept), Term = l[["term"]])
  printWaldTest(w)
  cat("\n")
}

# Methods, separately, hence in comparison to e-HLA
ltot1 <- list(list(hyp = "MAGPrediction = e-HLA", term = 2),
              list(hyp = "HLA*IMP:02 = e-HLA", term = 3),
              list(hyp = "HIBAG = e-HLA", term = 4))
for (l in ltot1) {
  cat("## Testing hypothesis: \n ### ", l[["hyp"]], "### \n")
  w <- wald.test(b = coef(mod.intercept), Sigma = vcov(mod.intercept), Term = l[["term"]])
  printWaldTest(w)
  cat("\n")
}

## Full Array, using the helpers
waldLocus <- sapply(5:6, function(i) sapply((i+1):7, function(j) printWaldSignif(mod.intercept, i, j)))
waldMethods <- sapply(2:3, function(i) sapply((i+1):4, function(j) printWaldSignif(mod.intercept, i, j)))
# For Europeans
waldLocus <- sapply(5:6, function(i) sapply((i+1):7, function(j) printWaldSignif(modEuro.intercept, i, j)))
waldMethods <- sapply(2:3, function(i) sapply((i+1):4, function(j) printWaldSignif(modEuro.intercept, i, j)))
# With all methods
waldLocus <- sapply(8:9, function(i) sapply((i+1):10, function(j) printWaldSignif(modAll.intercept, i, j)))
waldMethods <- sapply(2:6, function(i) sapply((i+1):7, function(j) printWaldSignif(modAll.intercept, i, j)))

cbind(unlist(waldMethods))


# SS ------------------------------------------------------------------

pct <- c(87.0, 88.0, 84.0, 86.0,
         58.5, 69.5, 65.5, 64.0,
         82.5, 94.0, 93.5, 92.0,
         42.5, 55.5, 49.5, 53.0)
methods <- c( "e-HLA", "HIBAG",  "HLA*IMP:02", "MAGPrediction")

matchData3SSA <- within(matchData3, {
  method <- rep(methods, 4)
  success <- round(pct*1.0*2)
  fail <- round((100-pct)*1.0*2)}
)

mod <- glm(data = matchData3SSA, cbind(success, fail) ~ locus + method, family = binomial(link=logit))
summary(mod)
confint(mod)
cbind(within(matchData3SSA, prop <- success / (success + fail)), p= plogit(predict(mod)), predict(mod), residuals(mod))

matchData4SSA <- booleaniseFactors(matchData3SSA, c("method", "locus"))
matchData4SSA <- matchData4SSA[c(1:2, 4:6, 3, 8:10, 7)]
MODSSA.intercept <- glm(data = matchData4SSA, cbind(success, fail) ~ ., family = binomial(logit))


# Random ------------------------------------------------------------------

sampleIds <- unique(matchData$sample.id)
matchData1R <- matchData[ matchData$sample.id %in% sampleIds[sample(length(sampleIds), 103)],]
matchData2R <- ddply(matchData1R, ~ locus + method, summarize, success = sum(success, na.rm = T), fail = length(sample.id)*2 - success)
within(matchData2R, sumsf <- success + fail) 

matchData3R <- matchData2R[ ! (matchData2$method %in% c("HIBAG2", "HLA*IMP:02Leg", "HLA*IMP:02Std")), ]
matchData3R$method <- factor(matchData3R$method)
levels(matchData3R$method)[levels(matchData3R$method)=="HLA*IMP:02Fast"] <- "HLA*IMP:02"

mod <- glm(data = matchData3R, cbind(success, fail) ~ locus + method, family = binomial(link=logit))
summary(mod)
confint(mod)
cbind(within(matchData3Euros, prop <- success / (success + fail)), p= plogit(predict(mod)), predict(mod), residuals(mod))

matchData4R <- booleaniseFactors(matchData3R, c("method", "locus"))
matchData4R <- matchData4R[c(1:2, 4:6, 3, 8:10, 7)]
MODR.intercept <- glm(data = matchData4R, cbind(success, fail) ~ ., family = binomial(logit))


# Plots -------------------------------------------------------------------

printGGplot(file=paste0(outputFolder, "Figure_Model_All"),
            plot=plotModelColors(MOD.intercept),
            w=8.9*2, h=7.2*2, units="cm", res=300/2)

printGGplot(file=paste0(outputFolder, "Figure_Model_Euro"), 
            plot=plotModelColors(MODEURO.intercept),
            w=8.9*2, h=7.2*2, units="cm", res=300/2)

ggmod <- ggplot_build(plotModelColors(MODEURO.intercept))
ylims <- c( min(ggmod$data[[2]]$ymin), max(ggmod$data[[2]]$ymax))
printGGplot(file=paste0(outputFolder, "Figure_Model_All_scaled"),
            plot=plotModelColors(MOD.intercept, ylims),
            w=8.9*2, h=7.2*2, units="cm", res=300/2)

ggsave(filename=paste0("/media/FD/Dropbox/IMMPUTE/Modeling/", "SupplFigure_ModelEuroReduced",".png"), 
       plot=plotModelColors(MODEURO.grouped),
       w=8.9*2, h=7.2*2, units="cm", dpi=300/2)

ggsave(filename=paste0("/media/FD/Dropbox/IMMPUTE/Modeling/", "SupplFigure_ModelSSaf",".png"), 
       plot=plotModelColors(MODSSA.intercept),
       w=8.9*2, h=7.2*2, units="cm", dpi=300/2)

ggsave(filename=paste0("/media/FD/Dropbox/IMMPUTE/Modeling/", "SupplFigure_ModelRandomSubset",".png"), 
       plot=plotModelColors(MODR.intercept),
       w=8.9*2, h=7.2*2, units="cm", dpi=300/2)


# Output model coefficients for the supplementary table --------------------

getCoeffs <- function(model) {
  cbind(round(coefficients(model),2), round(confint(model),2))
}

print(getCoeffs(MOD.intercept))
print(getCoeffs(MODEURO.intercept))


# Comparison of the different flavours within a method --------------------

## HIBAG
matchDataHB <- matchData2[matchData2$method %in% c("HIBAG", "HIBAG2"), ]
matchDataHB$method <- factor(matchDataHB$method)

matchData4HB <- booleaniseFactors(matchDataHB, c("method", "locus"))
matchData4HB <- matchData4HB[c(1:2, 4:3, 6:8, 5)]

MODHibags <- glm(data = matchData4HB, cbind(success, fail) ~ ., family = binomial(logit))
modHibags <- glm(data = matchDataHB, cbind(success, fail) ~ locus + method,  family = binomial(logit))
wald.test(b = coef(modHibags), Sigma = vcov(modHibags), Terms = 5)

## HLA*IMP:02
matchDataHLAIMP <- matchData2[matchData2$method %in% c("HLA*IMP:02Leg", "HLA*IMP:02Fast", "HLA*IMP:02Std"), ]
matchDataHLAIMP$method <- factor(matchDataHLAIMP$method)

matchData4HLAIMP <- booleaniseFactors(matchDataHLAIMP, c("method", "locus"))
matchData4HLAIMP <- matchData4HLAIMP[c(1:2, 4:5, 3, 7:9, 6)]

MODHLAIMP <- glm(data = matchData4HLAIMP, cbind(success, fail) ~ ., family = binomial(logit))
modHLAIMP <- glm(data = matchDataHLAIMP, cbind(success, fail) ~ locus + method,  family = binomial(logit))
wald.test(b = coef(modHLAIMP), Sigma = vcov(modHLAIMP), Terms = 5:6) # Significant
wald.test(b = coef(modHLAIMP), Sigma = vcov(modHLAIMP), Terms = 5)
wald.test(b = coef(modHLAIMP), Sigma = vcov(modHLAIMP), Terms = 6)


# Merge coefficients
coeffDf <- rbind(getCoeffDf(MODHibags), getCoeffDf(MODHLAIMP))
coeffDf <- coeffDf[coeffDf$type == "method",]
coeffDf$type <- c("HIBAG", "HIBAG", "HLA*IMP:02", "HLA*IMP:02", "HLA*IMP:02")
coeffDf$value <- relevel(coeffDf$value, "HLA*IMP:02Leg")

g.intercept.line <- geom_hline(y=1, color = "red", size = 0.8)
g.scale <- geom_blank()
ggMethods <- ggplot(coeffDf, aes(x=value, y=OR, ymin = OR_2.5, ymax = OR_97.5)) + theme_bw(base_size = 12) +
  g.intercept.line +
  geom_errorbar() + geom_point(aes(fill = value, shape = value), size = 4) +
  scale_fill_manual(values = rep("grey95", 5), guide = "none") +
  scale_shape_manual(values = c(21,21,22,22, 23), guide = "none")+
  facet_grid(~ type, scale = "free_x") + labs(x=NULL) +
  g.scale +
  theme(axis.text.x = element_text(angle = -12, hjust = 0.3),
        plot.margin = unit(c(0.5, 1, 0., 0.), "lines"),
        #           axis.title.x = element_blank(),
        panel.margin = unit(0.5, "lines"))

printGGplot(file=paste0(outputFolder, "Figure_Model_NewVsOld"),
            plot=ggMethods,
            w=8.9*2, h=7.2*2, units="cm", res=300/2)