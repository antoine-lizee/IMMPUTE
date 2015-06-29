
require(reshape2)
require(plyr)
require(aod)
require(ggplot2)

rm(list = ls())

source("Model_Figures/Code/helpers/Manuscript_palettes.R")
source("Model_Figures/Code/helpers/Manuscript_Utilities.R")
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
methods <- c("HIBAG", "HIBAG2", "MAGPrediction", "e-HLA", "HLA*IMP:02", "HLA*IMP:02Fast", "HLA*IMP:02Std")
##Method slow
# for (i in 1:length(institutes) ) {
#   institute <- intitutes[i]
#   institute.b <- matchData$method == institute 
#   matchData$method[institute.b] <- methods[i]
# }
##Method smarter (not working with factors)
# matchData$method[matchData$method %in% institutes] <- methods[match(matchData$method, institutes, nomatch=0)]
##Method simple, factor based -> change the levels
levels(matchData$method)[levels(matchData$method) %in% institutes] <- methods[match(levels(matchData$method), institutes, nomatch=0)]

# Quick look at the raw data
gp <- ggplot(matchData) + geom_bar(aes(x=success)) + facet_grid(method ~ locus) + theme_bw()



# Preprocessing and first models ------------------------------------------

# Prepare data
matchData2 <- ddply(matchData, ~ locus + method, summarize, success = sum(success, na.rm = T), fail = length(sample.id)*2 - sum(success))
within(matchData2, sumsf <- success + fail) # just to check

# First Quick models:
mod1 <- glm(data = matchData, cbind(success, 2-success) ~ locus + method, family = binomial(link=logit))
mod2 <- glm(data = matchData2, cbind(success, fail) ~ locus + method, family = binomial(link=logit)) # Slightly different than mod1, because each allele is supposed independent here.
# mod3 <- glm(data = matchData, cbind(success, 2-success) ~ ., family = binomial(link=logit)) # Very slow of course and resulting in more spread out coefficients for the ones which are not sample based

# Inverse-link function:
plogit <- function(t) {
  return(1/(1+exp(-t)))
}

matchData3 <- matchData2[ ! (matchData2$method %in% c("HIBAG2", "HLA*IMP:02", "HLA*IMP:02Std")), ]
matchData3$method <- factor(matchData3$method) 

mod <- glm(data = matchData3, cbind(success, fail) ~ locus + method, family = binomial(link=logit)) # Slightly different than mod1, because each allele is supposed independent here.
summary(mod)
confint(mod)
cbind(within(matchData3, prop <- success / (success + fail)), p= plogit(predict(mod)), predict(mod), residuals(mod))
anova(mod)

## Wald Tests : both are significant
# loci 
wald.test(b = coef(mod), Sigma = vcov(mod), Terms = 2:4, verbose = T)
# Methods
wald.test(b = coef(mod), Sigma = vcov(mod), Terms = 5:7)

## More subtle tests
printWaldTest <- function(w, digits = 2) {
  v <- w[["result"]][["chi2"]]
  cat("X2 = ", format(v["chi2"], digits = digits, nsmall = 1), 
      ", df = ", v["df"], ", P(> X2) = ", format(v["P"], digits = digits, 
                                                 nsmall = 1), "\n", sep = "")
}
# Methods, separately
ltot1 <- list(list(hyp = "MAGPrediction = e-HLA", term = 5),
              list(hyp = "HLA*IMP:02 = e-HLA", term = 6),
              list(hyp = "HIBAG = e-HLA", term = 7))
for (l in ltot1) {
  cat("## Testing hypothesis: \n ### ", l[["hyp"]], "### \n")
  w <- wald.test(b = coef(mod), Sigma = vcov(mod), Term = l[["term"]])
  printWaldTest(w)
  cat("\n")
}

# Methods, one against the others
ltot2 <- list(list(hyp = "MAGPrediction = HLA*IMP:02", coeffs = cbind(0, 0, 0, 0, 1, -1, 0)),
              list(hyp = "MAGPrediction = HIBAG", coeffs = cbind(0, 0, 0, 0, 1, 0, -1)),
              list(hyp = "HLA*IMP:02 = HIBAG", coeffs = cbind(0, 0, 0, 0, 0, 1, -1)),
              list(hyp = "locus B = locus DRB1", coeffs = cbind(0, 1, 0, -1, 0, 0, 0)),
              list(hyp = "locus A = locus B", coeffs = cbind(1, -1, 0, 0, 0, 0, 0)),
              list(hyp = "locus C = locus A", coeffs = cbind(-1, 0, 1, 0, 0, 0, 0)))
for (l in ltot2) {
  cat("## Testing hypothesis: \n ### ", l[["hyp"]], "### \n")
  w <- wald.test(b = coef(mod), Sigma = vcov(mod), L = l[["coeffs"]])
  printWaldTest(w,5)
  cat("\n")
}


# Europeans ---------------------------------------------------------------

pct <- c(96.1, 91.7, 89.8, 94.2, 
         89.3, 86.9, 84.0, 75.7, 
         95.6, 95.1, 90.3, 95.1, 
         58.3, 52.9, 53.4, 53.4)
methods <- c( "HIBAG", "MAGPrediction", "e-HLA", "HLA*IMP:02")

matchData3Euros <- within(matchData3, {
  method <- rep(methods, 4)
  success <- round(pct*1.03*2)
  fail <- round((100-pct)*1.03*2)}
)

mod <- glm(data = matchData3Euros, cbind(success, fail) ~ locus + method, family = binomial(link=logit))
summary(mod)
confint(mod)
cbind(within(matchData3Euros, prop <- success / (success + fail)), p= plogit(predict(mod)), predict(mod), residuals(mod))

anova(mod)
## Wald Tests
# loci 
wald.test(b = coef(mod), Sigma = vcov(mod), Terms = 2:4, verbose = T)
# Methods
wald.test(b = coef(mod), Sigma = vcov(mod), Terms = 5:7)

# Methods, separately
for (l in ltot1) {
  cat("## Testing hypothesis: \n ### ", l[["hyp"]], "### \n")
  w <- wald.test(b = coef(mod), Sigma = vcov(mod), Term = l[["term"]])
  printWaldTest(w)
  cat("\n")
}

# Methods, one against the others
for (l in ltot2) {
  cat("## Testing hypothesis: \n ### ", l[["hyp"]], "### \n")
  w <- wald.test(b = coef(mod), Sigma = vcov(mod), L = l[["coeffs"]])
  printWaldTest(w,4)
  cat("\n")
}


# Sum it up ---------------------------------------------------------------

# Change the data representation

booleaniseFactors <- function(df, factorNames) {
  for (cn in factorNames) {
    for (v in unique(c <- df[[cn]])) {
      df[paste0(cn, ".", v)] <- as.numeric((c == v))
    }
    df[cn] <- NULL
  }
  return(df)
}

swapColumns <- function(df, colName1, colName2, fixed = T) {
  ind1 <- grep(colName1, colnames(df), fixed = fixed)
  ind2 <- grep(colName2, colnames(df), fixed = fixed)
  stopifnot((length(ind1) == 1) & (length(ind2) == 1))
  if (ind1 < ind2) {
    indi <- ind2
    ind2 <- ind1
    ind1 <- indi
  }
  dfstart.b <- ind2 != 1
  dfend.b <- ind1 != (nc <- ncol(df))
  if (abs(ind1-ind2) == 1) {
    return(df[c(1:(ind2-1), ind1, ind2, (ind1+1):ncol(df))])
  }
  dfMiddle <- data.frame(df[ind1], df[(ind2+1):(ind1-1)], df[ind2],check.names = F)
  if (!dfstart.b) {
    if (!dfend.b) {
      return(dfMiddle)
    } else {
      dfEnd <- data.frame(df[(ind1+1):nc],check.names = F)
      return( data.frame(dfMiddle, dfEnd,check.names = F))
    }
  } else {
    dfStart <- data.frame(df[1:(ind2-1)],check.names = F)
    if (!dfend.b) {
      return(data.frame(dfStart, dfMiddle,check.names = F))
    } else {
      dfEnd <- data.frame(df[(ind1+1):nc],check.names = F)
      return(data.frame(dfStart, dfMiddle, dfEnd,check.names = F))
    }
  }
}

matchData4 <- booleaniseFactors(matchData3, c("method", "locus"))
matchData4 <- swapColumns(matchData4, "DRB", "locus.A")
matchData4 <- swapColumns(matchData4, "e-HLA", "HIBAG")
matchData4Euros <- booleaniseFactors(matchData3Euros, c("method", "locus"))
matchData4Euros <- swapColumns(matchData4Euros, "DRB", "locus.A")
matchData4Euros <- swapColumns(matchData4Euros, "e-HLA", "HLA*IMP")

MOD <- glm(data = matchData4, cbind(success, fail) ~ . + 0, family = binomial(logit))
MOD.intercept <- glm(data = matchData4, cbind(success, fail) ~ ., family = binomial(logit))
MOD.grouped <- glm(data = matchData4, cbind(success, fail) ~ locus.B + locus.C + locus.DRB1 + locus.A + method.HIBAG, family = binomial(logit))

mod <- glm(data = matchData3, cbind(success, fail) ~ method + locus + 0, family = binomial(logit))
mod.intercept <- glm(data = matchData3, cbind(success, fail) ~ method + locus, family = binomial(logit))
mod.grouped <- glm(data = matchData3, cbind(success, fail) ~ locus + (method=="HIBAG"), family = binomial(logit))

MODEURO <- glm(data = matchData4Euros, cbind(success, fail) ~ . + 0, family = binomial(logit))
MODEURO.intercept <- glm(data = matchData4Euros, cbind(success, fail) ~ ., family = binomial(logit))
modEuro.grouped <- glm(data = matchData3Euros, cbind(success, fail) ~ locus + (method=="HIBAG"), family = binomial(logit))

plotModel2 <- function(model) {
  c <- coefficients(model)
  nc <- gsub("`", "", names(c))
  coeffDf <- data.frame(name = nc, value = sub(pattern="(locus|method)(\\.)*", "", nc), type = sub(pattern="(locus|method)(.*)", "\\1", nc), OR = exp(c), exp(confint(model)))
  colnames(coeffDf)[5:6] <- c("OR_2.5", "OR_97.5")
  intercept.b <- any(intercept.bv <- coeffDf$name == "(Intercept)")
  if (intercept.b) {
    interceptData <- c(data.matrix(coeffDf)[intercept.bv,c("OR", "OR_2.5", "OR_97.5")])
    coeffDf <- coeffDf[!intercept.bv,]
    g.intercept.line <- geom_hline(y=interceptData/interceptData[1], color = "red", size = c(2, 1, 1), alpha = c(0.8, 0.5, 0.5)) 
  } else {
    g.intercept.line <- geom_blank()
  }
  values <- levels(coeffDf$value)
  coeffDf$value <- factor(coeffDf$value, levels = values[order(tolower(values))])
  coeffDf <- na.omit(coeffDf)
  ggplot(coeffDf, aes(x=value, y=OR, ymin = OR_2.5, ymax = OR_97.5)) + theme_bw() +
    g.intercept.line +
    geom_errorbar() + geom_point() +
    facet_grid(~ type, scale = "free_x") + labs(x="")
}

plotModel <- function(model) {
  c <- coefficients(model)
  nc <- gsub("`", "", names(c))
  coeffDf <- data.frame(name = nc, value = sub(pattern="(locus|method)(\\.)*", "", nc), type = sub(pattern="(locus|method)(.*)", "\\1", nc), OR = exp(c), exp(confint(model)))
  colnames(coeffDf)[5:6] <- c("OR_2.5", "OR_97.5")
  intercept.b <- any(intercept.bv <- coeffDf$name == "(Intercept)")
  if (intercept.b) {
    interceptData <- c(data.matrix(coeffDf)[intercept.bv,c("OR", "OR_2.5", "OR_97.5")])
    coeffDf <- coeffDf[!intercept.bv,]
    if (sum(missingCoeffs.bv <- is.na(coeffDf$OR) > 0)) {
      coeffDf[missingCoeffs.bv, 4:6] <- matrix(interceptData/interceptData[1], nr = sum(missingCoeffs.bv), nc=3, byrow=T)
    } else {
      
    }
  } else {
  }
  values <- levels(coeffDf$value)
  coeffDf$value <- factor(coeffDf$value, levels = values[order(tolower(values))])
  g.intercept.line <- geom_hline(y=1, color = "red", size = 0.8)
  ggplot(coeffDf, aes(x=value, y=OR, ymin = OR_2.5, ymax = OR_97.5)) + theme_bw() +
    g.intercept.line +
    geom_errorbar() + geom_point() +
    facet_grid(~ type, scale = "free_x") + labs(x="")
}

library(grid)
getCoeffDf <- function(model) {
  c <- coefficients(model)
  nc <- gsub("`", "", names(c))
  coeffDf <- data.frame(name = nc, value = sub(pattern="(locus|method)(\\.)*", "", nc), type = sub(pattern="(locus|method)(.*)", "\\1", nc), OR = exp(c), exp(confint(model)))
  colnames(coeffDf)[5:6] <- c("OR_2.5", "OR_97.5")
  intercept.b <- any(intercept.bv <- coeffDf$name == "(Intercept)")
  if (intercept.b) {
    interceptData <- c(data.matrix(coeffDf)[intercept.bv,c("OR", "OR_2.5", "OR_97.5")])
    coeffDf <- coeffDf[!intercept.bv,]
    if (sum(missingCoeffs.bv <- is.na(coeffDf$OR) > 0)) {
      coeffDf[missingCoeffs.bv, 4:6] <- matrix(interceptData/interceptData[1], nr = sum(missingCoeffs.bv), nc=3, byrow=T)
    } else {
    }
  } else {
  }
  values <- levels(coeffDf$value)
  coeffDf$value <- factor(coeffDf$value, levels = values[order(tolower(values))])
  return(coeffDf)
}

plotModelColors <- function(model, yLims = NULL) {
  coeffDf <- getCoeffDf(model)
  g.intercept.line <- geom_hline(y=1, color = "red", size = 0.8)
  if (!is.null(yLims)){
    g.scale <- scale_y_continuous(limits = yLims) 
  } else {
    g.scale <- geom_blank()
  }
  ggplot(coeffDf, aes(x=value, y=OR, ymin = OR_2.5, ymax = OR_97.5)) + theme_bw(base_size = 12) +
    g.intercept.line +
    geom_errorbar() + geom_point(aes(fill = value, shape = value), size = 4) +
    scale_fill_manual(values = c(rep("grey95", 4), colors_perso), guide = "none") +
    scale_shape_manual(values = c(21,22,23,24, rep(25,4)), guide = "none")+
    facet_grid(~ type, scale = "free_x") + labs(x=NULL) +
    g.scale +
    theme(axis.text.x = element_text(angle = -12, hjust = 0.3),
          plot.margin = unit(c(0.5, 1, 0., 0.), "lines"),
          #           axis.title.x = element_blank(),
          panel.margin = unit(0.5, "lines"))
}

ggsave(filename="ModelIntercept.png", 
       plot=plotModel(MOD.intercept),
       width = 8, height = 6, dpi = 150)

ggsave(filename="Model.png", 
       plot=plotModel(MOD),
       width = 8, height = 6, dpi = 150)

ggsave(filename="ModelGrouped.png", 
       plot=plotModel(mod.grouped),
       width = 8, height = 6, dpi = 150)

ggsave(filename="ModelEuroIntercept.png", 
       plot=plotModel(MODEURO.intercept),
       width = 8, height = 6, dpi = 150)

ggsave(filename="ModelEuro.png", 
       plot=plotModel(MODEURO),
       width = 8, height = 6, dpi = 150)

ggsave(filename="ModelEuroGrouped.png", 
       plot=plotModel(modEuro.grouped),
       width = 8, height = 6, dpi = 150)


# SS ------------------------------------------------------------------

pct <- c(87.0, 87.0, 85.5, 86.0,
         58.5, 70.0, 55.5, 64.0,
         82.5, 94.0, 88.5, 92.0,
         42.5, 56.0, 49.5, 53.0)
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

sampleIds <- unique(matchData1$sample.id)
matchData1R <- matchData1[ matchData1$sample.id %in% sampleIds[sample(length(sampleIds), 103)],]
matchData2R <- ddply(matchData1R, ~ locus + method, summarize, success = sum(success, na.rm = T), fail = length(sample.id)*2 - sum(success))
within(matchData2R, sumsf <- success + fail) 

mod <- glm(data = matchData2R, cbind(success, fail) ~ locus + method, family = binomial(link=logit))
summary(mod)
confint(mod)
cbind(within(matchData3Euros, prop <- success / (success + fail)), p= plogit(predict(mod)), predict(mod), residuals(mod))

matchData4R <- booleaniseFactors(matchData2R, c("method", "locus"))
matchData4R <- matchData4R[c(1:2, 4:6, 3, 8:10, 7)]
MODR.intercept <- glm(data = matchData4R, cbind(success, fail) ~ ., family = binomial(logit))

# Final output ------------------------------------------------------------

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

MODEURO.grouped <- glm(data = matchData4Euros, cbind(success, fail) ~ locus.B+ locus.DRB1 + (locus.C|locus.A)  + method.HIBAG + (`method.HLA*IMP:02`|`method.e-HLA`|method.MAGPrediction), family = binomial(logit))
names(MODEURO.grouped$coefficients)[c(4,6)] <- c("locus.AorC", "method.NON-HIBAG")
# names(MODEURO.grouped$coefficients)[4] <- "locus.AorC"
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

if (extraStuff.b <- FALSE) {
  
  # Conditional models ------------------------------------------------------
  
  stratifiedDataFileName <- "StratifiedData.Rdata"
  if (!any(grepl(stratifiedDataFileName, dir()))) {
    matchData5 <- NULL
    for (i in 1:nrow(matchData)) {
      matchData5 <- rbind(matchData5, 
                          rbind(within(matchData[i,], success <- floor(success/2)),
                                within(matchData[i,], success <- ceiling(success/2))))
    }
    save(file = stratifiedDataFileName, matchData5)
  } else {
    load(stratifiedDataFileName)
  }
  
  # write.table(file = "~/Dropbox/IMMPUTE/Modeling/DataForPA.txt", matchData5)
  # write.table(file = "~/Dropbox/IMMPUTE/Modeling/DataForPA.txt", matchData5, quote = F)
  
  stratifiedModels1FileName <- "StratifiedModels1.Rdata"
  require(survival)
  library(ordinal)
  if (!any(grepl(stratifiedModels1FileName, dir()))) {
    cmod <- clogit(data = matchData5, success ~ method + locus + strata(sample.id))
    ## Ordinal fit, different from the binomial one.
    fm1 <- clmm(factor(success) ~ method + locus + (1|sample.id), data = matchData1, Hess = TRUE)
    # fm2 <- clmm(factor(success) ~ method + locus + (1|sample.id), data = matchData, Hess = FALSE) # Just less info
    fm3 <- clmm(factor(success) ~ method + locus + (1|sample.id), data = matchData1, Hess = TRUE, nAGQ = 10)
    fm4 <- clm(factor(success) ~ method + locus, data = matchData1, Hess = TRUE) 
    # to be compared with
    fm4Binom <- glm(data = matchData1, cbind(success, 2-success) ~ method + locus, family = binomial(link=logit))
    save(file = stratifiedModels1FileName, cmod, fm1, fm2, fm3, fm4)
  } else {
    load(stratifiedModels1FileName)
  }
  
  
  # Full study of the stratified model --------------------------------------
  # See http://www.ats.ucla.edu/stat/r/dae/melogit.htm
  # and https://stat.ethz.ch/pipermail/r-sig-mixed-models/2013q1/020016.html for MCMCs
  
  # library(GGally) # helper for ggplot, not used here
  library(lme4)
  # estimate the model and store results in m
  m <- glmer(success ~ method + locus + (1|sample.id), data = matchData5, 
             family = binomial, control = glmerControl(optimizer = "bobyqa"))
  
  library(MCMCglmm)
  
}



