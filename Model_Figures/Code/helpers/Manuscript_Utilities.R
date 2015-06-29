# Script defining some simple utility function for teh creation of figures for the immpute manuscript.
#
# Copyright Antoine Lizee 09/2014 antoine.lizee@gmail.com


# Plotting helpers --------------------------------------------------------

printGGplot <- function(plot, file, formats = c("png", "tiff", "pdf"), res, ...) {
  # This helper creates multiple outputs from the same graphm for visualisation / publication purposes.
  for (format in formats) {
    ggsave(filename = paste(file, format, sep = "." ), plot = plot, dpi = res, ...)
  }
}

cmToInches <- function(x) {
  x / 2.54  
}


# Data Preprocessing for the models ---------------------------------------

booleaniseFactors <- function(df, factorNames) {
  # This creates booleans from categorical variables
  for (cn in factorNames) {
    for (v in unique(c <- df[[cn]])) {
      df[paste0(cn, ".", v)] <- as.numeric((c == v))
    }
    df[cn] <- NULL
  }
  return(df)
}

swapColumns <- function(df, colName1, colName2, fixed = T) {
  # This rather long and verbose function has for only purpose switching columns in a data.frame for the sake of modeling tweaks.
  # The colNames can be patterns
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


# Significance testing helpers --------------------------------------------

printWaldTest <- function(w, digits = 2) {
  v <- w[["result"]][["chi2"]]
  cat("X2 = ", format(v["chi2"], digits = digits, nsmall = 1), 
      ", df = ", v["df"], ", P(> X2) = ", format(v["P"], digits = digits, 
                                                 nsmall = 1), "\n", sep = "")
}

cVec <- function(i,j,n) {
  v <- t(rep(0, n))
  v[i] <- 1
  v[j] <- -1
  return(v)
}

printWaldSignif <- function(mod, i, j) {
  coeffs <- names(coef(mod))
  hyp <- paste(coeffs[i], coeffs[j], sep = " == ")
  waldCoeffs <- cVec(i, j, length(coeffs))
  cat("## Testing hypothesis with Wald-Test:\n #", hyp, "#\n ")
  w <- wald.test(b = coef(mod), Sigma = vcov(mod), L = waldCoeffs)
  printWaldTest(w, 5)
  cat("\n")
  return(structure(w[["result"]][["chi2"]]['P'], names = hyp) )
}


# Plotting models ---------------------------------------------------------

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

plotModelColors <- function(model, yLims = NULL, coeffDf = getCoeffDf(model)) {
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


