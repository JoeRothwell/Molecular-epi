# PLS to find signature of GRS
source("GRS_biocrates.R")

# Get signature of PRS as a PLS object
signatureGRS <- function() {
  # First adjust matrix for study, centre, sex, batch, BMI
  logmat <- select(controls, -(Match_Caseset:GRSgroup)) %>% as.matrix
  adj    <- function(x) residuals(lm(x ~ study + Country + Sex, data = controls))
  adjmat <- apply(logmat, 2, adj)
  
  # PLS model for metabolic signature of GRS
  # Subset GRS and bind to residuals-adjusted matrix
  grs <- data_frame(score = controls$GRS)
  plsdata <- cbind(grs, adjmat)
  
  library(pls)
  mod <- plsr(score ~ ., data = plsdata, validation = "CV")
  #summary(mod)
  
  # Find the number of dimensions with lowest cross validation error
  cv <- RMSEP(mod)
  plot(RMSEP(mod), legendpos = "topright")
  
  # Calculate optimal number of dimensions and rerun model
  best.dims <- which.min(cv$val[estimate = "adjCV", , ]) - 1
  mod <- plsr(score ~ ., data = plsdata, ncomp = best.dims)
  
  return(mod)
  
  # explained variances
  explvar(mod)
  
  # Coefficients and variable importance
  coefficients <- coef(mod)
  sum.coef <- sum(sapply(coefficients, abs))
  coefficients <- coefficients * 100 / sum.coef
  coefficients <- sort(coefficients[, 1 , 1])
  # plot(coefficients)
  
  # Get top and bottom deciles of compound coefficients
  df <- data.frame(as.list(coefficients)) %>% gather(Cmpd, VIP)
  
  qs <- quantile(coefficients, probs = seq(0, 1, 0.05))
  df1 <- df[df$VIP > qs[18], ]
  df2 <- df[df$VIP < qs[4], ]
  
  # if(modonly == T) return(list(df1, df2))
  # Vector of colours for plot points
  vec <- ifelse(df$VIP > qs[18] | df$VIP < qs[4], "black", "grey")
  
  # Now plot data, adding text
  plot(coefficients, pch = 17, col=vec, xlab = "", ylab = "Variable Importance",
       main = paste(nrow(plsdata), "subjects, optimal dimensions =", best.dims))
  text(nrow(df) - nrow(df1):1, df1$VIP, df1$Cmpd, pos=2, cex = 0.6)
  text(1:nrow(df2), df2$VIP, df2$Cmpd, pos=4, cex=0.6)
  abline(a=0, b=0, lty = "dotted")
}
mod <- signatureGRS()


# Need to get Idepics of 740 cases/controls with GWAS data and model them with signature instead of PRS

