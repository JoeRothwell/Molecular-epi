# Compute Biocrates and fatty acid signatures of WCRF score
source("CRC_data_prep.R")

get.Biocrates.sig <- function(fasting = T){
  
  library(tidyverse)
  library(lme4)
  library(zoo)
  library(pls)
  
  # Prepare controls matrix. Replace zero, impute with half mins, scale
  concs <- as.matrix(controls)
  concs[concs == 0] <- NA
  concs1 <- na.aggregate(concs, FUN = function(x) min(x)/2)
  logconcs <- log2(concs1) %>% scale
  
  # adjust matrix for study, centre, sex, batch, BMI
  adj5 <- function(x) residuals(lmer(x ~ Center + batch_no + Sex + Bmi_C + (1|Study), data = ctrl))
  adjmat <- apply(logconcs, 2, adj5)
  
  # Subset calibrated scores
  score <- data_frame(score = ctrl$Wcrf_C_Cal)
  
  # Data setup. Must be a df with Bind scores to log matrix
  plsdata <- cbind(score, adjmat) %>% filter(!is.na(score))
  
  # PLS model
  # For metabolic signature of WCRF score on controls dataset
  set.seed(111)
  mod <- plsr(score ~ ., data = plsdata, validation = "CV")
  
  # Find the number of dimensions with lowest cross validation error
  cv <- RMSEP(mod)
  plot(RMSEP(mod), legendpos = "topright")
  
  # Calculate optimal number of dimensions and rerun model
  best.dims <- which.min(cv$val[estimate = "adjCV", , ]) - 1
  mod <- plsr(score ~ ., data = plsdata, ncomp = best.dims)
  
  # explained variances
  # explvar(mod)
  
  # Prediction, scores, loadings plots
  # plot(mod)
  # plot(mod, plottype = "scores")
  # plot(mod, "loadings", legendpos = "topleft")
  
}
mod1 <- get.Biocrates.sig()

get.FA.sig  <- function(cor.data = F){
  
  # Fatty acid signatures of WCRF score
  library(tidyverse)
  library(lme4)
  library(zoo)
  library(pls)
  
  # PLS to get signature of fatty acid metabolism for high and low scorers
  # Set zeros to NA and impute with half miminum conc, log transform
  # Bind to data frame
  
  common_cols <- output[[2]]
  fa.ctrl <- output[[1]]
  
  # select common FAs in the same order. FAs is for correlation with Biocrates compounds
  FAs   <- read_dta("Database_Fatty acids.dta") %>% select(Idepic, one_of(common_cols)) 
  CRCfa <- FAs %>% select(-Idepic)
  concs <- fa.ctrl %>% select(one_of(common_cols))
  #identical(colnames(CRCfa), colnames(concs))
  
  if(cor.data == T) return(FAs)
  
  concs <- as.matrix(concs)
  concs[concs == 0] <- NA
  concs <- na.aggregate(concs, FUN = function(x) min(x)/2)
  logconcs <- log2(concs) %>% scale
  
  # adjust matrix for study, centre, sex, batch, BMI
  #getresiduals <- function(x) residuals(lm(x ~ STUDY + LABO + Country, data = fa.scores))
  getresiduals <- function(x) residuals(lmer(x ~ (1|STUDY) + LABO + Country, data = fa.ctrl))
  resmat <- apply(logconcs, 2, getresiduals)
  
  # Bind WCRF scores to log2 concs
  plsdata <- data.frame(score = fa.ctrl$Wcrf_C_Cal, resmat) %>% filter(!is.na(score))
  set.seed(111)
  mod <- plsr(score ~ ., data = plsdata, validation = "CV")
  
  # Find the number of dimensions with lowest cross validation error
  cv <- RMSEP(mod)
  plot(RMSEP(mod), legendpos = "topright")
  
  # Calculate optimal number of dimensions
  best.dims <- which.min(cv$val[estimate = "adjCV", , ]) - 1
  
  # Rerun the model with the optimum number of components
  mod <- plsr(score ~ ., data = plsdata, ncomp = best.dims)
  
  # explained variances
  # explvar(mod)
  
  # Plots: prediction, scores, loadings
  #plot(mod)
  #plot(mod, plottype = "scores")
  #plot(mod, "loadings", legendpos = "topleft")
  
}
mod2 <- get.FA.sig()

# Produce tables of important compounds, using compound metadata to get proper names

plot.Biocrates.sig <- function(mod, no.cmpds = 7){
  
  library(tidyverse)
  
  cmpd.meta <- read.csv("Biocrates_cmpd_metadata.csv")
  
  # Coefficients and variable importance. First subset one-matrix array to get matrix
  k <- 126
  coeff <- data.frame(value = round(coef(mod)[1:k, 1, 1], 3))
  dat <- coeff %>% 
    mutate(sm  = sum(abs(value)), 
           VIP = round((value*100)/sm, 2), 
           Compound = rownames(coeff))
  dat <- dat %>% left_join(cmpd.meta, by = "Compound")
  
  # Subset appropriate columns and print influential compounds
  dat <- dat %>% select(compound = displayname, Coefficient = value, VIP) %>% arrange(VIP)
  
  # choose number of influential compounds to plot
  n_top <- no.cmpds
  
  df1 <- top_n(dat, n_top)  %>% arrange(-VIP)
  df2 <- top_n(dat, -n_top) %>% arrange(VIP)
  
  print(df1)
  print(df2)
  
  # Vector of black and grey for plot points
  vec <- c( rep("black", n_top), rep("grey", k-(n_top*2)), rep("black", n_top) )
  
  # Now plot data, adding text
  plot(sort(coeff$value), pch = 17, col=vec, xlab = "", ylab = "Coefficient",
       main = paste(nrow(mod$scores), "fasted subjects, optimal dimensions =", mod$ncomp))
  # High and low labels
  text(nrow(dat) : (nrow(dat)-n_top), df1$Coefficient, df1$compound, pos=2, cex = 0.6)
  text(1:nrow(df2), df2$Coefficient, df2$compound, pos=4, cex=0.6)
  abline(a=0, b=0, lty = "dotted")
  
  output <- bind_rows(df1, df2)
  
}
table3a <- plot.Biocrates.sig(mod1)

plot.FA.sig  <- function(mod){
  
  library(tidyverse)
  
  cmpd.meta <- read.csv("FA_compound_data.csv")
  
  # Plot the variable importance
  k <- 34
  coeff <- data.frame(value = round(coef(mod)[1:k, 1, 1], 3))
  dat <- coeff %>% 
    mutate(sm  = sum(abs(value)), 
           VIP = round((value*100)/sm, 2), 
           Compound = rownames(coeff))
  
  dat <- dat %>% left_join(cmpd.meta, by = "Compound")
  
  # Subset appropriate columns and print influential compounds
  dat <- dat %>% select(compound = displayname2, Coefficient = value, VIP) %>% arrange(VIP)
  
  # choose number of influential compounds to plot
  n_top <- 5
  
  df1 <- top_n(dat, n_top)  %>% arrange(-VIP)
  df2 <- top_n(dat, -n_top) %>% arrange(VIP)
  
  print(df1)
  print(df2)
  
  # Vector of colours for plot points
  vec <- c( rep("black", n_top), rep("grey", k-(n_top*2)), rep("black", n_top) )
  
  # Now plot data, adding text
  plot(sort(coeff$value), pch = 17, col=vec, xlab = "", ylab = "Coefficient",
       main = paste("Fatty acids:", nrow(mod$scores), "fasted subjects, optimal dimensions =", mod$ncomp))
  # High and low compounds
  text(nrow(dat) : (nrow(dat)-n_top), df1$Coefficient, df1$compound, pos=2, cex = 0.6)
  text(1:nrow(df2), df2$Coefficient, df2$compound, pos=4, cex=0.6)
  abline(a=0, b=0, lty = "dotted")
  
  output <- bind_rows(df1, df2)
}
table3b <- plot.FA.sig(mod2)

# Save workspace (for .Rmd file)
# save.image(file="metabolic_signatures.Rdata")
