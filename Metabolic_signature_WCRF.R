# Find metabolic and fatty acids signatures of WCRF score

# Get biocrates and fatty acids controls dataset for PLS model

wcrf.crc <- function(study = c("large", "small"), fasting = T){


  library(tidyverse)
  
  # EPIC controls dataset for modelling (controls)
  ctrl <- readRDS("Biocrates data controls.rds")
  outliers <- c(3600, 3722)
  ob <- ctrl[-outliers, ]
  
  # make new numeric variable for batch and exclude colorectal controls and non-fasted samples
  ob <- ob %>% mutate(batch_no = as.numeric(flatten(str_extract_all(Batch_MetBio, "[0-9]+")))) %>% 
    filter(Fasting_C == "Yes", Study != "Colonrectum")
  
  # Subset biocrates compounds
  controls <- ob %>% select(Acylcarn_C0 : Sugars_H1)
  colnames(controls) %>% length
  # 147 variables
  
  
  # ----
  
  
  # Prepare CRC dataset (from Jelena, ~ 1200 case-control pairs), call it setA
  # Join scores and filter out non-fasted samples. Remove odd cases or controls
  
  library(haven)
  wcrf <- read_dta("Wcrf_Score.dta") %>% select(Idepic, Wcrf_C_Cal)
  
  crcA <- read_csv("biocrates_p150.csv")
  crcA <- crcA %>% left_join(wcrf, by = "Idepic")
    
  crcA <- if(fasting == T) crcA %>% filter(Fasting_C == 2) else crcA
  
  print(paste("number of observations = ", nrow(crcA)))
  
  
  setA <- crcA %>% select(Acylcarn_C0 : Sugars_H1, -starts_with("Outdq_"))
  colnames(setA) %>% length
  # 163 variables
  
  # Get common compound between controls and CRC CC
  common_cols <- intersect(colnames(controls), colnames(setA))
  #126 compounds in common
  
  
  # ----
  
  
  # Prepare CRC dataset (from Bertrand, ~ 490 case-control pairs), call it set2b
  # Remove unpaired observations
  crcB <- readRDS("CRC_smallerCC.rds") %>% group_by(Match_Caseset) %>% filter(n() == 2) %>% ungroup
  
  print(paste("number of observations = ", nrow(crcB)))
  
  setB <- crcB %>% select(Acylcarn_C0 : Sugars_H1, -Batch_MetBio)
  colnames(setB) %>% length
  
  #common_cols <- intersect(colnames(controls), colnames(setB))
  #145 compounds in common
  
  
  # ----
  
  
  # Get common cols between all three datasets (controls, small CC, big CC)
  common_cols2 <- intersect(common_cols, colnames(setB))
  # Glyceroph_Lysopc_A_C24_0 is removed
  
  
  # ----  
  
  # subset and reorder both datasets to get the same 126 compounds in the same order
  controls <- controls %>% select(one_of(common_cols2))
  setA <- setA %>% select(one_of(common_cols2))
  setB <- setB %>% select(one_of(common_cols2))
  #dim(controls)
  #dim(setA)
  
  
  # check colnames are the same for both sets
  identical(colnames(controls), colnames(setA))
  
  
  # Perform PLS to get metabolic signature of WCRF score on controls dataset
  
  concs <- as.matrix(controls)
  
  # replace zero with NA
  concs[concs == 0] <- NA
  
  # Impute missing concentrations with half minimum value
  # Scale data to give metabolites equal importance
  library(zoo)
  concs1 <- na.aggregate(concs, FUN = function(x) min(x)/2)
  logconcs <- log2(concs1) %>% scale
  
  # adjust matrix for study, centre, sex
  library(lme4)
  adj5 <- function(x) residuals(lmer(x ~ Center + batch_no + Sex + (1|Study), data = ob))
  adjmat <- apply(logconcs, 2, adj5)
  
  # Subset calibrated scores
  score <- data_frame(score = ob$Wcrf_C_Cal)
  
  # Data setup. Must be a df with Bind scores to log matrix
  plsdata <- cbind(score, adjmat) %>% filter(!is.na(score))
  
  
  # ---- PLS model
  
  
  library(pls)
  mod <- plsr(score ~ ., data = plsdata, validation = "CV")
  #summary(mod)
  
  # Find the number of dimensions with lowest cross validation error
  cv <- RMSEP(mod)
  plot(RMSEP(mod), legendpos = "topright")
  
  # Calculate optimal number of dimensions
  best.dims <- which.min(cv$val[estimate = "adjCV", , ]) - 1
  
  # Rerun the model with the optimum number of components
  mod <- plsr(score ~ ., data = plsdata, ncomp = best.dims)
  
  # explained variances
  explvar(mod)
  
  # Plots: prediction, scores, loadings
  plot(mod)
  plot(mod, plottype = "scores")
  plot(mod, "loadings", legendpos = "topleft")
  
  # Plot the variable importance
  coefficients <- coef(mod)
  sum.coef <- sum(sapply(coefficients, abs))
  coefficients <- coefficients * 100 / sum.coef
  coefficients <- sort(coefficients[, 1 , 1])
  plot(coefficients)
  
  # Most important compounds
  par(mai = c(1,2,0,0.5), mfrow = c(2,1))
  barplot(head(coefficients, 10), horiz = T, las=1, col="red", xlab="Variable importance")
  barplot(tail(coefficients, 10), horiz = T, las=1, col="dodgerblue", xlab="Variable importance")
  
  
  # ---- Score predictions and model of C/C status on score
  
  # Transform data to matrix, impute, log, scale
  if(study == "large") mat <- setA else mat <- setB
  
  mat[mat == 0] <- NA
  mat <- na.aggregate(mat, FUN = function(x) min(x)/2)
  obs2predict <- log2(mat) %>% scale
  
  # now use predict to predict the scores of new observations (ie case-control study)
  crcscores <- data.frame(predict(mod, ncomp = 2, obs2predict))
  
  if(study == "large") crc <- cbind(crcscores, crcA) else crc <- cbind(crcscores, crcB)

}
fa.crc <- function(){
  
  # Fatty acid signatures of WCRF score
  library(tidyverse)
  
  # Get test dataset (CRC case control from Elom)
  library(haven)
  wcrf <- read_dta("Wcrf_Score.dta") %>% select(Idepic, Wcrf_C_Cal)
  
  # Join scores to fatty acids data
  CRCfa1 <- read_dta("Database_Fatty acids.dta") %>% left_join(wcrf, by = "Idepic")
  CRCfa <- CRCfa1 %>% select(P14_0 : PCLA_9t_11c) 
  
  # Get dataset for PLS modelling (all EPIC controls). Exclude compounds with many missings
  fa.scores <- readRDS("FA_WCRF_scores.rds")
  concs <- fa.scores %>% select(P14_0 : PCLA_9t_11c, -P24_0, -P20_0)
  
  # Check missing values
  #library(Amelia)
  #missmap(concs)
  #missmap(CRCfa)
  
  length(colnames(concs))
  length(colnames(CRCfa))
  common_cols <- intersect(colnames(concs), colnames(CRCfa))
  
  # select common FAs in the same order
  CRCfa <- read_dta("Database_Fatty acids.dta") %>% select(one_of(common_cols)) 
  concs <- fa.scores %>% select(one_of(common_cols))
  
  identical(colnames(CRCfa), colnames(concs))
  
  # PLS to get signature of fatty acid metabolism for high and low scorers
  # Set zeros to NA and impute with half miminum conc, log transform
  # Bind to data frame
  
  concs <- as.matrix(concs)
  concs[concs == 0] <- NA
  library(zoo)
  concs <- na.aggregate(concs, FUN = function(x) min(x)/2)
  logconcs <- log2(concs) %>% scale
  
  # Bind WCRF scores to log2 concs
  plsdata <- data.frame(score = fa.scores$Wcrf_C_Cal, logconcs) %>% filter(!is.na(score))
  
  library(pls)
  mod <- plsr(score ~ ., data = plsdata, validation = "CV")
  # summary(mod)
  
  # Find the number of dimensions with lowest cross validation error
  cv <- RMSEP(mod)
  plot(RMSEP(mod), legendpos = "topright")
  
  # Calculate optimal number of dimensions
  best.dims <- which.min(cv$val[estimate = "adjCV", , ]) - 1
  
  # Rerun the model with the optimum number of components
  mod <- plsr(score ~ ., data = plsdata, ncomp = best.dims)
  
  # explained variances
  explvar(mod)
  
  # Plots: prediction, scores, loadings
  plot(mod)
  plot(mod, plottype = "scores")
  plot(mod, "loadings", legendpos = "topleft")
  
  # Plot the variable importance
  coefficients <- coef(mod)
  sum.coef <- sum(sapply(coefficients, abs))
  coefficients <- coefficients * 100 / sum.coef
  coefficients <- sort(coefficients[, 1 , 1])
  plot(coefficients)
  
  # Most important compounds
  par(mai = c(1,2,0,0.5), mfrow = c(2,1))
  par(mfrow = c(2,1))
  barplot(head(coefficients, 10), horiz = T, las=1, col="red", xlab="Variable importance")
  barplot(tail(coefficients, 10), horiz = T, las=1, col="dodgerblue", xlab="Variable importance")
  #dev.off()
  
  # Predict scores for CRC dataset ---------------------------------------------------------------------
  
  CRCfa <- as.matrix(CRCfa)
  CRCfa[CRCfa == 0] <- NA
  library(zoo)
  CRCfa <- na.aggregate(CRCfa, FUN = function(x) min(x)/2)
  obs2predict <- log2(CRCfa) %>% scale %>% data.frame
  
  
  #table(CRCfa$Match_Caseset, CRCfa$Cncr_Caco_Clrt)
  #write.csv(df, file = "CRC biocrates data.csv")
  
  # now use predict to predict the scores of new observations (ie case-control study)
  crcscores <- data.frame(predict(mod, ncomp = 2, obs2predict))
  
  # Put the predicted scores together with the original data, remove unpaired sample
  crc <- cbind(crcscores, CRCfa1) %>% group_by(Match_Caseset) %>% filter(n() == 2)
  
}

# Run functions to get datasets
large <- wcrf.crc()
small <- wcrf.crc(study = "small")
large.nofast <- wcrf.crc(fasting = F)
FAs <- fa.crc()

# Models: Biocrates metabolites ----

# Large study
library(survival)
fit1 <- clogit(Cncr_Caco_Clrt ~ score.2.comps + strata(Match_Caseset), data = large)
fit2 <- clogit(Cncr_Caco_Clrt ~ score.2.comps + strata(Match_Caseset), data = large.nofast)

# Calculated scores
fit3 <- clogit(Cncr_Caco_Clrt ~ Wcrf_C_Cal + strata(Match_Caseset), data = large)
fit4 <- clogit(Cncr_Caco_Clrt ~ Wcrf_C_Cal + strata(Match_Caseset), data = large.nofast)

# Small study
fit5 <- clogit(Cncr_Caco_Clrt ~ score.2.comps + strata(Match_Caseset), data = small)
fit6 <- clogit(Cncr_Caco_Clrt ~ Wcrf_C_Cal + strata(Match_Caseset), data = small)

# Fatty acids: metabolic and calculated scores
metabolic <- clogit(Cncr_Caco_Clrt ~ score.2.comps + strata(Match_Caseset), data = FAs)
calculated <- clogit(Cncr_Caco_Clrt ~ Wcrf_C_Cal + strata(Match_Caseset), data = FAs)

# Forest plots ----

# Biocrates
library(broom)
t1 <- bind_rows(tidy(fit1), tidy(fit2), tidy(fit3), tidy(fit4), tidy(fit5), tidy(fit6))
nvec <- rev(c(980, 980, 606, 2254, 614, 2370))
lvec <- c("WCRF score", "Metabolites", "Metabolites (F)", "Metabolites", "WCRF score (F)", "WCRF score")
study <- c(rep("Small C/C", 2), rep("Large C/C", 4))

library(metafor)
par(mar=c(5,4,2,2))
forest(x = exp(t1$estimate), ci.lb = exp(t1$conf.low), ci.ub = exp(t1$conf.high),
       refline = 1, xlab = "Odds ratio (per unit increase in WCRF score)", pch = 18, 
       psize = 1.5, slab = study, ilab = lvec, ilab.xpos = -0.2, xlim = c(-1, 2))

text(-1, 8, "Study", pos = 4)
text(-0.2, 8, "Parameter (F = fasting)", pos = 3)
text(2, 8, "OR [95% CI]", pos = 2)

# Fatty acids
t2 <- bind_rows(tidy(metabolic), tidy(calculated))

par(mar=c(5,4,2,2))
forest(x = exp(t2$estimate), ci.lb = exp(t2$conf.low), ci.ub = exp(t2$conf.high),
       refline = 1, xlab = "Odds ratio (per unit increase in WCRF score)", pch = 18, 
       psize = 1.5, 
       slab = c("High WCRF score", "Fatty acid signature"),
       ilab.xpos = -0.2, xlim = c(-0.5, 1.5))

text(nvec, x = 0.1, y=c(1:4))













