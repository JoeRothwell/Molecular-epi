# Find metabolic and fatty acids signatures of WCRF score
# Requires crc1, crc2, controls and PLS models to be prepared from CRC_data_prep.R
source("CRC_data_prep.R")

# Functions to get CC subjects with signature-derived WCRF scores
df.bioc <- function(study = c("large", "small"), fasting = T, scorecomp.only = F){

  library(tidyverse)
  # EPIC controls dataset for signature. Data prep is done in CRC_data_prep.R
  # First dataset, 3771 obs; updated November 2018 7191 obs
  ob <- ctrl
  print(paste(nrow(ob), "Controls read"))
  
  # Subset biocrates compounds
  #controls <- ob %>% select(Acylcarn_C0 : Sugars_H1)
  #controls <- ob %>% select(Glyceroph_Lysopc_A_C16_0 : Acylcarn_C5_M_Dc)
  
  controls <- ob %>% select(matches("Acylcarn_|Aminoacid_|Biogenic_|Glyceroph_|Sphingo_|Sugars_"), 
                          -starts_with("Outdq"))
  
  zerocols <- apply(controls, 2, function(x) sum(x, na.rm = T)) != 0
  controls <- controls[, zerocols]
  
  colnames(controls) %>% length # 147 variables
  
  # Large CRC metabolomics subset from Jelena, ~ 1200 case-control pairs), call it crcA
  # Join scores and filter out non-fasted samples
  crcA <- crc2
  crcA <- if(fasting == T) crcA %>% filter(Fasting_C == 2) else crcA
  print(paste("Subjects in larger case-control: ", nrow(crcA)))
  # Variables were converted to factors in CRC_data_prep.R
  
  setA <- crcA %>% select(Acylcarn_C0 : Sugars_H1, -starts_with("Outdq_"))
  colnames(setA) %>% length
  # 163 variables
  
  # Get common compound between controls and CRC CC
  common_cols <- intersect(colnames(controls), colnames(setA))   # 126 compounds in common
  
  # Small CRC metabolomics subset----
  # (from Bertrand, ~ 490 case-control pairs), call it crcB
  crcB <- crc1
  print(paste("Subjects in smaller case-control:", nrow(crcB)))
  
  setB <- crcB %>% select(Acylcarn_C0 : Sugars_H1, -Batch_MetBio)
  colnames(setB) %>% length
  
  # Convert variables to factors
  var.list <- c("Country", "Center", "Sex")
  crcB <- crcB %>% mutate_at(vars(var.list), as.factor)
  #common_cols <- intersect(colnames(controls), colnames(setB))   # 145 compounds in common

  # Get common cols between all three datasets (controls, small CC, big CC)
  common_cols2 <- intersect(common_cols, colnames(setB))   # Glyceroph_Lysopc_A_C24_0 is removed
  
  # subset and reorder both datasets to get the same 126 compounds in the same order
  controls <- controls %>% select(one_of(common_cols2))
  setA <- setA %>% select(one_of(common_cols2))
  setB <- setB %>% select(one_of(common_cols2))
  
  # check colnames are the same for both sets
  identical(colnames(controls), colnames(setA))
 
  # Score predictions ----
  # Get signature scores and bind to original dataset
  
  # Transform data to matrix, impute, log, scale
  if(study == "large") mat <- setA else mat <- setB
  
  mat[mat == 0] <- NA
  mat <- na.aggregate(mat, FUN = function(x) min(x)/2)
  obs2predict <- log2(mat) %>% scale
  
  # now use predict to predict the scores of new observations (ie case-control study)
  crcscores <- data.frame(predict(mod1, ncomp = 2, 
                                  obs2predict))
  #score.2.comps <- predict(mod1, obs2predict)
  #if(scorecomp.only == T) return(crcscores)
  #if(study == "large") output <- cbind(score.2.comps, crcA) else output <- cbind(score.2.comps, crcB)
  if(study == "large") output <- cbind(crcscores, crcA) else output <- cbind(crcscores, crcB)

}
df.FAs  <- function(){
  
  # Fatty acid signatures of WCRF score
  library(tidyverse)
  library(haven)
  
  # Get WCRF scores
  wcrf <- read_dta("Wcrf_Score.dta") %>% select(Idepic, Wcrf_C_Cal)
  
  # Get CRC dataset from Elom and join WCRF scores. Convert categorical co-variates to factors
  var.list <- c("L_School", "Smoke_Stat")
  CRCfa1 <- read_dta("Database_Fatty acids.dta") %>% left_join(wcrf, by = "Idepic") %>% mutate_at(vars(var.list), as.factor)
  
  # Subset Biocrates compounds
  CRCfa <- CRCfa1 %>% select(P14_0 : PCLA_9t_11c) 
  
  # Get dataset for PLS modelling (all EPIC controls). Exclude compounds with many missings
  # Note: new version from Carine received 18/11/2018 with technical covariates
  fa.scores <- readRDS("FA_WCRF_scores1.rds") %>% filter(STUDY != "Colorectum")
  fa.scores$N_Serie <- as.numeric(fa.scores$N_Serie)
  
  # convert categorical variables to factors
  var.list <- c("Country", "Center", "STUDY", "LABO")
  fa.scores <- fa.scores %>% mutate_at(vars(var.list), as.factor)
  
  concs <- fa.scores %>% select(P14_0 : PCLA_9t_11c, -P24_0, -P20_0)

  #library(Amelia)
  #missmap(concs)
  #missmap(CRCfa)
  
  length(colnames(concs))
  length(colnames(CRCfa))
  common_cols <- intersect(colnames(concs), colnames(CRCfa))
  
  # select common FAs in the same order
  CRCfa <- read_dta("Database_Fatty acids.dta") %>% select(one_of(common_cols)) 
  identical(colnames(CRCfa), colnames(concs))
  
  # Predict scores for CRC dataset
  CRCfa <- as.matrix(CRCfa)
  CRCfa[CRCfa == 0] <- NA
  library(zoo)
  CRCfa <- na.aggregate(CRCfa, FUN = function(x) min(x)/2)
  obs2predict <- log2(CRCfa) %>% scale %>% data.frame
  
  # now use predict to predict the scores of new observations (ie case-control study)
  crcscores <- data.frame(predict(mod2, ncomp = 2, 
                                  obs2predict))
  
  #score.2.comps <- predict(mod2, obs2predict)
  
  # Put the predicted scores together with the original data, remove unpaired sample
  output <- cbind(crcscores, CRCfa1) %>% group_by(Match_Caseset) %>% filter(n() == 2)
  #output <- cbind(score.2.comps, CRCfa1) %>% group_by(Match_Caseset) %>% filter(n() == 2)
  
}

small   <- df.bioc(study = "small")
large   <- df.bioc(study = "large", fasting = F)
large.F <- df.bioc(study = "large")
FAs     <- df.FAs()

# Model CC status from calculated or signature-predicted score for four datasets

library(survival)
base <- Cncr_Caco_Clrt ~ Qe_Energy + L_School + Smoke_Stat + strata(Match_Caseset)

# Biocrates large fasted
fit1 <- clogit(update(base, ~. + score.2.comps), data = large.F)
fit2 <- clogit(update(base, ~. + Wcrf_C_Cal), data = large.F)

# Biocrates large all
fit3 <- clogit(update(base, ~. + score.2.comps), data = large)
fit4 <- clogit(update(base, ~. + Wcrf_C_Cal), data = large)

# Biocrates small
fit5 <- clogit(update(base, ~. + score.2.comps), data = small)
fit6 <- clogit(update(base, ~. + Wcrf_C_Cal), data = small)

# Fatty acids
fit7 <- clogit(update(base, ~. + score.2.comps), data = FAs)
fit8 <- clogit(update(base, ~. + Wcrf_C_Cal), data = FAs)



# Forest plots ----

# Signature only for Biocrates small and large (all subjects in study) and fatty acids small
ll2 <- list(fit3, fit5, fit7)
library(broom)
library(tidyverse)
t2 <- map_df(ll2, tidy) %>% filter(term == "score.2.comps")
studies <- data.frame(CC = c("Large", rep("Small", 2)), nvec = map_int(ll2, 10),
  metabolites = c(rep("Biocrates", 2), "Fatty acids"))

par(mar=c(5,4,1,2))
library(metafor)
forest(t2$estimate, ci.lb = t2$conf.low, ci.ub = t2$conf.high, refline = 1, 
       xlab = "Odds ratio CRC (per unit increase in score)", pch = 18, 
       transf = exp, psize = 1.5, slab = studies$CC, ilab = studies[, 2:3], 
       ylim = c(0, 6),
       ilab.pos = 4, ilab.xpos = c(-0.8, -0.4), xlim = c(-1.2, 2))

text(c(-1.2, -0.8, -0.4), 5, c("Study", "n", "Metabolites"), pos = 4)
text(2, 5, "OR [95% CI]", pos = 2)


# Meta-analysis
par(mar=c(5,4,0,2))
ma1 <- rma(estimate, sei = std.error, data=t2, method="FE", subset = 1:2)
forest(ma1, transf = exp, refline = 1, slab = c("Large", "Small"), xlab = "OR", efac = 4)
hh <- par("usr")

text(hh[1], 4, "Study", pos = 4)
text(hh[2], 4, "OR [95% CI]", pos = 2)

ma2 <- rma(estimate, sei = std.error, data=t2, method="REML", subset = 1:2)
forest(ma2, transf = exp, refline = 1, slab = c("Large", "Small"), xlab = "OR", efac = 4)

text(hh[1], 4, "Study", pos = 4)
text(hh[2], 4, "OR [95% CI]", pos = 2)



# All models ----

library(broom)
# Old: data for forest plot: all 8 models above for Biocrates and fatty acids
ll <- list(fit7, fit8, fit1, fit2, fit3, fit4, fit5, fit6)
t1 <- map_df(ll, tidy) %>% filter(term == "score.2.comps" | term == "Wcrf_C_Cal")
studies <- data.frame(
  CC = c(rep("Small", 2), rep("Large, fast", length(ll)/3), rep("Large, all", length(ll)/3), rep("Small", 2)),
  nvec = map_int(ll, 10),
  metabolites = c(rep("Fatty acids", 2), rep("Biocrates", 6)),
  mod = rep(c("Signature", "WCRF score"), 4)
  )

library(metafor)
par(mar=c(5,4,1,2))
forest(t1$estimate, ci.lb = t1$conf.low, ci.ub = t1$conf.high,
       refline = 1, xlab = "Odds ratio CRC (per unit increase in score)", pch = 18, 
       transf = exp, psize = 1.5, slab = studies$CC, ilab = studies[, 2:4], 
       rows = c(1:2, 4:5, 7:8, 10:11),
       ylim = c(0, 14),
       ilab.pos = 4, ilab.xpos = c(-0.8, -0.6, -0.2), 
       xlim = c(-1.2, 2))

text(c(-1.2, -0.8, -0.6, -0.2), 13, c("Study", "n", "Metabolites", "Predictor"), pos = 4)
text(2, 13, "OR [95% CI]", pos = 2)









