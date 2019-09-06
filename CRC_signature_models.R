# Find metabolic and fatty acids signatures of WCRF score
# Requires crc1, crc2, controls and PLS models to be prepared from CRC_data_prep.R
source("CRC_prep_data.R")

# Functions to get CC subjects with signature-derived WCRF scores
get.scores.bioc <- function(crc, dat, mod, scorecomp.only = F){

  library(tidyverse)
  colnames(dat) %>% length
  
  # ----
  
  print(paste("Subjects in case-control: ", nrow(crc)))
  # Variables were converted to factors in CRC_data_prep.R
  
  # Put CRC CC compounds in same order as in controls dataset
  cols <- colnames(dat)
  cmpds <- crc %>% select(one_of(cols))
  
  colnames(cmpds) %>% length
  # 163 variables
  
  # Convert variables to factors
  var.list <- c("Country", "Center", "Sex")
  crc <- crc %>% mutate_at(vars(var.list), as.factor)
  
  # check colnames are the same for both sets
  if(identical(colnames(dat), colnames(cmpds)) == T) print("Identical colnames OK")
 
  # Score predictions ----
  # Get signature scores and bind to original dataset
  
  # Transform data to matrix, impute, log, scale
  mat <- cmpds
  mat[mat == 0] <- NA
  mat <- na.aggregate(mat, FUN = function(x) min(x)/2)
  obs2predict <- log2(mat) %>% scale
  
  # now use predict to predict the scores of new observations (ie case-control study)
  crcscores <- data.frame(predict(mod, obs2predict))
  #score.2.comps <- predict(mod1, obs2predict)
  #if(scorecomp.only == T) return(crcscores)
  #output <- cbind(score.2.comps, crc)
  output <- cbind(crcscores, crc)

}
get.scores.FA  <- function(){
  
  # Fatty acid signatures of WCRF score
  library(tidyverse)
  library(haven)
  
  # Convert categorical co-variates to factors
  var.list <- c("L_School", "Smoke_Stat")
  
  length(colnames(concs))
  length(colnames(CRCfa))
  
  # select common FAs in the same order
  CRCfa <- CRCfa %>% select(one_of(common.cols))
  identical(colnames(CRCfa), colnames(concs))
  
  # Predict scores for CRC dataset
  CRCfa <- as.matrix(CRCfa)
  CRCfa[CRCfa == 0] <- NA
  library(zoo)
  CRCfa <- na.aggregate(CRCfa, FUN = function(x) min(x)/2)
  obs2predict <- log2(CRCfa) %>% scale %>% data.frame
  
  # now use predict to predict the scores of new observations (ie case-control study)
  crcscores <- data.frame(predict(mod2, obs2predict))
  
  #score.2.comps <- predict(mod2, obs2predict)
  
  # Put the predicted scores together with the original data, remove unpaired sample
  output <- cbind(crcscores, CRCfa1) %>% group_by(Match_Caseset) %>% filter(n() == 2)
  #output <- cbind(score.2.comps, CRCfa1) %>% group_by(Match_Caseset) %>% filter(n() == 2)
  
}

# Biocrates compounds overlap individual CC/control only
small <- get.scores.bioc(crc1, ctrlA, mod1a)
large <- get.scores.bioc(crc2, ctrlB, mod1b)

# Biocrates compounds overlap all datasets
small <- get.scores.bioc(crc1, ctrls0, mod0)
large <- get.scores.bioc(crc2, ctrls0, mod0)

FAs  <- get.scores.FA()

#large.F <- df.bioc(study = "large", fasting = T)

# Model CC status from calculated or signature-predicted score for four datasets

library(survival)
base <- Cncr_Caco_Clrt ~ Qe_Energy + L_School + Smoke_Stat + strata(Match_Caseset)

# Biocrates small
fit5 <- clogit(update(base, ~. + score.2.comps), data = small)
fit6 <- clogit(update(base, ~. + Wcrf_C_Cal), data = small)

# Biocrates large
fit3 <- clogit(update(base, ~. + score.2.comps), data = large)
fit4 <- clogit(update(base, ~. + Wcrf_C_Cal), data = large)

# Fatty acids
fit7 <- clogit(update(base, ~. + score.2.comps), data = FAs)
fit8 <- clogit(update(base, ~. + Wcrf_C_Cal), data = FAs)

# Biocrates large fasted (not used, excluding non-fasted doesn't improve OR)
#fit1 <- clogit(update(base, ~. + score.2.comps), data = large.F)
#fit2 <- clogit(update(base, ~. + Wcrf_C_Cal), data = large.F)



# Forest plots ----

# Signature only for Biocrates small and large (all subjects in study) and fatty acids small
library(broom)
library(tidyverse)

ll2 <- list(fit3, fit5, fit7)
t2 <- map_df(ll2, tidy) %>% filter(str_detect(term, "score."))

studies <- data.frame(CC = c("B", rep("A", 2)), nvec = #map_int(ll2, 10),
                      c(2330, 978, 922),
  metabolites = c(rep("Endogenous", 2), "Fatty acids"))

par(mar=c(5,4,1,2))
library(metafor)
forest(t2$estimate, ci.lb = t2$conf.low, ci.ub = t2$conf.high, refline = 1, 
       rows = c(5,6,1),
       xlab = "Odds ratio (per category increase in score)", pch = 18, 
       transf = exp, psize = 1.5, slab = studies$CC, ilab = studies[, 2:3], 
       ylim = c(0, 9), xlim = c(-1.2, 1.8),
       ilab.pos = 4, ilab.xpos = c(-0.7, -0.3))
par("usr")

text(c(-1.2, -0.7, -0.3), 8, c("Case-control", "n", "Metabolite signature"), pos = 4)
text(1.8, 8, "OR [95% CI]", pos = 2)

# Perform meta-analysis of Biocrates and add to line 3
ma1 <- rma(estimate, sei = std.error, data=t2, method="FE", subset = 1:2)

addpoly(ma1, row = 3.5, transf = exp, mlab = "", efac = 2)
text(-1.2, 3.5, "Fixed effects meta-analysis of A and B", pos = 4, cex = 0.9)
text(-1.2, 2.5, bquote(paste(I^2," = 0, ",italic(p),"-heterogeneity = 0.32")), pos = 4, cex = 0.9)

#text(-1.2, 3, bquote(paste("Fixed effects meta-analysis of A and B\n(p-heterogeneity = 0.32,",I^2," = 0)")), pos = 4, cex = 0.9)

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









