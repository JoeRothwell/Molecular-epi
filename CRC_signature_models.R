# Find metabolic and fatty acids signatures of WCRF score
# Requires crc1, crc2, controls and PLS models to be prepared from CRC_data_prep.R
source("CRC_get_signatures.R")

# Functions to get CC subjects with signature-derived WCRF scores
predict.scores <- function(crc, dat, mod, scorecomp.only = F){

  library(tidyverse)
  library(zoo)
  
  print(paste("Subjects in case-control: ", nrow(crc)))
  # Variables were converted to factors in CRC_data_prep.R
  
  # Put CRC CC compounds in same order as in controls dataset (fatty acids, biocrates)
  if(nrow(crc) == 877) { 
    cols <- common.cols
    var.list <- c("L_School", "Smoke_Stat")
    } else { 
    cols <- colnames(dat)
    var.list <- c("Country", "Center", "Sex")
    }
  
  cmpds <- crc %>% select(one_of(cols))
  
  # Convert variables to factors
  crc <- crc %>% mutate_at(vars(var.list), as.factor)
  
  # check colnames are the same for both sets
  if(identical(colnames(dat), colnames(cmpds)) == T) print("Identical colnames OK")
  
  # Transform data to matrix, impute, log, scale
  mat <- cmpds
  mat[mat == 0] <- NA
  mat <- na.aggregate(mat, FUN = function(x) min(x)/2)
  obs2predict <- log2(mat) %>% scale
  
  # now use predict to predict the scores of new observations (ie case-control study)
  crcscores <- data.frame(predict(mod, obs2predict))
  if(scorecomp.only == T) return(crcscores)
  
  output <- cbind(crcscores, crc) %>% group_by(Match_Caseset) %>% filter(n() == 2)

}

# Biocrates compounds overlap individual CC/control only
small <- predict.scores(crc1, ctrlA, mod1a)
large <- predict.scores(crc2, ctrlB, mod1b)
FAs <- predict.scores(CRCfa1, concs, mod2)

# Join small and large together for pooled questionnaire model
common.vars <- c("Cncr_Caco_Clrt", "Qe_Energy", "L_School", "Smoke_Stat", "Match_Caseset", 
                 "Wcrf_C_Cal", "Height_C")
small0 <- select(small, common.vars)
large0 <- select(large, common.vars)
all <- bind_rows(small0, large0)

# Biocrates compounds overlap all datasets no longer used
# Model CC status from calculated or signature-predicted score for four datasets

library(survival)
base <- Cncr_Caco_Clrt ~ Qe_Energy + L_School + Smoke_Stat + Height_C + strata(Match_Caseset)

# Biocrates small, large, fatty acids (large fasted subset did not improve OR)
fit1 <- clogit(update(base, ~. + score.1.comps), data = small)
fit2 <- clogit(update(base, ~. + Wcrf_C_Cal), data = small)
fit3 <- clogit(update(base, ~. + score.1.comps), data = large)
fit4 <- clogit(update(base, ~. + Wcrf_C_Cal), data = large)
fit5 <- clogit(update(base, ~. + score.2.comps), data = FAs)
fit6 <- clogit(update(base, ~. + Wcrf_C_Cal), data = FAs)
# Pooled questionnaire
fit0 <- clogit(update(base, ~. + Wcrf_C_Cal), data = all)

# Table for manuscript
library(broom)
t1 <- map_df(list(fit1, fit2, fit3, fit4, fit5, fit6, fit0), tidy) %>% 
  filter(str_detect(term, "score.|Wcrf")) %>% select(-(std.error : p.value)) %>%
  mutate_at(.vars = c("estimate", "conf.low", "conf.high"), exp)

score <- t1 %>% filter(term == "Wcrf_C_Cal")
sig   <- t1 %>% filter(term != "Wcrf_C_Cal")
  

# Forest plots ----
# Signature only for Biocrates small and large (all subjects in study) and fatty acids small
library(broom)
library(tidyverse)

ll2 <- list(fit3, fit1, fit5)
t2 <- map_df(ll2, tidy) %>% filter(str_detect(term, "score."))

studies <- data.frame(CC = c("B", rep("A", 2)), nvec = #map_int(ll2, 10),
            c(2330, 978, 922), metabolites = c(rep("Endogenous", 2), "Fatty acids"))

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
ma1$QEp # p-value for heterogeneity
ma1$I2 # I2

addpoly(ma1, row = 3.5, transf = exp, mlab = "", efac = 2)
text(-1.2, 3.5, "Fixed effects meta-analysis of A and B", pos = 4, cex = 0.9)
text(-1.2, 2.5, bquote(paste(I^2," = 0, ",italic(p),"-heterogeneity = 0.32")), pos = 4, cex = 0.9)

#text(-1.2, 3, bquote(paste("Fixed effects meta-analysis of A and B\n(p-heterogeneity = 0.32,",I^2," = 0)")), pos = 4, cex = 0.9)


# All models ----

library(broom)
# Old: data for forest plot: all 8 models above for Biocrates and fatty acids
ll <- list(fit7, fit8, fit3, fit4, fit5, fit6)
t1 <- map_df(ll, tidy) %>% filter(term == "score.2.comps" | term == "Wcrf_C_Cal")
studies <- data.frame(
  CC = c(rep("Small", 2), 
         #rep("Large, fast", length(ll)/3), 
         rep("Large, all", length(ll)/3), 
         rep("Small", 2)),
  nvec = map_int(ll, 10),
  metabolites = c(rep("Fatty acids", 2), rep("Biocrates", 4)),
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









