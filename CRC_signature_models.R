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

# Predict for colon (proximal and distal) cancer only
small.colon <- predict.scores(colon1, ctrlA, mod1a)
large.colon <- predict.scores(colon2, ctrlB, mod1b)
large.rectal <- predict.scores(rectal2, ctrlB, mod1b)

nrow(small); nrow(large); nrow(FAs); nrow(small.colon); nrow(large.colon); nrow(large.rectal)

# Join small and large together for pooled questionnaire model
common.vars <- c("Cncr_Caco_Clrt", "Qe_Energy", "L_School", "Smoke_Stat", "Match_Caseset", 
                 "Wcrf_C_Cal", "Height_C")
small0 <- select(small, common.vars)
large0 <- select(large, common.vars)
all <- bind_rows(small0, large0)

small.col <- select(small.colon, common.vars)
large.col <- select(large.colon, common.vars)
all.colon <- bind_rows(small.col, large.col)
large.rec <- select(large.rectal, common.vars)

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

# Biocrates colon only
fit7 <- clogit(update(base, ~. + score.1.comps), data = small.colon)
fit8 <- clogit(update(base, ~. + Wcrf_C_Cal), data = small.colon)
fit9 <- clogit(update(base, ~. + score.1.comps), data = large.colon)
fit10 <- clogit(update(base, ~. + Wcrf_C_Cal), data = large.colon)

# Rectal only
fit11 <- clogit(update(base, ~. + score.1.comps), data = large.rectal)
fit12 <- clogit(update(base, ~. + Wcrf_C_Cal), data = large.rectal)

# Pooled questionnaire (all and colon only)
fit0 <- clogit(update(base, ~. + Wcrf_C_Cal), data = all)
fit0col <- clogit(update(base, ~. + Wcrf_C_Cal), data = all.colon)

# Table for manuscript
library(broom)
t1 <- map_df(list(fit1, fit2, fit3, fit4, fit5, fit6, fit7, fit8, fit9, fit10, fit11, fit12, fit0, fit0col), tidy) %>% 
  filter(str_detect(term, "score.|Wcrf")) %>% select(-(std.error : p.value)) %>%
  mutate_at(.vars = c("estimate", "conf.low", "conf.high"), exp)

t1$model <- c("Bioc CRC small sig", "Bioc CRC small score", "Bioc CRC large sig", "Bioc CRC large score", "CRC FA sig", "CRC FA score",
              "Colon small sig", "Colon small score", "Colon large sig", "Colon large score", "Rectal sig",
              "Rectal score", "Bioc all CRC score", "Bioc all colon score")

# Tables of ORs for scores and signatures
score <- t1 %>% filter(term == "Wcrf_C_Cal")
sig   <- t1 %>% filter(term != "Wcrf_C_Cal")

# Forest plots ----
# Signature only for Biocrates small and large (all subjects in study) and fatty acids small
library(broom)
library(tidyverse)

ll2 <- list(fit1, fit3, fit7, fit9)
t2 <- map_df(ll2, tidy) %>% filter(str_detect(term, "score."))

studies <- data.frame(CC = c("CRC small", "CRC large", "Colon small", "Colon large"),
            c(934, 2282, 850, 1670), metabolites = rep("Endogenous", 4))

par(mar=c(5,4,1,2))
library(metafor)
forest(t2$estimate, ci.lb = t2$conf.low, ci.ub = t2$conf.high, refline = 1, 
       rows = c(2,3,6,7),
       xlab = "Odds ratio (per category increase in score)", pch = 18, 
       transf = exp, psize = 1.5, 
       slab = studies$CC, ilab = studies[, 2], 
       ylim = c(0, 10),
       #xlim = c(-1.2, 1.8),
       ilab.pos = 4, ilab.xpos = c(0.1))
par("usr")

#text(c(-1.2, -0.7, -0.3), 8, c("Case-control", "n", "Metabolite signature"), pos = 4)
#text(1.8, 8, "OR [95% CI]", pos = 2)

# Perform meta-analysis for Biocrates, CRC and colon
ma1 <- rma(estimate, sei = std.error, data=t2, method="FE", subset = 1:2)
ma1$QEp # p-value for heterogeneity
ma1$I2 # I2

ma2 <- rma(estimate, sei = std.error, data=t2, method="FE", subset = 3:4)
ma2$QEp
ma2$I2

addpoly(ma1, row = 1, transf = exp, mlab = "", efac = 2)
addpoly(ma2, row = 5, transf = exp, mlab = "", efac = 2)

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









