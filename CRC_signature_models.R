# Find metabolic and fatty acids signatures of WCRF score
# Requires crc1, crc2, controls and PLS models to be prepared from CRC_data_prep.R
source("CRC_get_signatures.R")

# Functions to get signature-derived WCRF scores for case-controls and pooled controls
predict.scores <- function(crc, dat, mod, bioc = T, scorecomp.only = F){

  library(tidyverse)
  library(zoo)
  
  print(paste("Subjects in case-control: ", nrow(crc)))
  
  # Put CRC CC compounds in same order as in controls dataset (fatty acids, biocrates)
  if(bioc == F)  cols <- common.cols else cols <- colnames(dat)
  cmpds <- crc %>% select(one_of(cols))

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

# Predict scores for different datasets
pred.crc1 <- predict.scores(crc1, ctrlA, mod1a)
pred.crc2 <- predict.scores(crc2, ctrlB, mod1b)
pred.fa   <- predict.scores(crc1fa, concs, mod2, bioc = F)
pred.colon1 <- predict.scores(colon1, ctrlA, mod1a)
pred.colon2 <- predict.scores(colon2, ctrlB, mod1b)
pred.rectal <- predict.scores(rectal2, ctrlB, mod1b)

# Predict sex-specific scores (not used)
pred.crc1.m <- predict.scores(crc1.ma, ctrlAm, mod1m)
pred.crc1.f <- predict.scores(crc1.fe, ctrlAf, mod1f)
pred.crc2.m <- predict.scores(crc2.ma, ctrlAm, mod1m)
pred.crc2.f <- predict.scores(crc2.fe, ctrlAf, mod1f)
pred.FAf <- predict.scores(crc1faf, concs.f, modFAf, bioc = F)

nrow(pred.crc1); nrow(pred.crc2); nrow(pred.fa); nrow(pred.colon1); nrow(pred.colon2); nrow(pred.rectal)

# Join small and large together for pooled questionnaire model
common.vars <- c("Cncr_Caco_Clrt", "Qe_Energy", "L_School", "Smoke_Intensity", "Match_Caseset", 
                 "Wcrf_C_Cal", "Height_C")
pred.crc1.vars <- select(pred.crc1, all_of(common.vars))
pred.crc2.vars <- select(pred.crc2, all_of(common.vars))
pred.all       <- bind_rows(pred.crc1.vars, pred.crc2.vars)

pred.colon1.vars <- select(pred.colon1, all_of(common.vars))
pred.colon2.vars <- select(pred.colon2, all_of(common.vars))
pred.colon.all  <- bind_rows(pred.colon1.vars, pred.colon2.vars)
pred.rectal.vars   <- select(pred.rectal, all_of(common.vars))

# Model CC status from calculated or signature-predicted score for four datasets

# Load predicted score tables for CRC
load("predicted_score_tables.Rdata")

library(survival)
# Note: Smoke intensity has too many levels for the by sex analysis
base <- Cncr_Caco_Clrt ~ Qe_Energy + L_School + Smoke_Stat + #Smoke_Intensity + 
  Height_C + strata(Match_Caseset)

# Biocrates small, large, fatty acids (large fasted subset did not improve OR)
# CRC A: score and signature
fit1 <- clogit(update(base, ~. + score.1.comps), data = pred.crc1)
fit2 <- clogit(update(base, ~. + Wcrf_C_Cal), data = pred.crc1)
# By sex
fit1m <- clogit(update(base, ~. + score.1.comps), data = pred.crc1, subset = Sex == 1)
fit1f <- clogit(update(base, ~. + score.1.comps), data = pred.crc1, subset = Sex == 2)
fit2m <- clogit(update(base, ~. + Wcrf_C_Cal), data = pred.crc1, subset = Sex == 1)
fit2f <- clogit(update(base, ~. + Wcrf_C_Cal), data = pred.crc1, subset = Sex == 2)
# By sex, sex-specific signatures
fit1ms <- clogit(update(base, ~. + score.1.comps), data = pred.crc1.m)
fit1fs <- clogit(update(base, ~. + score.1.comps), data = pred.crc1.f)


# CRC B Score, signature and by sex
fit3 <- clogit(update(base, ~. + score.1.comps), data = pred.crc2)
fit4 <- clogit(update(base, ~. + Wcrf_C_Cal), data = pred.crc2)
# By sex
fit3m <- clogit(update(base, ~. + score.1.comps), data = pred.crc2, subset = Sex == 1)
fit3f <- clogit(update(base, ~. + score.1.comps), data = pred.crc2, subset = Sex == 2)
fit4m <- clogit(update(base, ~. + Wcrf_C_Cal), data = pred.crc2, subset = Sex == 1)
fit4f <- clogit(update(base, ~. + Wcrf_C_Cal), data = pred.crc2, subset = Sex == 2)
# By sex, sex-specific signatures
fit3ms <- clogit(update(base, ~. + score.1.comps), data = pred.crc2.m)
fit3fs <- clogit(update(base, ~. + score.1.comps), data = pred.crc2.f)


# CRC A Fatty acids score, signature and by sex
fit5 <- clogit(update(base, ~. + score.2.comps), data = pred.fa)
fit6 <- clogit(update(base, ~. + Wcrf_C_Cal), data = pred.fa)
# By sex
fit5m <- clogit(update(base, ~. + score.2.comps), data = pred.fa, subset = Sex == 1)
fit5f <- clogit(update(base, ~. + score.2.comps), data = pred.fa, subset = Sex == 2)
fit6m <- clogit(update(base, ~. + Wcrf_C_Cal), data = pred.fa, subset = Sex == 1)
fit6f <- clogit(update(base, ~. + Wcrf_C_Cal), data = pred.fa, subset = Sex == 2)
# By sex, sex-specific signature (female only)
fit5fs <- clogit(update(base, ~. + score.2.comps), data = pred.FAf)


# By subsite
# Biocrates colon only
fit7 <- clogit(update(base, ~. + score.1.comps), data = pred.colon1)
fit8 <- clogit(update(base, ~. + Wcrf_C_Cal), data = pred.colon1)
fit9 <- clogit(update(base, ~. + score.1.comps), data = pred.colon2)
fit10 <- clogit(update(base, ~. + Wcrf_C_Cal), data = pred.colon2)

# Rectal only
fit11 <- clogit(update(base, ~. + score.1.comps), data = pred.rectal)
fit12 <- clogit(update(base, ~. + Wcrf_C_Cal), data = pred.rectal)

# Pooled questionnaire (all and colon only)
fit0 <- clogit(update(base, ~. + Wcrf_C_Cal), data = pred.all)
fit0col <- clogit(update(base, ~. + Wcrf_C_Cal), data = pred.colon.all)

# Table for manuscript
library(broom)
t1 <- map_df(list(fit1, fit2, fit3, fit4, fit5, fit6, fit7, fit8, fit9, fit10, fit11, fit12, fit0, fit0col), tidy) %>% 
  filter(str_detect(term, "score.|Wcrf")) %>% select(-(std.error : p.value)) %>%
  #mutate_at(.vars = c("estimate", "conf.low", "conf.high"), exp)
  mutate_at(c("estimate", "conf.low", "conf.high"), ~ round(exp(.), 2))

t1$model <- c("Colorectal A bioc sig", "Colorectal A bioc score", "Colorectal B bioc sig", 
              "Colorectal B bioc score", "Colorectal A FA sig", "Colorectal A FA score",
              "Colon A sig", "Colon A score", "Colon B sig", "Colon B score", "Rectal sig",
              "Rectal B score", "Colorectal all bioc score", "Colon all bioc score")

# Tables of ORs for scores and signatures
score <- t1 %>% filter(term == "Wcrf_C_Cal")
sig   <- t1 %>% filter(term != "Wcrf_C_Cal")

# Meta analysis of Biocrates studies
library(broom)
library(tidyverse)

ll2 <- list(fit1m, fit1f, fit3m, fit3f, fit5m, fit5f)
t2 <- map_df(ll2, tidy) %>% filter(str_detect(term, "score."))

# Perform meta-analysis for Biocrates, CRC and colon
ma1 <- rma(estimate, sei = std.error, data=t2, method="FE", subset = 1:2)
ma1$QEp # p-value for heterogeneity
ma1$I2 # I2

ma2 <- rma(estimate, sei = std.error, data=t2, method="FE", subset = 3:4)
ma2$QEp
ma2$I2

# Forest plot (removed from manuscript)
# Signature only for Biocrates small and large (all subjects in study) and fatty acids small

studies <- data.frame(CC = c("CRC A", "CRC B", "Colon A", "Colon B"),
            c(934, 2282, 850, 1670), metabolites = rep("Endogenous", 4))

par(mar=c(5,4,1,2))
library(metafor)
forest(t2$estimate, ci.lb = t2$conf.low, ci.ub = t2$conf.high, refline = 1, 
       rows = c(2,3,6,7),
       xlab = "Odds ratio per category increase in score", #pch = 18, 
       transf = exp, #psize = 1.5, 
       slab = studies$CC, ilab = studies[, 2], 
       ylim = c(0, 10), ilab.pos = 4, ilab.xpos = c(0.1))
par("usr")

addpoly(ma1, row = 1, transf = exp, mlab = "", efac = 2)
addpoly(ma2, row = 5, transf = exp, mlab = "", efac = 2)
text(-1.2, 3.5, "Fixed effects meta-analysis of A and B", pos = 4, cex = 0.9)
text(-1.2, 2.5, bquote(paste(I^2," = 0, ",italic(p),"-heterogeneity = 0.32")), pos = 4, cex = 0.9)
#text(c(-1.2, -0.7, -0.3), 8, c("Case-control", "n", "Metabolite signature"), pos = 4)
#text(1.8, 8, "OR [95% CI]", pos = 2)
#text(-1.2, 3, bquote(paste("Fixed effects meta-analysis of A and B\n(p-heterogeneity = 0.32,",I^2," = 0)")), pos = 4, cex = 0.9)


# Test for heterogeneity by sex for signature 
ll2 <- list(fit1, fit3, fit7, fit9)
t2 <- map_df(ll2, tidy) %>% filter(str_detect(term, "score."))

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









