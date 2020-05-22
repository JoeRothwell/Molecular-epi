# Model CRC status from WCRF score or signature
# Requires crc1, crc2, controls and PLS models to be prepared from CRC_data_prep.R
source("CRC_get_signatures.R")
load("predicted_score_tables1.Rdata")
# Define basic model. Note: Smoke intensity has too many levels for the by sex analysis

library(survival)
base <- Cncr_Caco_Clrt ~ Qe_Energy + L_School + Smoke_Int + Smoke_Stat + #Smoke_Intensity + 
  Height_C + strata(Match_Caseset)
# Dairy product intake = QGE05

# Biocrates small, large, fatty acids (large fasted subset did not improve OR)
# CRC A: score and signature
fit1 <- clogit(update(base, ~. + comp1), data = crc1.ph)
fit2 <- clogit(update(base, ~. + Wcrf_C_Cal), data = crc1.ph)

# By sex
fit1f <- clogit(update(base, ~. + comp1), data = crc1.ph, subset = Sex == 1)
fit1m <- clogit(update(base, ~. + comp1), data = crc1.ph, subset = Sex == 2)
fit2f <- clogit(update(base, ~. + Wcrf_C_Cal), data = crc2.ph, subset = Sex == 1)
fit2m <- clogit(update(base, ~. + Wcrf_C_Cal), data = crc2.ph, subset = Sex == 2)

# CRC B Score, signature and by sex
fit3 <- clogit(update(base, ~. + comp1), data = crc2.ph)
fit4 <- clogit(update(base, ~. + Wcrf_C_Cal), data = crc2.ph)

# By sex
fit3f <- clogit(update(base, ~. + comp1), data = crc2.ph, subset = Sex == 1)
fit3m <- clogit(update(base, ~. + comp1), data = crc2.ph, subset = Sex == 2)
fit4f <- clogit(update(base, ~. + Wcrf_C_Cal), data = crc2.ph, subset = Sex == 1)
fit4m <- clogit(update(base, ~. + Wcrf_C_Cal), data = crc2.ph, subset = Sex == 2)

# CRC A Fatty acids score, signature and by sex
fit5 <- clogit(update(base, ~. + comp2), data = crc3.ph)
fit6 <- clogit(update(base, ~. + Wcrf_C_Cal), data = crc3.ph)

# By sex
fit5f <- clogit(update(base, ~. + comp2), data = crc3.ph, subset = Sex == 1)
fit5m <- clogit(update(base, ~. + comp2), data = crc3.ph, subset = Sex == 2)
fit6f <- clogit(update(base, ~. + Wcrf_C_Cal), data = crc3.ph, subset = Sex == 1)
fit6m <- clogit(update(base, ~. + Wcrf_C_Cal), data = crc3.ph, subset = Sex == 2)

# By subsite
# Biocrates colon only
fit7 <- clogit(update(base, ~. + comp1), data = col1.ph)
fit8 <- clogit(update(base, ~. + Wcrf_C_Cal), data = col1.ph)
fit9 <- clogit(update(base, ~. + comp1), data = col2.ph)
fit10 <- clogit(update(base, ~. + Wcrf_C_Cal), data = col2.ph)

# Rectal only
fit11 <- clogit(update(base, ~. + comp1), data = rec2.ph)
fit12 <- clogit(update(base, ~. + Wcrf_C_Cal), data = rec2.ph)

# Pooled questionnaire (all and colon only)
fit0 <- clogit(update(base, ~. + Wcrf_C_Cal), data = crc.both)
fit0col <- clogit(update(base, ~. + Wcrf_C_Cal), data = col.both)

ll <- list(fit1, fit2, fit3, fit4, fit5, fit6, fit7, fit8, fit9, fit10, fit11, fit12, fit0, fit0col)
ord <- c(6,5,2,1,4,3,13,8,7,10,9,14,12,11)

# Table for manuscript
library(broom)
t1 <- map_df(ll[ord], tidy) %>% 
  filter(str_detect(term, "comp|Wcrf")) %>% select(-(std.error : p.value)) %>%
  mutate_at(c("estimate", "conf.low", "conf.high"), ~ round(exp(.), 2))

# Model names
models <- c("Colorectal A bioc,sig", "Colorectal A bioc,score", "Colorectal B bioc,sig", 
              "Colorectal B bioc,score", "Colorectal A FA,sig", "Colorectal A FA,score",
              "Colon A,sig", "Colon A,score", "Colon B,sig", "Colon B,score", "Rectal,sig",
              "Rectal B,score", "Colorectal all bioc,score", "Colon all bioc,score")

t1$model <- models[ord]

# Tables of ORs for scores and signatures
t1 %>% separate(model, into = c("Site", "Predictor"), sep = ",")
score <- t1 %>% filter(term == "Wcrf_C_Cal")
sig   <- t1 %>% filter(term != "Wcrf_C_Cal")

# Meta analysis of Biocrates studies for CRC and colon
library(metafor)

t2 <- list(fit1, fit3, fit7, fit9) %>% map_df(tidy) %>% filter(str_detect(term, "score."))
ma.crc <- rma(estimate, sei = std.error, data=t2, method="FE", subset = 1:2)
ma.col <- rma(estimate, sei = std.error, data=t2, method="FE", subset = 3:4)

# p-value for heterogeneity and I squared (percentage)
ma.crc$QEp; ma.crc$I2; ma.col$QEp; ma.col$I2

# Forest plot (removed from manuscript)
# Signature only for Biocrates small and large (all subjects in study) and fatty acids small

studies <- data.frame(CC = c("CRC A", "CRC B", "Colon A", "Colon B"),
            c(934, 2282, 850, 1670), metabolites = rep("Endogenous", 4))

par(mar=c(5,4,1,2))
library(metafor)
forest(t2$estimate, ci.lb = t2$conf.low, ci.ub = t2$conf.high, refline = 1, 
       rows = c(2,3,6,7),
       xlab = "Odds ratio per category increase in score", pch = 15, 
       transf = exp, #psize = 1.5, 
       slab = studies$CC, ilab = studies[, 2], header = T,
       ylim = c(0, 10), ilab.pos = 4, ilab.xpos = c(0.1))

addpoly(ma.crc, row = 1, transf = exp, mlab = "", efac = 3)
addpoly(ma.col, row = 5, transf = exp, mlab = "", efac = 3)
text(-0.5, 1, "Fixed effects meta-analysis of A and B", pos = 4, cex = 0.9)
text(-1.2, 2.5, bquote(paste(I^2," = 0, ",italic(p),"-heterogeneity = 0.32")), pos = 4, cex = 0.9)
#text(c(-1.2, -0.7, -0.3), 8, c("Case-control", "n", "Metabolite signature"), pos = 4)
#text(-1.2, 3, bquote(paste("Fixed effects meta-analysis of A and B\n(p-heterogeneity = 0.32,",I^2," = 0)")), pos = 4, cex = 0.9)


# Subgroup analysis by sex
ll2 <- list(fit1m, fit1f, fit2m, fit2f, fit3m, fit3f, fit4m, fit4f, fit5m, fit5f, fit6m, fit6f)
t2 <- map_df(ll2, tidy) %>% filter(str_detect(term, "score.|Wcrf")) %>% 
  select(-(std.error : p.value)) %>%
  mutate_at(c("estimate", "conf.low", "conf.high"), ~ round(exp(.), 2))

t2$models2 <- c("CRC A male, sig", "CRC A female, sig", "CRC A male, score", "CRC A female, score", 
             "CRC B male, sig", "CRC B female, sig", "CRC B male, score", "CRC B female, score", 
             "CRC FA male, sig", "CRC FA female, sig", "CRC FA male, score", "CRC FA female, score")

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









