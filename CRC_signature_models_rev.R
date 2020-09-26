# Model CRC status from WCRF score or signature
# Requires crc1, crc2, controls and PLS models to be prepared from CRC_data_prep.R
#source("CRC_get_signatures_rev.R")
load("pred_score_tables_rev.Rdata")
library(tidyverse)

# For smoke intensity, categories 8, 9 and 10 are collapsed into other (Smoke_Int)
# Define basic model. Note: need to remove variables for rectal subset

### Revised submission to CGH: crc1 and crc2 now merged ###
# Models for score, LH column of table (fatty acids unchanged for resubmission)
# Models for WCRF score LH column of table (corresponding subsets)

library(survival)
base <- Cncr_Caco_Clrt ~ Qe_Energy + L_School + Smoke_Stat + Smoke_Int + Height_C + strata(Match_Caseset)

fit1 <- clogit(update(base, ~. + Wcrf_C_Cal), data = crc3.ph)
fit3 <- clogit(update(base, ~. + Wcrf_C_Cal), data = crc3f.ph)
fit5 <- clogit(update(base, ~. + Wcrf_C_Cal), data = crc3m.ph)
fit7 <- clogit(update(base, ~. + Wcrf_C_Cal), data = crc.ph)
fit9 <- clogit(update(base, ~. + Wcrf_C_Cal), data = crcF.ph)
fit11 <- clogit(update(base, ~. + Wcrf_C_Cal), data = crcM.ph)
fit13 <- clogit(update(base, ~. + Wcrf_C_Cal), data = col.ph)
fit15 <- clogit(update(base, ~. + Wcrf_C_Cal), data = colF.ph)
fit17 <- clogit(update(base, ~. + Wcrf_C_Cal), data = colM.ph)
fit19 <- clogit(update(base, ~. + Wcrf_C_Cal), data = rec.ph)
fit21 <- clogit(update(base, ~. - Smoke_Stat + Wcrf_C_Cal), data = recF.ph)
fit23 <- clogit(update(base, ~. - L_School + Wcrf_C_Cal), data = recM.ph)

# Models for signature, RH column of table (fatty acids unchanged for resubmission)
fit2 <- clogit(update(base, ~. + comp2), data = crc3.ph)
fit4 <- clogit(update(base, ~. + comp2), data = crc3f.ph)
fit6 <- clogit(update(base, ~. + comp2), data = crc3m.ph)
fit8 <- clogit(update(base, ~. + comp1), data = crc.ph)
fit10 <- clogit(update(base, ~. + comp1), data = crcF.ph)
fit12 <- clogit(update(base, ~. + comp1), data = crcM.ph)
fit14 <- clogit(update(base, ~. + comp1), data = col.ph)
fit16 <- clogit(update(base, ~. + comp1), data = colF.ph)
fit18 <- clogit(update(base, ~. + comp1), data = colM.ph)
fit20 <- clogit(update(base, ~. + comp1), data = rec.ph)
fit22 <- clogit(update(base, ~. - Smoke_Stat + comp1), data = recF.ph)
fit24 <- clogit(update(base, ~. - Smoke_Stat - L_School + comp1), data = recM.ph)


# Put score and signature models in separate lists (separate columns in table)
scomodlist <- list(fit1, fit3, fit5, fit7, fit9, fit11, fit13, fit15, fit17, fit19, fit21, fit23)
sigmodlist <- list(fit2, fit4, fit6, fit8, fit10, fit12, fit14, fit16, fit18, fit20, fit22, fit24)

scorenames <- c("CRC FA, all", "CRC FA, women", "CRC FA men", 
                "CRC endogenous all", "CRC endogenous women", "CRC endogenous men", 
                "Colon endogenous all", "Colon endogenous women", "Colon endogenous men", 
                "Rectal endogenous, all", "Rectal endogenous, women", "Rectal endogenous, men")

library(broom)
scomods <- map_df(scomodlist, ~tidy(., exponentiate = T)) %>% filter(term == "Wcrf_C_Cal") %>%
  mutate_if(is.numeric, ~round(., 2)) %>% unite(OR.CI, estimate, conf.low, conf.high, sep = "-") %>%
  add_column(model = scorenames, .before = T)

sigmods <- map_df(sigmodlist, ~tidy(., exponentiate = T)) %>% filter(term == "comp2" | term == "comp1") %>%
  mutate_if(is.numeric, ~round(., 2)) %>% unite(OR.CI, estimate, conf.low, conf.high, sep = "-")# %>%
  #add_column(model = scorenames[-c(10, 13)], .before = T)


# Additional models for revision: colon proximal and distal by sex, colon and distal for fatty acids
fit1 <- clogit(update(base, ~. + comp1), data = prox.ph)
fit2 <- clogit(update(base, ~. + comp1), data = proxF.ph)
fit3 <- clogit(update(base, ~. + comp1), data = proxM.ph)

fit4 <- clogit(update(base, ~. + comp1), data = dist.ph)
fit5 <- clogit(update(base, ~. + comp1), data = distF.ph)
fit6 <- clogit(update(base, ~. + comp1), data = distM.ph)

fit7 <- clogit(update(base, ~. + comp2), data = col3.ph)
fit8 <- clogit(update(base, ~. + comp2), data = col3f.ph)
fit9 <- clogit(update(base, ~. + comp2), data = col3m.ph)

fit10 <- clogit(update(base, ~. + comp2), data = dist3.ph)
fit11 <- clogit(update(base, ~. + comp2), data = prox3.ph)

modlist <- list(fit1, fit2, fit3, fit4, fit5, fit6, fit7, fit8, fit9, fit10, fit11)

library(broom)
mods <- map_df(modlist, ~tidy(., exponentiate = T)) %>% filter(term == "comp2" | term == "comp1") %>%
  mutate_if(is.numeric, ~round(., 2)) %>% unite(OR.CI, estimate, conf.low, conf.high, sep = "-")


# Score
fit1 <- clogit(update(base, ~. + Wcrf_C_Cal), data = prox.ph)
fit2 <- clogit(update(base, ~. + Wcrf_C_Cal), data = proxF.ph)
fit3 <- clogit(update(base, ~. + Wcrf_C_Cal), data = proxM.ph)

fit4 <- clogit(update(base, ~. + Wcrf_C_Cal), data = dist.ph)
fit5 <- clogit(update(base, ~. + Wcrf_C_Cal), data = distF.ph)
fit6 <- clogit(update(base, ~. + Wcrf_C_Cal), data = distM.ph)

fit7 <- clogit(update(base, ~. + Wcrf_C_Cal), data = col3.ph)
fit8 <- clogit(update(base, ~. + Wcrf_C_Cal), data = col3f.ph)
fit9 <- clogit(update(base, ~. + Wcrf_C_Cal), data = col3m.ph)

fit10 <- clogit(update(base, ~. + Wcrf_C_Cal), data = dist3.ph)
fit11 <- clogit(update(base, ~. + Wcrf_C_Cal), data = prox3.ph)

modlistA <- list(fit1, fit2, fit3, fit4, fit5, fit6, fit7, fit8, fit9, fit10, fit11)

library(broom)
modsA <- map_df(modlistA, ~tidy(., exponentiate = T)) %>% filter(term == "Wcrf_C_Cal") %>%
  mutate_if(is.numeric, ~round(., 2)) %>% unite(OR.CI, estimate, conf.low, conf.high, sep = "-")

modscoresig <- cbind(modsA, mods)










