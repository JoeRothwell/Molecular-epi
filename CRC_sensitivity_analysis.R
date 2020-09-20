load("pred_score_tables_rev.Rdata")
library(survival)
library(broom)
base <- Cncr_Caco_Clrt ~ Qe_Energy + 
  L_School + Smoke_Stat + Smoke_Int + Height_C + strata(Match_Caseset)

### CRC1 and 2 merged for revision ###

# For smoke intensity, categories 8, 9 and 10 are collapsed into other (Smoke_Int)
# Replace smoking duration NAs with 0
crc.ph$Dur_Smok[is.na(crc.ph$Dur_Smok)] <- 0
col.ph$Dur_Smok[is.na(col.ph$Dur_Smok)] <- 0
rec.ph$Dur_Smok[is.na(rec.ph$Dur_Smok)] <- 0
crc3.ph$Dur_Smok[is.na(crc3.ph$Dur_Smok)] <- 0

# Dairy product intake = QGE05

# CRC A Fatty acids score, signature and by sex
fit5 <- clogit(update(base, ~. + comp2), data = crc3.ph)
fit6 <- clogit(update(base, ~. + Wcrf_C_Cal), data = crc3.ph)

fit5si <- clogit(update(base, ~. + comp2 + Smoke_Int), data = crc3.ph)
fit6si <- clogit(update(base, ~. + Wcrf_C_Cal + Smoke_Int), data = crc3.ph)

fit5sd <- clogit(update(base, ~. + comp2 + Dur_Smok), data = crc3.ph)
fit6sd <- clogit(update(base, ~. + Wcrf_C_Cal + Dur_Smok), data = crc3.ph)

fit5dp <- clogit(update(base, ~. + comp2 + Qge05), data = crc3.ph)
fit6dp <- clogit(update(base, ~. + Wcrf_C_Cal + Qge05), data = crc3.ph)

fit5t <- clogit(update(base, ~. + comp2), data = crc3t.ph)
fit6t <- clogit(update(base, ~. + Wcrf_C_Cal), data = crc3t.ph)

# Extra for revision: Signature model adjusted for score
fit1 <- clogit(update(base, ~. + Wcrf_C_Cal + comp2), data = crc3.ph)
summary(fit)

sigmodsC <- map_df(list(fit5, fit5si, fit5sd, fit5dp, fit5t), ~tidy(., exponentiate = T)) %>% filter(term == "comp2")
scomodsC <- map_df(list(fit6, fit6si, fit6sd, fit6dp, fit6t), ~tidy(., exponentiate = T)) %>% filter(term == "Wcrf_C_Cal")


# CRC Biocrates
fit1 <- clogit(update(base, ~. + comp1), data = crc.ph)
fit2 <- clogit(update(base, ~. + Wcrf_C_Cal), data = crc.ph)

fit1si <- clogit(update(base, ~. + comp1 + Smoke_Int), data = crc.ph)
fit2si <- clogit(update(base, ~. + Wcrf_C_Cal + Smoke_Int), data = crc.ph)

fit1sd <- clogit(update(base, ~. + comp1 + Dur_Smok), data = crc.ph)
fit2sd <- clogit(update(base, ~. + Wcrf_C_Cal + Dur_Smok), data = crc.ph)

fit1dp <- clogit(update(base, ~. + comp1 + Qge05), data = crc.ph)
fit2dp <- clogit(update(base, ~. + Wcrf_C_Cal + Qge05), data = crc.ph)

fit1t <- clogit(update(base, ~. + comp1), data = crcT.ph)
fit2t <- clogit(update(base, ~. + Wcrf_C_Cal), data = crcT.ph)

# Signature model adjusted for score
fit1 <- clogit(update(base, ~. + Wcrf_C_Cal + comp1), data = crc.ph)
summary(fit)

sigmodsA <- map_df(list(fit1, fit1si, fit1sd, fit1dp, fit1t), ~tidy(., exponentiate = T)) %>% filter(term == "comp1")
scomodsA <- map_df(list(fit2, fit2si, fit2sd, fit2dp, fit2t), ~tidy(., exponentiate = T)) %>% filter(term == "Wcrf_C_Cal")


# Colon Biocrates
fit3 <- clogit(update(base, ~. + comp1), data = col.ph)
fit4 <- clogit(update(base, ~. + Wcrf_C_Cal), data = col.ph)

fit3si <- clogit(update(base, ~. + comp1 + Smoke_Int), data = col.ph)
fit4si <- clogit(update(base, ~. + Wcrf_C_Cal + Smoke_Int), data = col.ph)

fit3sd <- clogit(update(base, ~. + comp1 + Dur_Smok), data = col.ph)
fit4sd <- clogit(update(base, ~. + Wcrf_C_Cal + Dur_Smok), data = col.ph)

fit3dp <- clogit(update(base, ~. + comp1 + Qge05), data = col.ph)
fit4dp <- clogit(update(base, ~. + Wcrf_C_Cal + Qge05), data = col.ph)

fit3t <- clogit(update(base, ~. + comp1), data = colT.ph)
fit4t <- clogit(update(base, ~. + Wcrf_C_Cal), data = colT.ph)

# Signature model adjusted for score
fit1 <- clogit(update(base, ~. + Wcrf_C_Cal + comp1), data = col.ph)
summary(fit1)

sigmodsA <- map_df(list(fit3, fit3si, fit3sd, fit3dp, fit3t), ~tidy(., exponentiate = T)) %>% filter(term == "comp1")
scomodsA <- map_df(list(fit4, fit4si, fit4sd, fit4dp, fit4t), ~tidy(., exponentiate = T)) %>% filter(term == "Wcrf_C_Cal")


# Rectal Biocrates
fit7 <- clogit(update(base, ~. + comp1), data = rec.ph)
fit8 <- clogit(update(base, ~. + Wcrf_C_Cal), data = rec.ph)

fit7si <- clogit(update(base, ~. + comp1 + Smoke_Int), data = rec.ph)
fit8si <- clogit(update(base, ~. + Wcrf_C_Cal + Smoke_Int), data = rec.ph)

fit7sd <- clogit(update(base, ~. + comp1 + Dur_Smok), data = rec.ph)
fit8sd <- clogit(update(base, ~. + Wcrf_C_Cal + Dur_Smok), data = rec.ph)

fit7dp <- clogit(update(base, ~. + comp1 + Qge05), data = rec.ph)
fit8dp <- clogit(update(base, ~. + Wcrf_C_Cal + Qge05), data = rec.ph)

fit7t <- clogit(update(base, ~. + comp1), data = recT.ph)
fit8t <- clogit(update(base, ~. + Wcrf_C_Cal), data = recT.ph)

fit1 <- clogit(update(base, ~. + Wcrf_C_Cal + comp1), data = rec.ph)
summary(fit1)

sigmodsD <- map_df(list(fit7, fit7si, fit7sd, fit7dp, fit7t), ~tidy(., exponentiate = T)) %>% filter(term == "comp1")
scomodsD <- map_df(list(fit8, fit8si, fit8sd, fit8dp, fit8t), ~tidy(., exponentiate = T)) %>% filter(term == "Wcrf_C_Cal")


# Make final table to copy and paste into manuscript
sigmods <- rbind(sigmodsC, sigmodsA, sigmodsB) %>% mutate_if(is.numeric, ~round(., 2)) %>% 
  unite(OR.CI, estimate, conf.low, conf.high, sep = "-")

scoremods <- rbind(scomodsC, scomodsA, scomodsB) %>% mutate_if(is.numeric, ~round(., 2)) %>%
  unite(OR.CI, estimate, conf.low, conf.high, sep = "-")


# By subsite only for supplemental data

base <- Cncr_Caco_Clrt ~ Qe_Energy + L_School + Smoke_Int + Height_C + strata(Match_Caseset)

fit2 <- clogit(update(base, ~. + Wcrf_C_Cal), data = crc.ph)

# Fatty acids
fit1 <- clogit(update(base, ~. + Wcrf_C_Cal), data = crc3.ph)
fit2 <- clogit(update(base, ~. + Wcrf_C_Cal), data = prox3.ph)
fit3 <- clogit(update(base, ~. + Wcrf_C_Cal), data = dist3.ph)

# Biocrates
fit4 <- clogit(update(base, ~. + Wcrf_C_Cal), data = crc.ph)
fit5 <- clogit(update(base, ~. + Wcrf_C_Cal), data = col.ph)
fit6 <- clogit(update(base, ~. + Wcrf_C_Cal), data = prox.ph)
fit7 <- clogit(update(base, ~. + Wcrf_C_Cal), data = dist.ph)
fit8 <- clogit(update(base, ~. + Wcrf_C_Cal), data = rec.ph)


# Models for signature, RH column of table (fatty acids unchanged for resubmission)
# Fatty acids
fit1 <- clogit(update(base, ~. + comp2), data = crc3.ph)
fit2 <- clogit(update(base, ~. + comp2), data = prox3.ph)
fit3 <- clogit(update(base, ~. + comp2), data = dist3.ph)

# Biocrates
fit4 <- clogit(update(base, ~. + Wcrf_C_Cal), data = crc.ph)
fit5 <- clogit(update(base, ~. + Wcrf_C_Cal), data = col.ph)
fit6 <- clogit(update(base, ~. + Wcrf_C_Cal), data = prox.ph)
fit7 <- clogit(update(base, ~. + Wcrf_C_Cal), data = dist.ph)
fit8 <- clogit(update(base, ~. + Wcrf_C_Cal), data = rec.ph)

# By score of individual components: colorectal, colon, rectal
library(broom)
library(tidyverse)
fit1 <- clogit(update(base, ~. + Wcrf_C_Cal), data = crc3.ph)


# Variable names: Wcrf_Bmi", "Wcrf_Pa", "Wcrf_Fwg_Cal", "Wcrf_Pf_Cal", "Wcrf_Meat_Cal", "Wcrf_Alc", "Wcrf_C_Cal"
fit1 <- clogit(update(base, ~. + Wcrf_Bmi), data = crc.ph)
fit2 <- clogit(update(base, ~. + Wcrf_Pa), data = crc.ph)
fit3 <- clogit(update(base, ~. + Wcrf_Fwg_Cal), data = crc.ph)
fit4 <- clogit(update(base, ~. + Wcrf_Pf_Cal), data = crc.ph)
fit5 <- clogit(update(base, ~. + Wcrf_Meat_Cal), data = crc.ph)
fit6 <- clogit(update(base, ~. + Wcrf_Alc), data = crc.ph)

ll <- list(fit1, fit2, fit3, fit4, fit5, fit6)
mods <- map_df(ll, ~tidy(., exponentiate = T)) %>% filter(grepl("Wcrf_", term)) %>%
  mutate_if(is.numeric, ~round(., 2))

fit1 <- clogit(update(base, ~. + Wcrf_Bmi), data = rec.ph)
fit2 <- clogit(update(base, ~. + Wcrf_Pa), data = rec.ph)
fit3 <- clogit(update(base, ~. + Wcrf_Fwg_Cal), data = rec.ph)
fit4 <- clogit(update(base, ~. + Wcrf_Pf_Cal), data = rec.ph)
fit5 <- clogit(update(base, ~. + Wcrf_Meat_Cal), data = rec.ph)
fit6 <- clogit(update(base, ~. + Wcrf_Alc), data = rec.ph)

ll <- list(fit1, fit2, fit3, fit4, fit5, fit6)
mods <- map_df(ll, ~tidy(., exponentiate = T)) %>% filter(grepl("Wcrf_", term)) %>%
  mutate_if(is.numeric, ~round(., 2))



# Metabolite-CRC associations for manuscript Table 2
# First need to run CRC_data_prep to get crc and crc3 objects
library(ggplot)
# Biocrates compounds. Split into quartiles with cut_number
df2 <- crc[, cols] %>%
  select(Glyceroph_Lysopc_A_C17_0, Glyceroph_Pc_Ae_C40_6, Glyceroph_Pc_Ae_C36_2, Glyceroph_Pc_Ae_C38_2, Aminoacid_Ser, 
         Glyceroph_Pc_Ae_C40_6, Glyceroph_Lysopc_A_C18_2, Aminoacid_Gly, Glyceroph_Pc_Ae_C40_3, Glyceroph_Pc_Aa_C32_1, 
         Glyceroph_Pc_Aa_C38_4, Glyceroph_Pc_Aa_C36_4, Aminoacid_Glu, Glyceroph_Pc_Aa_C34_4, Glyceroph_Pc_Aa_C40_4, 
         Glyceroph_Pc_Ae_C38_3) %>%
  mutate_all(~cut_number(., n = 4, labels = 1:4))

# Define function to apply across quartiles
clrfun <- function(x)  { clogit(Cncr_Caco_Clrt ~ x + Qe_Energy + L_School + Smoke_Stat + Smoke_Int + Height_C + lab +
                               strata(Match_Caseset), data = crc) }
fits <- apply(df2, 2, clrfun)
mods1 <- map_df(fits, ~tidy(., exponentiate = T)) %>% filter(grepl("x4", term)) %>%
  mutate_if(is.numeric, ~round(., 2))

# Fatty acids. Modelled continuously because few values.
df5 <- crc3[, colnames(ctrlC)] %>% 
  select(P17_0, P15_0, P22_5n_6, P18_1n_9c, P16_1n_7c_n_9c, P22_1n_9, P20_3n_9, P24_1n_9) %>% scale()

clrfun1 <- function(x)  { clogit(Cncr_Caco_Clrt ~ x + Qe_Energy + L_School + Smoke_Stat + Smoke_Int + Height_C + 
                                  strata(Match_Caseset), data = crc3) }
fits <- apply(df5, 2, clrfun1)
mods <- map_df(fits, ~tidy(., exponentiate = T)) %>% filter(grepl("x", term)) %>% mutate_if(is.numeric, ~round(., 2))
