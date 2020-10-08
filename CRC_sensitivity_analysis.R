library(survival)
library(broom)
load("predscore_df_subsite.Rdata")
base <- Cncr_Caco_Clrt ~ Qe_Energy + L_School + Smoke_Stat +  Height_C + strata(Match_Caseset)

### CRC1 and 2 merged for revision ###

# For smoke intensity, categories 8, 9 and 10 are collapsed into other (Smoke_Int)
# Replace smoking duration NAs with 0
crc.ph$Dur_Smok[is.na(crc.ph$Dur_Smok)] <- 0
col.ph$Dur_Smok[is.na(col.ph$Dur_Smok)] <- 0
rec.ph$Dur_Smok[is.na(rec.ph$Dur_Smok)] <- 0
crc3.ph$Dur_Smok[is.na(crc3.ph$Dur_Smok)] <- 0

# Dairy product intake = QGE05

# CRC A Fatty acids score, signature and by sex
fit1 <- clogit(update(base, ~. + comp2), data = crc3.ph)
fit2 <- clogit(update(base, ~. + Wcrf_C_Cal), data = crc3.ph)
fit3 <- clogit(update(base, ~. + comp2 + Smoke_Int), data = crc3.ph)
fit4 <- clogit(update(base, ~. + Wcrf_C_Cal + Smoke_Int), data = crc3.ph)
fit5 <- clogit(update(base, ~. + comp2 + Dur_Smok), data = crc3.ph)
fit6 <- clogit(update(base, ~. + Wcrf_C_Cal + Dur_Smok), data = crc3.ph)
fit7 <- clogit(update(base, ~. + comp2 + Qge05), data = crc3.ph)
fit8 <- clogit(update(base, ~. + Wcrf_C_Cal + Qge05), data = crc3.ph)
fit9 <- clogit(update(base, ~. + comp2), data = crc3t.ph)
fit10 <- clogit(update(base, ~. + Wcrf_C_Cal), data = crc3t.ph)

# Extra for revision: Signature model adjusted for score
fit11 <- clogit(update(base, ~. + Wcrf_C_Cal + comp2), data = crc3.ph)

# Fatty acids
fit12 <- clogit(update(base, ~. + Smoke_Int + comp2), data = crc3N.ph)
fit13 <- clogit(update(base, ~. + Smoke_Int + comp2), data = crc3O.ph)
fit14 <- clogit(update(base, ~. + Smoke_Int + comp2), data = crc3L.ph)
fit15 <- clogit(update(base, ~. + Smoke_Int + comp2), data = crc3H.ph)

sigmodsC <- map_df(list(fit1, fit3, fit5, fit7, fit11, fit12, fit13, fit14, fit15, fit9), ~tidy(., exponentiate = T)) %>% 
  filter(term == "comp2") %>% mutate_if(is.numeric, ~round(., 2)) %>% 
  unite(OR.CI, estimate, conf.low, conf.high, sep = "-") 

scomodsC <- map_df(list(fit2, fit4, fit6, fit8, fit10), ~tidy(., exponentiate = T)) %>% filter(term == "Wcrf_C_Cal") %>% 
  mutate_if(is.numeric, ~round(., 2)) %>% unite(OR.CI, estimate, conf.low, conf.high, sep = "-")

rm(list = ls(pattern = "fit"))
# CRC Biocrates
fit1 <- clogit(update(base, ~. + lab + comp1), data = crc.ph)
fit2 <- clogit(update(base, ~. + Wcrf_C_Cal), data = crc.ph)
fit3 <- clogit(update(base, ~. + lab + comp1 + Smoke_Int), data = crc.ph)
fit4 <- clogit(update(base, ~. + Wcrf_C_Cal + Smoke_Int), data = crc.ph)
fit5 <- clogit(update(base, ~. + lab + comp1 + Dur_Smok), data = crc.ph)
fit6 <- clogit(update(base, ~. + Wcrf_C_Cal + Dur_Smok), data = crc.ph)
fit7 <- clogit(update(base, ~. + lab + comp1 + Qge05), data = crc.ph)
fit8 <- clogit(update(base, ~. + Wcrf_C_Cal + Qge05), data = crc.ph)
fit9 <- clogit(update(base, ~. + lab + comp1), data = crcT.ph)
fit10 <- clogit(update(base, ~. + Wcrf_C_Cal), data = crcT.ph)

# Signature model adjusted for score
fit11 <- clogit(update(base, ~. + Wcrf_C_Cal + comp1), data = crc.ph)

# Models for high and low BMI and WCRF score
fit12 <- clogit(update(base, ~. + lab + comp1), data = crcN.ph)
fit13 <- clogit(update(base, ~. + lab + comp1), data = crcO.ph)
fit14 <- clogit(update(base, ~. + lab + comp1), data = crcL.ph)
fit15 <- clogit(update(base, ~. + lab + comp1), data = crcH.ph)

# Summarise in tables
sigmodsA <- map_df(list(fit1, fit3, fit5, fit7, fit11, fit12, fit13, fit14, fit15, fit9), ~tidy(., exponentiate = T)) %>% 
  filter(term == "comp1") %>% mutate_if(is.numeric, ~round(., 2)) %>% 
  unite(OR.CI, estimate, conf.low, conf.high, sep = "-")  

scomodsA <- map_df(list(fit2, fit4, fit6, fit8, fit10), ~tidy(., exponentiate = T)) %>% filter(term == "Wcrf_C_Cal") %>% 
  mutate_if(is.numeric, ~round(., 2)) %>% unite(OR.CI, estimate, conf.low, conf.high, sep = "-") 

rm(list = ls(pattern = "fit"))
### Colon Biocrates
fit1 <- clogit(update(base, ~. + lab + comp1), data = col.ph)
fit2 <- clogit(update(base, ~. + Wcrf_C_Cal), data = col.ph)
fit3 <- clogit(update(base, ~. + lab + comp1 + Smoke_Int), data = col.ph)
fit4 <- clogit(update(base, ~. + Wcrf_C_Cal + Smoke_Int), data = col.ph)
fit5 <- clogit(update(base, ~. + lab + comp1 + Dur_Smok), data = col.ph)
fit6 <- clogit(update(base, ~. + Wcrf_C_Cal + Dur_Smok), data = col.ph)
fit7 <- clogit(update(base, ~. + lab + comp1 + Qge05), data = col.ph)
fit8 <- clogit(update(base, ~. + Wcrf_C_Cal + Qge05), data = col.ph)
fit9 <- clogit(update(base, ~. + lab + comp1), data = colT.ph)
fit10 <- clogit(update(base, ~. + Wcrf_C_Cal), data = colT.ph)

# Signature model adjusted for score
fit11 <- clogit(update(base, ~. + lab + Wcrf_C_Cal + comp1), data = col.ph)

sigmodsA <- map_df(list(fit1, fit3, fit5, fit7, fit11, fit10), ~tidy(., exponentiate = T)) %>% filter(term == "comp1") %>% 
  mutate_if(is.numeric, ~round(., 2)) %>% unite(OR.CI, estimate, conf.low, conf.high, sep = "-")

scomodsA <- map_df(list(fit2, fit4, fit6, fit8, fit10), ~tidy(., exponentiate = T)) %>% filter(term == "Wcrf_C_Cal") %>% 
  mutate_if(is.numeric, ~round(., 2)) %>% unite(OR.CI, estimate, conf.low, conf.high, sep = "-")

rm(list = ls(pattern = "fit"))
# Rectal Biocrates
fit1 <- clogit(update(base, ~. + lab + comp1), data = rec.ph)
fit2 <- clogit(update(base, ~. + Wcrf_C_Cal), data = rec.ph)
fit3 <- clogit(update(base, ~. + lab + comp1 + Smoke_Int), data = rec.ph)
fit4 <- clogit(update(base, ~. + Wcrf_C_Cal + Smoke_Int), data = rec.ph)
fit5 <- clogit(update(base, ~. + lab + comp1 + Dur_Smok), data = rec.ph)
fit6 <- clogit(update(base, ~. + Wcrf_C_Cal + Dur_Smok), data = rec.ph)
fit7 <- clogit(update(base, ~. + lab + comp1 + Qge05), data = rec.ph)
fit8 <- clogit(update(base, ~. + Wcrf_C_Cal + Qge05), data = rec.ph)
fit9 <- clogit(update(base, ~. + lab + comp1), data = recT.ph)
fit10 <- clogit(update(base, ~. + Wcrf_C_Cal), data = recT.ph)
fit11 <- clogit(update(base, ~. + Wcrf_C_Cal + comp1), data = rec.ph)

sigmodsD <- map_df(list(fit1, fit3, fit5, fit7, fit11, fit9), ~tidy(., exponentiate = T)) %>% filter(term == "comp1") %>% 
  mutate_if(is.numeric, ~round(., 2)) %>% unite(OR.CI, estimate, conf.low, conf.high, sep = "-")
scomodsD <- map_df(list(fit2, fit4, fit6, fit8, fit10), ~tidy(., exponentiate = T)) %>% filter(term == "Wcrf_C_Cal") %>% 
  mutate_if(is.numeric, ~round(., 2)) %>% unite(OR.CI, estimate, conf.low, conf.high, sep = "-")


# Make final table to copy and paste into manuscript
sigmods <- rbind(sigmodsC, sigmodsA, sigmodsB) %>% mutate_if(is.numeric, ~round(., 2)) %>% 
  unite(OR.CI, estimate, conf.low, conf.high, sep = "-")

scoremods <- rbind(scomodsC, scomodsA, scomodsB) %>% mutate_if(is.numeric, ~round(., 2)) %>%
  unite(OR.CI, estimate, conf.low, conf.high, sep = "-")


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
df2 <- crc %>% #crc[, cols] %>%
  select(Glyceroph_Lysopc_A_C17_0, Glyceroph_Pc_Ae_C40_6, Glyceroph_Pc_Ae_C36_2, Glyceroph_Pc_Ae_C38_2, Aminoacid_Ser, 
         Glyceroph_Lysopc_A_C18_2, Aminoacid_Gly, Glyceroph_Pc_Ae_C40_3, 
         # low:
         Glyceroph_Pc_Aa_C32_1, Glyceroph_Pc_Aa_C38_4, Glyceroph_Pc_Aa_C36_4, Aminoacid_Glu, Glyceroph_Pc_Aa_C34_4, 
         Glyceroph_Pc_Aa_C40_4, Glyceroph_Pc_Ae_C38_3) %>%
  mutate_all(~cut_number(., n = 4, labels = 1:4))

# Define function to apply across quartiles (already matched by lab)
clrfun <- function(x)  { clogit(Cncr_Caco_Clrt ~ x + Bmi_C + Qe_Energy + L_School + Smoke_Stat + Smoke_Int + Height_C +
                                  Qe_Alc + Qge0701 + Qge0704 + strata(Match_Caseset), data = crc) }
fits <- apply(df2, 2, clrfun)
mods1 <- map_df(fits, ~tidy(., exponentiate = T)) %>% filter(grepl("x4", term)) %>% mutate_if(is.numeric, ~round(., 2)) %>%
  unite(OR.CI, estimate, conf.low, conf.high, sep = "-")

# Fatty acids. Modelled continuously because few values.
df5 <- crc3[, colnames(ctrlC)] %>% 
  select(P17_0, P15_0, P15_1, P22_5n_6, P18_1n_9c, P16_1n_7c_n_9c, P16_0, P20_3n_9, P22_1n_9) %>% scale()

clrfun1 <- function(x)  { 
  clogit(Cncr_Caco_Clrt ~ x + Bmi_C + Qe_Energy + L_School + Smoke_Stat + Smoke_Int + Height_C + Qe_Alc +
                                  Qge0701 + Qge0704 + strata(Match_Caseset), data = crc3) 
  }
fits <- apply(df5, 2, clrfun1)
mods <- map_df(fits, ~tidy(., exponentiate = T)) %>% filter(grepl("x", term)) %>% mutate_if(is.numeric, ~round(., 2)) %>%
  unite(OR.CI, estimate, conf.low, conf.high, sep = "-")
