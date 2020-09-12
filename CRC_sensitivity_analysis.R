load("pred_score_tables_rev.Rdata")
library(survival)
library(broom)
base <- Cncr_Caco_Clrt ~ Qe_Energy + 
  L_School + Smoke_Stat + Height_C + strata(Match_Caseset)

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

sigmodsD <- map_df(list(fit7, fit7si, fit7sd, fit7dp, fit7t), ~tidy(., exponentiate = T)) %>% filter(term == "comp1")
scomodsD <- map_df(list(fit8, fit8si, fit8sd, fit8dp, fit8t), ~tidy(., exponentiate = T)) %>% filter(term == "Wcrf_C_Cal")


# Make final table to copy and paste into manuscript
sigmods <- rbind(sigmodsC, sigmodsA, sigmodsB) %>% mutate_if(is.numeric, ~round(., 2)) %>% 
  unite(OR.CI, estimate, conf.low, conf.high, sep = "-")

scoremods <- rbind(scomodsC, scomodsA, scomodsB) %>% mutate_if(is.numeric, ~round(., 2)) %>%
  unite(OR.CI, estimate, conf.low, conf.high, sep = "-")

