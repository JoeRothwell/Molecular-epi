# Remove unneeded variables from workspace
rm(list = ls()[!str_detect(ls(), ".ph|.both")])
rm(list = ls()[str_detect(ls(), "vars")])

# Load the predicted score tables for modelling
load("predicted_score_tables_sex.Rdata")

# Heterogeneity test for sex. Biocrates A and B and Fatty acids A
# Matching factors were age, sex, study centre, follow-up time since blood collection, fasting 
# status, menopausal status and phase of menstrual cycle at blood collection.

### Revised manuscript: studies A and B merged ###

library(tidyverse)
# Normal GLM with matching factors. Follow up time is meaningless for controls
# Menopause variables are incomplete
base <- Cncr_Caco_Clrt ~ Age_Blood + #Tfollowup + #Phase_Mnscycle + #Menopause + #Qe_Energy + 
  Center + Fasting_C + L_School + Smoke_Stat + Smoke_Int + Height_C 

# Biocrates small, large, fatty acids (large fasted subset did not improve OR)

# CRC endogenous: score and signature
library(lmtest)
fit1h <- glm(update(base, ~. + comp1 + Sex), data = crc.ph, family = "binomial")
fit1i <- glm(update(base, ~. + comp1 * Sex), data = crc.ph, family = "binomial")
lrtest(fit1h, fit1i)
# pHET = 0.022

library(broom)
s2 <- map_df(list(fit1h, fit1i), ~ tidy(., conf.int = T)) %>% filter(str_detect(term, "score."))
forest(s2$estimate, ci.lb = s2$conf.low, ci.ub = s2$conf.high, refline = 1, transf = exp) 


fit2h <- glm(update(base, ~. + Wcrf_C_Cal + Sex), data = crc.ph, family = "binomial")
fit2i <- glm(update(base, ~. + Wcrf_C_Cal * Sex), data = crc.ph, family = "binomial")
lrtest(fit2h, fit2i)
# pHET = 0.022

s2 <- map_df(list(fit2h, fit2i), ~ tidy(., conf.int = T)) %>% filter(str_detect(term, "Wcrf_"))
forest(s2$estimate, ci.lb = s2$conf.low, ci.ub = s2$conf.high, refline = 1, transf = exp) 


# CRC A fatty acids: score and signature
fit5h <- glm(update(base, ~. + comp2 + Sex), data = crc3.ph, family = "binomial")
fit5i <- glm(update(base, ~. + comp2 * Sex), data = crc3.ph, family = "binomial")
lrtest(fit5h, fit5i)
# pHET = 0.054

s2 <- map_df(list(fit5h, fit5i), ~ tidy(., conf.int = T)) %>% filter(str_detect(term, "score."))
forest(s2$estimate, ci.lb = s2$conf.low, ci.ub = s2$conf.high, refline = 1, transf = exp) 


fit6h <- glm(update(base, ~. + Wcrf_C_Cal + Sex), data = crc3.ph, family = "binomial")
fit6i <- glm(update(base, ~. + Wcrf_C_Cal * Sex), data = crc3.ph, family = "binomial")
lrtest(fit6h, fit6i)
# p = 0.35

s2 <- map_df(list(fit6h, fit6i), ~ tidy(., conf.int = T)) %>% filter(str_detect(term, "Wcrf_"))
forest(s2$estimate, ci.lb = s2$conf.low, ci.ub = s2$conf.high, refline = 1, transf = exp) 

# Colon endogenous: score and signature
fit7h <- glm(update(base, ~. + comp1 + Sex), data = col.ph, family = "binomial")
fit7i <- glm(update(base, ~. + comp1 * Sex), data = col.ph, family = "binomial")
lrtest(fit7h, fit7i)
# pHET = 0.025

fit8h <- glm(update(base, ~. + Wcrf_C_Cal + Sex), data = col.ph, family = "binomial")
fit8i <- glm(update(base, ~. + Wcrf_C_Cal * Sex), data = col.ph, family = "binomial")
lrtest(fit8h, fit8i)
# pHET = 0.002

# Rectal endogenous: score and signature
fit3h <- glm(update(base, ~. + comp1 + Sex), data = rec.ph, family = "binomial")
fit3i <- glm(update(base, ~. + comp1 * Sex), data = rec.ph, family = "binomial")
lrtest(fit3h, fit3i)
# pHET = 0.49

fit4h <- glm(update(base, ~. + Wcrf_C_Cal + Sex), data = rec.ph, family = "binomial")
fit4i <- glm(update(base, ~. + Wcrf_C_Cal * Sex), data = rec.ph, family = "binomial")
lrtest(fit4h, fit4i)
# pHET = 0.8346

