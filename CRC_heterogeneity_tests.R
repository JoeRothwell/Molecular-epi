# Remove unneeded variables from workspace
rm(list = ls()[!str_detect(ls(), ".ph|.both")])
rm(list = ls()[str_detect(ls(), "vars")])

# Load the predicted score tables for modelling
load("predicted_score_tables.Rdata")

# Heterogeneity test for sex. Biocrates A and B and Fatty acids A
# Matching factors were age, sex, study centre, follow-up time since blood collection, fasting 
# status, menopausal status and phase of menstrual cycle at blood collection.

library(tidyverse)
# Normal GLM with matching factors. Follow up time is meaningless for controls
# Menopause variables are incomplete
base <- Cncr_Caco_Clrt ~ Age_Blood + #Tfollowup + #Phase_Mnscycle + #Menopause + 
  Center + Fasting_C + Qe_Energy + L_School + Smoke_Stat + Smoke_Int + Height_C 

# Biocrates small, large, fatty acids (large fasted subset did not improve OR)

# CRC A: score and signature
library(lmtest)
fit1h <- glm(update(base, ~. + comp1 + Sex), data = crc1.ph, family = "binomial")
fit1i <- glm(update(base, ~. + comp1 * Sex), data = crc1.ph, family = "binomial")
lrtest(fit1h, fit1i)
# pHET = 0.021

library(broom)
s2 <- map_df(list(fit1h, fit1i), ~ tidy(., conf.int = T)) %>% filter(str_detect(term, "score."))
forest(s2$estimate, ci.lb = s2$conf.low, ci.ub = s2$conf.high, refline = 1, transf = exp) 


fit2h <- glm(update(base, ~. + Wcrf_C_Cal + Sex), data = crc1.ph, family = "binomial")
fit2i <- glm(update(base, ~. + Wcrf_C_Cal * Sex), data = crc1.ph, family = "binomial")
lrtest(fit2h, fit2i)
# pHET = 0.64

s2 <- map_df(list(fit2h, fit2i), ~ tidy(., conf.int = T)) %>% filter(str_detect(term, "Wcrf_"))
forest(s2$estimate, ci.lb = s2$conf.low, ci.ub = s2$conf.high, refline = 1, transf = exp) 



# CRC B: score and signature
fit3h <- glm(update(base, ~. + comp1 + Sex), data = crc2.ph, family = "binomial")
fit3i <- glm(update(base, ~. + comp1 * Sex), data = crc2.ph, family = "binomial")
lrtest(fit3h, fit3i)
# p = 0.13

s2 <- map_df(list(fit3h, fit3i), ~ tidy(., conf.int = T)) %>% filter(str_detect(term, "score."))
forest(s2$estimate, ci.lb = s2$conf.low, ci.ub = s2$conf.high, refline = 1, transf = exp) 


fit4h <- glm(update(base, ~. + Wcrf_C_Cal + Sex), data = crc2.ph, family = "binomial")
fit4i <- glm(update(base, ~. + Wcrf_C_Cal * Sex), data = crc2.ph, family = "binomial")
lrtest(fit4h, fit4i)
# p = 0.008

s2 <- map_df(list(fit4h, fit4i), ~ tidy(., conf.int = T)) %>% filter(str_detect(term, "Wcrf_"))
forest(s2$estimate, ci.lb = s2$conf.low, ci.ub = s2$conf.high, refline = 1, transf = exp) 


# CRC A fatty acids: score and signature
fit5h <- glm(update(base, ~. + comp2 + Sex), data = crc3.ph, family = "binomial")
fit5i <- glm(update(base, ~. + comp2 * Sex), data = crc3.ph, family = "binomial")
lrtest(fit5h, fit5i)
# pHET = 0.06

s2 <- map_df(list(fit5h, fit5i), ~ tidy(., conf.int = T)) %>% filter(str_detect(term, "score."))
forest(s2$estimate, ci.lb = s2$conf.low, ci.ub = s2$conf.high, refline = 1, transf = exp) 


fit6h <- glm(update(base, ~. + Wcrf_C_Cal + Sex), data = crc3.ph, family = "binomial")
fit6i <- glm(update(base, ~. + Wcrf_C_Cal * Sex), data = crc3.ph, family = "binomial")
lrtest(fit6h, fit6i)
# p = 0.35

s2 <- map_df(list(fit6h, fit6i), ~ tidy(., conf.int = T)) %>% filter(str_detect(term, "Wcrf_"))
forest(s2$estimate, ci.lb = s2$conf.low, ci.ub = s2$conf.high, refline = 1, transf = exp) 

