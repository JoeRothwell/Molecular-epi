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
fit1h <- glm(update(base, ~. + lab + comp1 + Sex), data = crc.ph, family = "binomial")
fit1i <- glm(update(base, ~. + lab + comp1 * Sex), data = crc.ph, family = "binomial")
lrtest(fit1h, fit1i)
# pHET = 0.029

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
# pHET = 0.072

s2 <- map_df(list(fit5h, fit5i), ~ tidy(., conf.int = T)) %>% filter(str_detect(term, "score."))
forest(s2$estimate, ci.lb = s2$conf.low, ci.ub = s2$conf.high, refline = 1, transf = exp) 


fit6h <- glm(update(base, ~. + Wcrf_C_Cal + Sex), data = crc3.ph, family = "binomial")
fit6i <- glm(update(base, ~. + Wcrf_C_Cal * Sex), data = crc3.ph, family = "binomial")
lrtest(fit6h, fit6i)
# p = 0.36

s2 <- map_df(list(fit6h, fit6i), ~ tidy(., conf.int = T)) %>% filter(str_detect(term, "Wcrf_"))
forest(s2$estimate, ci.lb = s2$conf.low, ci.ub = s2$conf.high, refline = 1, transf = exp) 

# Colon endogenous: score and signature
fit7h <- glm(update(base, ~. + lab + comp1 + Sex), data = col.ph, family = "binomial")
fit7i <- glm(update(base, ~. + lab + comp1 * Sex), data = col.ph, family = "binomial")
lrtest(fit7h, fit7i)
# pHET = 0.03

fit8h <- glm(update(base, ~. + Wcrf_C_Cal + Sex), data = col.ph, family = "binomial")
fit8i <- glm(update(base, ~. + Wcrf_C_Cal * Sex), data = col.ph, family = "binomial")
lrtest(fit8h, fit8i)
# pHET = 0.002

# Prox colon endogenous: score and signature
fit7h <- glm(update(base, ~. + lab + comp1 + Sex), data = prox.ph, family = "binomial")
fit7i <- glm(update(base, ~. + lab + comp1 * Sex), data = prox.ph, family = "binomial")
lrtest(fit7h, fit7i)
# pHET = 0.21

fit8h <- glm(update(base, ~. + Wcrf_C_Cal + Sex), data = prox.ph, family = "binomial")
fit8i <- glm(update(base, ~. + Wcrf_C_Cal * Sex), data = prox.ph, family = "binomial")
lrtest(fit8h, fit8i)
# pHET = 0.12

# Dist endogenous: score and signature
fit7h <- glm(update(base, ~. + lab + comp1 + Sex), data = dist.ph, family = "binomial")
fit7i <- glm(update(base, ~. + lab + comp1 * Sex), data = dist.ph, family = "binomial")
lrtest(fit7h, fit7i)
# pHET = 0.12

fit8h <- glm(update(base, ~. + Wcrf_C_Cal + Sex), data = dist.ph, family = "binomial")
fit8i <- glm(update(base, ~. + Wcrf_C_Cal * Sex), data = dist.ph, family = "binomial")
lrtest(fit8h, fit8i)
# pHET = 0.005

# Rectal endogenous: score and signature
fit3h <- glm(update(base, ~. + lab + comp1 + Sex), data = rec.ph, family = "binomial")
fit3i <- glm(update(base, ~. + lab + comp1 * Sex), data = rec.ph, family = "binomial")
lrtest(fit3h, fit3i)
# pHET = 0.46

fit4h <- glm(update(base, ~. + Wcrf_C_Cal + Sex), data = rec.ph, family = "binomial")
fit4i <- glm(update(base, ~. + Wcrf_C_Cal * Sex), data = rec.ph, family = "binomial")
lrtest(fit4h, fit4i)
# pHET = 0.8346

# Colon fatty acids: score and signature
fit9h <- glm(update(base, ~. + Wcrf_C_Cal + Sex), data = col3.ph, family = "binomial")
fit9i <- glm(update(base, ~. + Wcrf_C_Cal * Sex), data = col3.ph, family = "binomial")
lrtest(fit9h, fit9i)
# pHET = 0.2836

fit10h <- glm(update(base, ~. + comp2 + Sex), data = col3.ph, family = "binomial")
fit10i <- glm(update(base, ~. + comp2 * Sex), data = col3.ph, family = "binomial")
lrtest(fit10h, fit10i)
# pHET = 0.1148

# Prox fatty acids: score and signature
fit11h <- glm(update(base, ~. + Wcrf_C_Cal + Sex), data = prox3.ph, family = "binomial")
fit11i <- glm(update(base, ~. + Wcrf_C_Cal * Sex), data = prox3.ph, family = "binomial")
lrtest(fit11h, fit11i)
# pHET = 0.4379

fit12h <- glm(update(base, ~. + comp2 + Sex), data = prox3.ph, family = "binomial")
fit12i <- glm(update(base, ~. + comp2 * Sex), data = prox3.ph, family = "binomial")
lrtest(fit12h, fit12i)
# pHET = 0.4328

# Dist fatty acids: score and signature
fit13h <- glm(update(base, ~. + Wcrf_C_Cal + Sex), data = dist3.ph, family = "binomial")
fit13i <- glm(update(base, ~. + Wcrf_C_Cal * Sex), data = dist3.ph, family = "binomial")
lrtest(fit13h, fit13i)
# pHET = 0.4864

fit14h <- glm(update(base, ~. + comp2 + Sex), data = dist3.ph, family = "binomial")
fit14i <- glm(update(base, ~. + comp2 * Sex), data = dist3.ph, family = "binomial")
lrtest(fit14h, fit14i)
# pHET = 0.1846
