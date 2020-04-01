# Remove unneeded variables from workspace
rm(list = ls()[!str_detect(ls(), "pred.")])
rm(list = ls()[str_detect(ls(), "vars")])
save.image("predicted_score_tables.Rdata")

# Run CRC_signature_models

# Heterogeneity test for sex. Biocrates A and B and Fatty acids A
# Matching factors were age, sex, study centre, follow-up time since blood collection, fasting 
# status, menopausal status and phase of menstrual cycle at blood collection.

library(tidyverse)

# With matching factors (could not find phase menstrual cycle)
base <- Cncr_Caco_Clrt ~ Age_Blood + Tfollowup + #Phase_Mnscycle + #Menopause + 
  Center + Fasting_C + Qe_Energy + L_School + Smoke_Int + Height_C 

# Biocrates small, large, fatty acids (large fasted subset did not improve OR)

# CRC A: score and signature
fit1h <- glm(update(base, ~. + score.1.comps + Sex), data = pred.crc1, family = "binomial")
fit1i <- glm(update(base, ~. + score.1.comps * Sex), data = pred.crc1, family = "binomial")

library(broom)
library(lmtest)
s2 <- map_df(list(fit1h, fit1i), ~ tidy(., conf.int = T)) %>% filter(str_detect(term, "score."))
forest(s2$estimate, ci.lb = s2$conf.low, ci.ub = s2$conf.high, refline = 1, transf = exp) 
lrtest(fit1h, fit1i)
# pHET = 0.02

fit2h <- glm(update(base, ~. + Wcrf_C_Cal + Sex), data = pred.crc1, family = "binomial")
fit2i <- glm(update(base, ~. + Wcrf_C_Cal * Sex), data = pred.crc1, family = "binomial")

s2 <- map_df(list(fit2h, fit2i), ~ tidy(., conf.int = T)) %>% filter(str_detect(term, "Wcrf_"))
forest(s2$estimate, ci.lb = s2$conf.low, ci.ub = s2$conf.high, refline = 1, transf = exp) 
lrtest(fit2h, fit2i)
# p = 0.60



# CRC B: score and signature
fit3h <- glm(update(base, ~. + score.1.comps + Sex), data = pred.crc2, family = "binomial")
fit3i <- glm(update(base, ~. + score.1.comps * Sex), data = pred.crc2, family = "binomial")

s2 <- map_df(list(fit3h, fit3i), ~ tidy(., conf.int = T)) %>% filter(str_detect(term, "score."))
forest(s2$estimate, ci.lb = s2$conf.low, ci.ub = s2$conf.high, refline = 1, transf = exp) 
lrtest(fit3h, fit3i)
# p = 0.01

fit4h <- glm(update(base, ~. + Wcrf_C_Cal + Sex), data = pred.crc2, family = "binomial")
fit4i <- glm(update(base, ~. + Wcrf_C_Cal * Sex), data = pred.crc2, family = "binomial")

s2 <- map_df(list(fit4h, fit4i), ~ tidy(., conf.int = T)) %>% filter(str_detect(term, "Wcrf_"))
forest(s2$estimate, ci.lb = s2$conf.low, ci.ub = s2$conf.high, refline = 1, transf = exp) 
lrtest(fit4h, fit4i)
# p = 0.005


# CRC A fatty acids: score and signature
fit5h <- glm(update(base, ~. + score.1.comps + Sex), data = pred.fa, family = "binomial")
fit5i <- glm(update(base, ~. + score.1.comps * Sex), data = pred.fa, family = "binomial")

s2 <- map_df(list(fit5h, fit5i), ~ tidy(., conf.int = T)) %>% filter(str_detect(term, "score."))
forest(s2$estimate, ci.lb = s2$conf.low, ci.ub = s2$conf.high, refline = 1, transf = exp) 
lrtest(fit5h, fit5i)
# pHET = 0.26

fit6h <- glm(update(base, ~. + Wcrf_C_Cal + Sex), data = pred.fa, family = "binomial")
fit6i <- glm(update(base, ~. + Wcrf_C_Cal * Sex), data = pred.fa, family = "binomial")

s2 <- map_df(list(fit6h, fit6i), ~ tidy(., conf.int = T)) %>% filter(str_detect(term, "Wcrf_"))
forest(s2$estimate, ci.lb = s2$conf.low, ci.ub = s2$conf.high, refline = 1, transf = exp) 
lrtest(fit6h, fit6i)
# p = 0.52

