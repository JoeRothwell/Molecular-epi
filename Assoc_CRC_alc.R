#load("pred_score_tables_rev.Rdata")
load("alc_mediation.Rdata")
library(haven)

# Alcohol-CRC associations
library(survival)
library(sjmisc)
base <- Cncr_Caco_Clrt ~ Bmi_C + Qe_Energy + L_School + Smoke_Stat + Smoke_Int + Height_C + 
  Qge0701 + strata(Match_Caseset)

# Make dichotomised "gate" variables
crc.ph$Qe_Alc_cat <- dicho(crc.ph$Qe_Alc, dich.by = 0.5)
crc3.ph$Qe_Alc_cat <- dicho(crc3.ph$Qe_Alc, dich.by = 0.5)
crc1$Qe_Alc_cat <- dicho(crc1$Qe_Alc, dich.by = 0.5)

# Both case-controls without and with gate variable
fit1 <- clogit(update(base, ~. + I(Qe_Alc/12)), data = crc.ph) #1.06 (1.00-1.11)
fit2 <- clogit(update(base, ~. + I(Qe_Alc/12) + Qe_Alc_cat), data = crc.ph) #1.07 (1.01-1.13)


# Small case-control without and with gate variable
fit3 <- clogit(update(base, ~. + I(Qe_Alc/12)), data = crc1) #1.14 (1.00-1.28)
fit4 <- clogit(update(base, ~. + I(Qe_Alc/12) + Qe_Alc_cat), data = crc1) #1.18 (1.03-1.34)


# Fatty acids case-control without and with gate variable
fit5 <- clogit(update(base, ~. + I(Qe_Alc/12)), data = crc3.ph) #1.15 (1.00-1.29)
fit6 <- clogit(update(base, ~. + I(Qe_Alc/12) + Qe_Alc_cat), data = crc3.ph) #1.19 (1.04-1.35)

modlist <- list(fit1, fit2, fit3, fit4, fit5, fit6)

library(broom)
scomods <- map_df(modlist, ~tidy(., exponentiate = T)) %>% filter(str_detect(term, "Alc/12")) %>%
  mutate_if(is.numeric, ~round(., 2)) %>% unite(OR.CI, estimate, conf.low, conf.high, sep = "-")

# Scaled and categorical variables (no associations)
#fit2 <- clogit(update(base, ~. + lab + scale(Qe_Alc)), data = crc.ph)
#fit3 <- clogit(update(base, ~. + lab + cut_number(Qe_Alc, n = 4)), data = crc.ph)

# Prepare data for PLS. Subset controls only
expr <- "(carn|oacid|genic|roph|ingo|Sugars)[_]"
ctlmat <- crc1 %>% filter(Cncr_Caco_Clrt == 0) %>%
  select(matches(expr), -contains("tdq")) %>% select_if(~ sum(., na.rm = T) != 0) %>% 
  log2 %>% scale

# Control subset all data
dat <- crc1[crc1$Cncr_Caco_Clrt == 0, ]

library(lme4)
adj <- function(x) residuals(lm(x ~ Center + Batch_MetBio, data = dat))
adjmat <- apply(ctlmat, 2, adj)

library(pls)
# Bind log transformed alcohol intakes to adjusted metabolite matrix
#plsdat <- bind_cols(alc = log2(dat$Qe_Alc + 0.1), adjmat) #%>% filter(!is.na(score))
plsdat <- bind_cols(alc = log(dat$Qe_Alc/12 + 0.01), adjmat)

set.seed(111)
mod <- plsr(alc ~ ., ncomp = 10, data = plsdat, validation = "CV")
# Find the number of dimensions with lowest cross validation error
cv <- RMSEP(mod)
plot(RMSEP(mod), legendpos = "topright")

# Get components with lowest RMSEP, min, "one SE" and "permutation" methods
which.min(cv$val[estimate = "adjCV", , ]) - 1
selectNcomp(mod, method = "onesigma", plot = F)
selectNcomp(mod, method = "randomization", plot = T)

# Use 2 components
mod1 <- plsr(alc ~ ., data = plsdat, ncomp = 2)
# Get coefficients
coeff <- data.frame(value = round(coef(mod)[ , 1, 1], 3))

# Get cases for prediction
mat <- crc1 %>% #filter(Cncr_Caco_Clrt == 1) %>%
  select(matches(expr), -contains("tdq")) %>% select_if(~ sum(., na.rm = T) != 0) %>% 
  log %>% scale

crc1$pred.alc <- predict(mod1, mat)[,,1]
fit3a <- clogit(update(base, ~. + exp(pred.alc)), data = crc1) #2.63 (1.43-4.81)
summary(fit3a)

# NDE (direct effect, adjusting for mediator)
modY <- clogit(update(base, ~. + I(Qe_Alc/12) + exp(pred.alc)), data = crc1) #1.08 (0.95-1.22)
summary(fit3b)

# NIE (indirect effect, mediator modelled by exposure)
# Is the mediator predicted by the exposure
modM <- lm(exp(pred.alc) ~ Bmi_C + Qe_Energy + L_School + Smoke_Stat + Smoke_Int + Height_C + 
              Qge0701 + I(Qe_Alc/12), data = crc1)

library(mediation)
mediate(modM, modY, treat = "Qe_Alc", mediator = "pred.alc", boot = T, sims = 100)


