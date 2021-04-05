source("Make_alc_mediators.R")
library(haven)
library(survival)
library(sjmisc)

# Alcohol-CRC associations
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

# Summarise OR (CI)
library(broom)
scomods <- map_df(modlist, ~tidy(., exponentiate = T)) %>% filter(str_detect(term, "Alc/12")) %>%
  mutate_if(is.numeric, ~round(., 2)) %>% unite(OR.CI, estimate, conf.low, conf.high, sep = "-")

# Scaled and categorical variables (no associations)
#fit2 <- clogit(update(base, ~. + lab + scale(Qe_Alc)), data = crc.ph)
#fit3 <- clogit(update(base, ~. + lab + cut_number(Qe_Alc, n = 4)), data = crc.ph)

### Calculation of NDE, NIE and RD ratio from Van Steenland 2010 (pg 1342)
# Run Make_alc_mediators.R

# Mediator 1 - controls from small case-control (as Laura did)
# Logistic model adjusting for mediator. Need theta1 (exposure coef) and theta2 
# (mediator coef) from this
modY <- clogit(update(base, ~. + I(Qe_Alc/12) + M.alc), data = crc1) 

# Linear model of outcome and mediator (OR from beta coefficients). Need beta1 (exposure coef)
modM <- lm(M.alc ~ Bmi_C + Qe_Energy + L_School + Smoke_Stat + Smoke_Int + Height_C + 
              Qge0701 + I(Qe_Alc/12), data = crc1)

# Natural direct effect is given as exp(coeff for 1 unit change in exposure) theta1
theta1s <- tidy(modY, exponentiate = T, conf.int = T)[20, -1] 

#NDE = 1.08 (0.95-1.22)
#      1.07 (0.94-1.22) for new mediator derived without zero alc values

# Natural indirect effect is given as exp(coeff of theta for mediator x beta for exposure)
theta2s <- tidy(modY, conf.int = T)[21, -1]
theta2s <- tidy(modY, conf.int = T)[21, ]

# For betas need to get CI separately
beta1s <- tidy(modM, conf.int = T)[21, -1] %>% as.numeric
ci     <- confint(modM, "I(Qe_Alc/12)")
beta1s <- c(beta1s, ci)
exp(theta2s*beta1s)
# NIE = 1.08 (1.02-1.17)
#       1.08 (1.02-1.18) for new mediator derived without zero alc values

# Summary
# TE = 1.14 (1.00-1.28)
# NDE = 1.07 (0.94-1.22)
# NIE = 1.08 (1.02-1.18)

# Risk difference ratio is defined as log(TE) - log(NDE) / log (NDE)  *100
logTE <- coef(fit3)[20] 
logNDE <- coef(modY)[20]
RDR <- 100 * (logTE - logNDE) / logNDE # 74.0%; 86.9% with new mediator


### Mediator 1 with gate variable
# Natural direct effect is given as exp(coeff for 1 unit change in exposure) theta1
modY <- clogit(update(base, ~. + I(Qe_Alc/12) + Qe_Alc_cat + M.alc), data = crc1) 
theta1s <- tidy(modY, exponentiate = T, conf.int = T)[20, -1] #NDE = 1.11 (0.97-1.28)
theta2s <- tidy(modY, conf.int = T)[21, -1]

# Linear model of outcome and mediator (OR from beta coefficients). Need beta1 (exposure coef)
modM <- lm(M.alc ~ Bmi_C + Qe_Energy + L_School + Smoke_Stat + Smoke_Int + Height_C + 
             Qge0701 + I(Qe_Alc/12) + Qe_Alc_cat, data = crc1)
# Need to get CIs separately
beta1s <- tidy(modM, conf.int = T)[21, -1] %>% as.numeric
ci     <- confint(modM, "I(Qe_Alc/12)")
beta1s <- c(beta1s, ci)
exp(theta2s*beta1s)
# NIE = 0.98 (0.96-1.01)

# Summary
# TE = 1.18 (1.03-1.34)
# NDE = 1.11 (0.97-1.28)

logTE <- coef(fit4)[20] 
logNDE <- coef(modY)[20]
RDR <- 100 * (logTE - logNDE) / logNDE # 51.2%


### Mediator 2 - EPIC pooled controls
modY <- clogit(update(base, ~. + I(Qe_Alc/12) + M.alc1), data = crc1) 

# Linear model of outcome and mediator (OR from beta coefficients). Need beta1 (exposure coef)
modM <- lm(M.alc1 ~ Bmi_C + Qe_Energy + L_School + Smoke_Stat + Smoke_Int + Height_C + 
             Qge0701 + I(Qe_Alc/12), data = crc1)

# Natural direct effect is given as exp(coeff for 1 unit change in exposure) theta1
theta1s <- tidy(modY, exponentiate = T)[20, -1] #NDE = 1.16 (1.02-1.31)

# Natural indirect effect is given as exp(coeff of theta for mediator x beta for exposure)
theta2s <- tidy(modY, conf.int = T)[21, -1]

# For betas need to get CI separately
beta1s <- tidy(modM)[21, -1] %>% as.numeric
ci     <- confint(modM, "I(Qe_Alc/12)")
beta1s <- c(beta1s, ci)
exp(theta2s*beta1s)
# NIE = 0.98 (0.97-1.01)

# Summary
# TE = 1.14 (1.00-1.28)
# NDE = 1.16 (1.02-1.31)

# Risk difference ratio is defined as log(TE) - log(NDE) / log (NDE)  *100
logTE <- coef(fit3)[20] 
logNDE <- coef(modY)[20]
RDR <- 100 * (logTE - logNDE) / logNDE # -12.9% (-13.4 for new mediator with zeros excluded)

# With gate variable
modY <- clogit(update(base, ~. + I(Qe_Alc/12) + Qe_Alc_cat + M.alc1), data = crc1) 

# Linear model of outcome and mediator (OR from beta coefficients). Need beta1 (exposure coef)
modM <- lm(M.alc1 ~ Bmi_C + Qe_Energy + L_School + Smoke_Stat + Smoke_Int + Height_C + 
             Qge0701 + I(Qe_Alc/12) + Qe_Alc_cat, data = crc1)

# Natural direct effect is given as exp(coeff for 1 unit change in exposure) theta1
theta1s <- tidy(modY, conf.int = T)[20, -1] #NDE = 1.16 (1.02-1.31)

# Natural indirect effect is given as exp(coeff of theta for mediator x beta for exposure)
theta2s <- tidy(modY, conf.int = T)[21, -1]

# For betas need to get CI separately
beta1s <- tidy(modM)[21, -1] %>% as.numeric
ci     <- confint(modM, "I(Qe_Alc/12)")
beta1s <- c(beta1s, ci)
exp(theta2s*beta1s)
# NIE = 0.99 (0.98-1.00)

# Summary
# TE = 1.18 (1.00-1.34)
# NDE = 1.16 (1.02-1.31)
# NIE = 0.99 (0.98-1.00)

# Risk difference ratio is defined as log(TE) - log(NDE) / log (NDE)  *100
logTE <- coef(fit4)[20] 
logNDE <- coef(modY)[20]
RDR <- 100 * (logTE - logNDE) / logNDE # -10.8%
