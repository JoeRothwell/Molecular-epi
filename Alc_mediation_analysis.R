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
#fit1 <- clogit(update(base, ~. + I(Qe_Alc/12)), data = crc.ph) #1.06 (1.00-1.11)
#fit2 <- clogit(update(base, ~. + I(Qe_Alc/12) + Qe_Alc_cat), data = crc.ph) #1.07 (1.01-1.13)

# Small case-control without and with gate variable
fit3 <- clogit(update(base, ~. + I(Qe_Alc/12)), data = crc1) #1.14 (1.00-1.28)
fit4 <- clogit(update(base, ~. + I(Qe_Alc/12) + Qe_Alc_cat), data = crc1) #1.18 (1.03-1.34)

# Fatty acids case-control without and with gate variable
fit5 <- clogit(update(base, ~. + I(Qe_Alc/12)), data = crc3.ph) #1.15 (1.00-1.29)
fit6 <- clogit(update(base, ~. + I(Qe_Alc/12) + Qe_Alc_cat), data = crc3.ph) #1.19 (1.04-1.35)

#modlist <- list(fit1, fit2, fit3, fit4, fit5, fit6)

# Summarise OR (CI)
library(broom)
#scomods <- map_df(modlist, ~tidy(., exponentiate = T)) %>% filter(str_detect(term, "Alc/12")) %>%
#  mutate_if(is.numeric, ~round(., 2)) %>% unite(OR.CI, estimate, conf.low, conf.high, sep = "-")

# Scaled and categorical variables (no associations)
#fit2 <- clogit(update(base, ~. + lab + scale(Qe_Alc)), data = crc.ph)
#fit3 <- clogit(update(base, ~. + lab + cut_number(Qe_Alc, n = 4)), data = crc.ph)

### Calculation of NDE, NIE and RD ratio from Van Steenland 2010 (pg 1342)
# Run Make_alc_mediators.R

# M1: controls from small case-control, all biocrates compounds (as Laura did)
# modY is a CLR adjusting for the mediator. Need exposure and mediator coefs theta1, theta2 
modY <- clogit(update(base, ~. + I(Qe_Alc/12) + M1), data = crc1) 

# modM is a linear model of outcome and mediator. Need exposure coef beta1
modM <- lm(M1 ~ Bmi_C + Qe_Energy + L_School + Smoke_Stat + Smoke_Int + Height_C + 
              Qge0701 + I(Qe_Alc/12), data = crc1)

# Natural direct effect is given as exp(coeff for 1 unit change in exposure), theta1
theta1 <- tidy(modY, exponentiate = T, conf.int = T) %>% filter(term == "I(Qe_Alc/12)") %>% select(-1)

# Natural indirect effect is exp(coeff of theta for mediator x beta for exposure)
# Get Theta coefficients of mod Y and beta coefficients of modM:
theta2 <- tidy(modY, conf.int = T) %>% filter(term == "M1") %>% select(-1)
beta1  <- tidy(modM, conf.int = T) %>% filter(term == "I(Qe_Alc/12)") %>% select(-1)
NIE <- exp(theta2*beta1)

# TE = 1.14 (1.00-1.28), NDE = 1.07 (0.94-1.22), NIE = 1.10 (1.02-1.23)

# Risk difference ratio is defined as log(TE) - log(NDE) / log (NDE)  *100
logTE <- coef(fit3)[20] 
logNDE <- coef(modY)[20]
RDR <- 100 * (logTE - logNDE) / logNDE # 90.6

# With gate variable
# Natural direct effect is given as exp(coeff for 1 unit change in exposure) theta1
modY <- clogit(update(base, ~. + I(Qe_Alc/12) + Qe_Alc_cat + M1), data = crc1) 
theta1 <- tidy(modY, exponentiate = T, conf.int = T) %>% filter(term == "I(Qe_Alc/12)") %>% select(-1)
theta2 <- tidy(modY, conf.int = T) %>% filter(term == "M1") %>% select(-1)

# Linear model of outcome and mediator (OR from beta coefficients). Need beta1 (exposure coef)
modM <- lm(M1 ~ Bmi_C + Qe_Energy + L_School + Smoke_Stat + Smoke_Int + Height_C + 
             Qge0701 + I(Qe_Alc/12) + Qe_Alc_cat, data = crc1)

beta1 <- tidy(modM, conf.int = T) %>% filter(term == "I(Qe_Alc/12)") %>% select(-1)
NIE <- exp(theta2*beta1)

# TE = 1.18 (1.03-1.34), NDE = 1.11 (0.96-1.27), NIE = 1.10 (1.01-1.24)

logTE <- coef(fit4)[20] 
logNDE <- coef(modY)[20]
RDR <- 100 * (logTE - logNDE) / logNDE # 61.6%


# M2: controls from small study, amino acids only
modY <- clogit(update(base, ~. + I(Qe_Alc/12) + M2), data = crc1) 
modM <- lm(M2 ~ Bmi_C + Qe_Energy + L_School + Smoke_Stat + Smoke_Int + Height_C + 
             Qge0701 + I(Qe_Alc/12), data = crc1)

# Natural direct effect is given as exp(coeff for 1 unit change in exposure) theta1
theta1 <- tidy(modY, exponentiate = T, conf.int = T) %>% filter(term == "M2") %>% select(-1)
theta2 <- tidy(modY, conf.int = T) %>% filter(term == "M2") %>% select(-1)
beta1  <- tidy(modM, conf.int = T) %>% filter(term == "I(Qe_Alc/12)") %>% select(-1)
NIE <- exp(theta2 * beta1)
# TE = 1.14 (1.00-1.28), NDE = 1.14 (1.01-1.29), NIE = 0.99 (1.00-1.00)
logTE <- coef(fit3)[20] 
logNDE <- coef(modY)[20]
RDR <- 100 * (logTE - logNDE) / logNDE # -3.9

# M2 + gate
modY <- clogit(update(base, ~. + I(Qe_Alc/12) + Qe_Alc_cat + M2), data = crc1) 
modM <- lm(M2 ~ Bmi_C + Qe_Energy + L_School + Smoke_Stat + Smoke_Int + Height_C + 
             Qge0701 + I(Qe_Alc/12) + Qe_Alc_cat, data = crc1)

# Natural direct effect is given as exp(coeff for 1 unit change in exposure) theta1
theta1 <- tidy(modY, exponentiate = T, conf.int = T) %>% filter(term == "I(Qe_Alc/12)") %>% select(-1)
theta2 <- tidy(modY, conf.int = T) %>% filter(term == "M2") %>% select(-1)
beta1  <- tidy(modM, conf.int = T) %>% filter(term == "I(Qe_Alc/12)") %>% select(-1)
NIE <- exp(theta2 * beta1)
# TE = 1.18 (1.03-1.34), NDE = 1.19 (1.04-1.35), NIE = 0.99 (1.00-1.00)
logTE <- coef(fit4)[20] 
logNDE <- coef(modY)[20]
RDR <- 100 * (logTE - logNDE) / logNDE # -3.9

# M3: controls from small study, lipid metabolites only
modY <- clogit(update(base, ~. + I(Qe_Alc/12) + M3), data = crc1) 
modM <- lm(M3 ~ Bmi_C + Qe_Energy + L_School + Smoke_Stat + Smoke_Int + Height_C + 
             Qge0701 + I(Qe_Alc/12), data = crc1)

# Natural direct effect is given as exp(coeff for 1 unit change in exposure) theta1
theta1s <- tidy(modY, exponentiate = T, conf.int = T) %>% filter(term == "I(Qe_Alc/12)") %>% select(-1)
theta2s <- tidy(modY, conf.int = T) %>% filter(term == "M3") %>% select(-1)
beta1s  <- tidy(modM, conf.int = T) %>% filter(term == "I(Qe_Alc/12)") %>% select(-1)
NIE <- exp(theta2s*beta1s)
# TE = 1.14 (1.00-1.28), NDE = 1.06 (0.93-1.21), NIE = 1.10 (1.02-1.21)
logTE <- coef(fit3)[20] 
logNDE <- coef(modY)[20]
RDR <- 100 * (logTE - logNDE) / logNDE # 118.6

# M3 + gate
modY <- clogit(update(base, ~. + I(Qe_Alc/12) + Qe_Alc_cat + M3), data = crc1) 
modM <- lm(M3 ~ Bmi_C + Qe_Energy + L_School + Smoke_Stat + Smoke_Int + Height_C + 
             Qge0701 + I(Qe_Alc/12) + Qe_Alc_cat, data = crc1)

# Natural direct effect is given as exp(coeff for 1 unit change in exposure) theta1
theta1 <- tidy(modY, exponentiate = T, conf.int = T) %>% filter(term == "M3") %>% select(-1)
theta2 <- tidy(modY, conf.int = T) %>% filter(term == "M3") %>% select(-1)
beta1  <- tidy(modM, conf.int = T) %>% filter(term == "I(Qe_Alc/12)") %>% select(-1)
NIE <- exp(theta2*beta1)
# TE = 1.18 (1.03-1.34), NDE = 1.06 (0.93-1.21), NIE = 1.10 (1.02-1.22)
logTE <- coef(fit4)[20] 
logNDE <- coef(modY)[20]
RDR <- 100 * (logTE - logNDE) / logNDE # 78.0



### Mediator derived from EPIC non-CRC pooled controls
modY <- clogit(update(base, ~. + I(Qe_Alc/12) + M4), data = crc1) 
modM <- lm(M4 ~ Bmi_C + Qe_Energy + L_School + Smoke_Stat + Smoke_Int + Height_C + 
             Qge0701 + I(Qe_Alc/12), data = crc1)

# Natural direct effect is given as exp(coeff for 1 unit change in exposure) theta1
theta1 <- tidy(modY, exponentiate = T, conf.int = T) %>% filter(term == "I(Qe_Alc/12)") %>% select(-1)
theta2 <- tidy(modY, conf.int = T) %>% filter(term == "M4") %>% select(-1)
beta1 <- tidy(modM, conf.int = T) %>% filter(term == "I(Qe_Alc/12)") %>% select(-1)
NIE <- exp(theta2*beta1)

#TE = 1.14 (1.00-1.28), NDE = 1.16 (1.02-1.31), NIE = 0.99 (0.98-1.01)
logTE <- coef(fit3)[20] 
logNDE <- coef(modY)[20]
RDR <- 100 * (logTE - logNDE) / logNDE # -11.8%

# With gate variable
modY <- clogit(update(base, ~. + I(Qe_Alc/12) + Qe_Alc_cat + M4), data = crc1) 
modM <- lm(M4 ~ Bmi_C + Qe_Energy + L_School + Smoke_Stat + Smoke_Int + Height_C + 
             Qge0701 + I(Qe_Alc/12) + Qe_Alc_cat, data = crc1)

# Natural direct effect is given as exp(coeff for 1 unit change in exposure) theta1
theta1 <- tidy(modY, exponentiate = T, conf.int = T) %>% filter(term == "I(Qe_Alc/12)") %>% select(-1)
theta2 <- tidy(modY, conf.int = T) %>% filter(term == "M4") %>% select(-1)
beta1 <- tidy(modM, conf.int = T) %>% filter(term == "I(Qe_Alc/12)") %>% select(-1)
NIE <- exp(theta2 * beta1)
# TE = 1.18 (1.00-1.34), NDE = 1.20 (1.05-1.37), NIE = 0.99 (0.98-1.01)

logTE <- coef(fit4)[20] 
logNDE <- coef(modY)[20]
RDR <- 100 * (logTE - logNDE) / logNDE # -9.5%

### Mediator derived from fatty acids data
modY <- clogit(update(base, ~. + I(Qe_Alc/12) + M5), data = crc3) 
modM <- lm(M5 ~ Bmi_C + Qe_Energy + L_School + Smoke_Stat + Smoke_Int + Height_C + 
             Qge0701 + I(Qe_Alc/12), data = crc3)
# Natural direct effect is given as exp(coeff for 1 unit change in exposure) theta1
theta1 <- tidy(modY, exponentiate = T, conf.int = T) %>% filter(term == "M5") %>% select(-1)
theta2 <- tidy(modY, conf.int = T) %>% filter(term == "M5") %>% select(-1)
beta1 <- tidy(modM, conf.int = T) %>% filter(term == "I(Qe_Alc/12)") %>% select(-1)
NIE <- exp(theta2 * beta1)
# TE = 1.18 (1.00-1.34), NDE = 1.21 (1.06-1.38), NIE = 0.95 (0.93-1.00)
logTE <- coef(fit5)[19] %>% as.numeric
logNDE <- coef(modY)[19] %>% as.numeric
RDR <- 100 * (logTE - logNDE) / logNDE # -29.1%

# Gate variable
modY <- clogit(update(base, ~. + I(Qe_Alc/12) + Qe_Alc_cat + M5), data = crc3) 
modM <- lm(M5 ~ Bmi_C + Qe_Energy + L_School + Smoke_Stat + Smoke_Int + Height_C + 
             Qge0701 + I(Qe_Alc/12) + Qe_Alc_cat, data = crc3)
# Natural direct effect is given as exp(coeff for 1 unit change in exposure) theta1
theta1 <- tidy(modY, exponentiate = T, conf.int = T) %>% filter(term == "I(Qe_Alc/12)") %>% select(-1)
theta2 <- tidy(modY, conf.int = T) %>% filter(term == "M5") %>% select(-1)
beta1 <- tidy(modM, conf.int = T) %>% filter(term == "I(Qe_Alc/12)") %>% select(-1)
NIE <- exp(theta2 * beta1)
# TE = 1.18 (1.00-1.34), NDE = 1.21 (1.06-1.38), NIE = 0.95 (0.93-1.00)
logTE <- coef(fit6)[19] 
logNDE <- coef(modY)[19]
RDR <- 100 * (logTE - logNDE) / logNDE # -9.5%
