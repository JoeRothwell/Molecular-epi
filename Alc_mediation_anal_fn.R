source("Make_alc_mediators.R")
library(haven)
library(survival)
library(sjmisc)

# Alcohol-CRC associations
base <- Cncr_Caco_Clrt ~ Bmi_C + Qe_Energy + L_School + Smoke_Stat + Smoke_Int + Height_C + 
  Qge0701 + strata(Match_Caseset)

# Make dichotomised "gate" variables
crc.ph$Qe_Alc_cat  <- dicho(crc.ph$Qe_Alc, dich.by = 0.5)
crc3.ph$Qe_Alc_cat <- dicho(crc3.ph$Qe_Alc, dich.by = 0.5)
crc1$Qe_Alc_cat    <- dicho(crc1$Qe_Alc, dich.by = 0.5)

# Total effects for alcohol, both case-controls without and with gate variable
#fit1 <- clogit(update(base, ~. + I(Qe_Alc/12)), data = crc.ph) # OR 1.06 (1.00-1.11)
#fit2 <- clogit(update(base, ~. + I(Qe_Alc/12) + Qe_Alc_cat), data = crc.ph) # OR 1.07 (1.01-1.13)

# Small case-control without and with gate variable
fit3 <- clogit(update(base, ~. + I(Qe_Alc/12)), data = crc1) # OR 1.14 (1.00-1.28)
fit4 <- clogit(update(base, ~. + I(Qe_Alc/12) + Qe_Alc_cat), data = crc1) # OR 1.18 (1.03-1.34)

# Fatty acids case-control without and with gate variable
fit5 <- clogit(update(base, ~. + I(Qe_Alc/12)), data = crc3.ph) # OR 1.15 (1.00-1.29)
fit6 <- clogit(update(base, ~. + I(Qe_Alc/12) + Qe_Alc_cat), data = crc3.ph) # OR 1.19 (1.04-1.35)

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
mediate <- function(fit, dat, M, gate = F){
  
  # M1: controls from small case-control, all biocrates compounds (as Laura did)
  # First get total effect (OR and CIs and just coefficient)
  TE <- tidy(fit, conf.int = T, exponentiate = T) %>% filter(term == "I(Qe_Alc/12)")
  ln.te <- tidy(fit) %>% filter(term == "I(Qe_Alc/12)") %>% select(2)
  
  # modY is a CLR adjusting for the mediator. Need exposure and mediator coefs theta1, theta2 
  if(gate == F) {
  modY <- clogit(Cncr_Caco_Clrt ~ Bmi_C + Qe_Energy + L_School + Smoke_Stat + Smoke_Int + 
                   Height_C + Qge0701 + strata(Match_Caseset) + I(Qe_Alc/12) + M, data = dat)
  } else {
  modY <- clogit(Cncr_Caco_Clrt ~ Bmi_C + Qe_Energy + L_School + Smoke_Stat + Smoke_Int + 
                   Height_C + Qge0701 + strata(Match_Caseset) + I(Qe_Alc/12) + Qe_Alc_cat + M, data = dat)   
  }
  
  # modM is a linear model of outcome and mediator. Need exposure coef beta1
  if(gate == F) {
  modM <- lm(M ~ Bmi_C + Qe_Energy + L_School + Smoke_Stat + Smoke_Int + Height_C + 
               Qge0701 + I(Qe_Alc/12), data = dat)
  } else {
  modM <- lm(M ~ Bmi_C + Qe_Energy + L_School + Smoke_Stat + Smoke_Int + Height_C + 
                 Qge0701 + I(Qe_Alc/12) + Qe_Alc_cat, data = dat)
  }
  
  # Natural direct effect is given as exp(coeff for 1 unit change in exposure), theta1
  NDE <- tidy(modY, exponentiate = T, conf.int = T) %>% filter(term == "I(Qe_Alc/12)")
  theta1 <- NDE %>% select(-1)
  
  # Natural indirect effect is exp(coeff of theta for mediator x beta for exposure)
  # Get Theta coefficients of mod Y and beta coefficients of modM:
  theta2 <- tidy(modY, conf.int = T) %>% filter(term == "M") %>% select(-1)
  beta1  <- tidy(modM, conf.int = T) %>% filter(term == "I(Qe_Alc/12)") %>% select(-1)
  NIE <- exp(theta2*beta1)
  ln.nie <- (theta2*beta1)[1]

  # Risk difference ratio is defined as log(TE) - log(NDE) / log (NDE)  *100
  logTE <- coef(fit)["I(Qe_Alc/12)"] 
  logNDE <- coef(modY)["I(Qe_Alc/12)"]
  RDR <- 100 * (logTE - logNDE) / logNDE # 90.6
  med <- (100 * ln.nie)/ln.te
  
  #return(list(te = TE, nde = NDE, nie = NIE, rdr = RDR))
  return(list(te.nde = bind_rows(TE, NDE), nie = NIE, rdr = RDR, perc = med))
}

# M1: controls from small case-control, all biocrates compounds (as Laura did)
mediate(fit3, crc1, M1, gate = F)
mediate(fit4, crc1, M1, gate = T)
# TE = 1.14 (1.00-1.28), NDE = 1.07 (0.94-1.22), NIE = 1.10 (1.02-1.23), RDR = 90.6%
# TE = 1.18 (1.03-1.34), NDE = 1.11 (0.96-1.27), NIE = 1.10 (1.01-1.24), RDR = 61.6%

# M2: controls from small study, amino acids only
mediate(fit3, crc1, M2, gate = F)
mediate(fit4, crc1, M2, gate = T)
# TE = 1.14 (1.00-1.28), NDE = 1.14 (1.01-1.29), NIE = 0.99 (1.00-1.00), RDR = -3.90%
# TE = 1.18 (1.03-1.34), NDE = 1.19 (1.04-1.35), NIE = 0.99 (1.00-1.00), RDR = -3.97%

# M3: controls from small study, lipid metabolites only
mediate(fit3, crc1, M3, gate = F)
mediate(fit4, crc1, M3, gate = T)
# TE = 1.14 (1.00-1.28), NDE = 1.06 (0.93-1.21), NIE = 1.10 (1.02-1.21), RDR = 118.6
# TE = 1.18 (1.03-1.34), NDE = 1.10 (0.95-1.26), NIE = 1.10 (1.02-1.22), RDR = 78.0

# M4: derived from EPIC non-CRC pooled controls
mediate(fit3, crc1, M4, gate = F)
mediate(fit4, crc1, M4, gate = T)
# TE = 1.14 (1.00-1.28), NDE = 1.16 (1.02-1.31), NIE = 0.99 (0.98-1.01), RDR = -11.8
# TE = 1.18 (1.00-1.34), NDE = 1.20 (1.05-1.37), NIE = 0.99 (0.98-1.01), RDR = -9.5

# M5: derived from fatty acids data
mediate(fit5, crc3.ph, M5, gate = F)
mediate(fit6, crc3.ph, M5, gate = T)
# TE = 1.15 (1.02-1.29), NDE = 1.21 (1.06-1.38), NIE = 0.95 (0.93-1.00), -29.1
# TE = 1.19 (1.04-1.35), NDE = 1.27 (1.10-1.46), NIE = 0.95 (0.93-1.00), -28.6

# Single compound mediators
# LysoPC 16:1 and 17:0, PC aa 32:1, 34:1, ae30:2, 36:2, 38:3, SM 14:1, 16:1, 22:2
mediate(fit3, crc1, M6, gate = F)
mediate(fit3, crc1, M7, gate = F)
mediate(fit3, crc1, M8, gate = F)
mediate(fit3, crc1, M9, gate = F)
mediate(fit3, crc1, M10, gate = F)
mediate(fit3, crc1, M11, gate = F)
mediate(fit3, crc1, M12, gate = F)
mediate(fit3, crc1, M13, gate = F)
mediate(fit3, crc1, M14, gate = F)
mediate(fit3, crc1, M15, gate = F)

# Oxidative stress mediators
mediate(fit3, crc1, M16, gate = F)
mediate(fit3, crc1, M17, gate = F)


