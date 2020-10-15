load("pred_score_tables_rev.Rdata")


# Alcohol-CRC associations
library(survival)
library(sjmisc)
base <- Cncr_Caco_Clrt ~ Bmi_C + Qe_Energy + L_School + Smoke_Stat + Smoke_Int + Height_C + 
  Qge0701 + strata(Match_Caseset)

# Both case-control sets, biocrates set, fatty acids set
fit1 <- clogit(update(base, ~. + lab + I(Qe_Alc/12)), data = crc.ph)
fit3 <- clogit(update(base, ~. + I(Qe_Alc/12)), data = crc1)
#fit5 <- clogit(update(base, ~. + I(Qe_Alc/12)), data = crc3.ph)

crc.ph$Qe_Alc_cat <- dicho(crc.ph$Qe_Alc, dich.by = 0.5)
crc1$Qe_Alc_cat <- dicho(crc1$Qe_Alc, dich.by = 0.5)
crc3$Qe_Alc_cat <- dicho(crc3$Qe_Alc, dich.by = 0.5)

# Binary variable for alcohol intake, cutoff 0.5
fit2 <- clogit(update(base, ~. + lab + Qe_Alc_cat), data = crc.ph)
fit4 <- clogit(update(base, ~. + Qe_Alc_cat), data = crc1)
#fit6 <- clogit(update(base, ~. + Qe_Alc_cat), data = crc3)

# Excluding zero consumption made association weaker

modlist <- list(fit1, fit2, fit3, fit4, fit5, fit6)
modlist <- list(fit1, fit2, fit3, fit4)

library(broom)
scomods <- map_df(modlist, ~tidy(., exponentiate = T)) %>% filter(str_detect(term, "Alc")) %>%
  mutate_if(is.numeric, ~round(., 2)) %>% unite(OR.CI, estimate, conf.low, conf.high, sep = "-")

#fit2 <- clogit(update(base, ~. + scale(Qe_Alc)), data = crc1)
#fit3 <- clogit(update(base, ~. + cut_number(Qe_Alc, n = 4)), data = crc1)
#fit2 <- clogit(update(base, ~. + lab + scale(Qe_Alc)), data = crc.ph)
#fit3 <- clogit(update(base, ~. + lab + cut_number(Qe_Alc, n = 4)), data = crc.ph)

# Prepare data for PLS
library(pls)
expr <- "(carn|oacid|genic|roph|ingo|Sugars)[_]"
crc1p <- crc1 %>% select(matches(expr), -contains("tdq")) %>% select_if(~ sum(., na.rm = T) != 0)
adj   <- function(x) residuals(lmer(x ~ Center + batch_no + Sex + (1|Study), data = ctrl))