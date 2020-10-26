# CRC lipids study
# 112 lipids in 3223 participants, no missing values


# Continuous models (needs imputation and log transformation)
# Remove non-matched subjects
crc1 <- crc %>% group_by(Match_Caseset) %>% filter(n() == 2) %>% ungroup()

scalemat <- crc1 %>%
  #group_by(Match_Caseset) %>% filter(n() == 2) %>% ungroup() %>%
  select(Acylcarn_C10:Acylcarn_C8, Glyceroph_Lysopc_A_C16_0:Sphingo_Sm_C26_1) %>%
  log2 %>% scale

missmap(as_tibble(scalemat))

# Continuous models per SD increase
base <- Cncr_Caco_Clrt ~ Qe_Energy + L_School + Smoke_Stat + Smoke_Int + Height_C + 
  strata(Match_Caseset)

library(survival)
clr <- function(x) { clogit(Cncr_Caco_Clrt ~ x + Height_C + Qe_Energy + L_School + Smoke_Stat + Smoke_Int + 
                            strata(Match_Caseset), 
                        data = crc1) }

apply(scalemat, 2, clr)




# Categorical associations (doesn't need imputation or scaling)
# First need to run CRC_data_prep to get crc and crc3 objects
library(ggplot)
library(broom)

# Biocrates compounds. Split into quartiles with cut_number
df3 <- crc %>% select(Acylcarn_C10:Acylcarn_C8, Glyceroph_Lysopc_A_C16_0:Sphingo_Sm_C26_1) %>% 
  mutate_all(~cut_number(., n = 4, labels = 1:4))

# Define function to apply across quartiles (already matched by lab)
clrfun <- function(x)  { clogit(Cncr_Caco_Clrt ~ x + Bmi_C + Qe_Energy + L_School + Smoke_Stat + Smoke_Int + Height_C +
                                  Qe_Alc + Qge0701 + Qge0704 + strata(Match_Caseset), data = crc) }

fits2 <- apply(df3, 2, clrfun)
mods2 <- map_df(fits2, ~tidy(., exponentiate = T)) %>% filter(grepl("x", term)) %>% 
  mutate_if(is.numeric, ~round(., 2)) %>%
  unite(OR.CI, estimate, conf.low, conf.high, sep = "-") %>%
  cbind(cmpd = colnames(df3)) %>% group_by(cmpd) %>%
  filter(min(p.value) < 0.05)
  arrange(cmpd)



