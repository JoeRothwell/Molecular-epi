# Separate studies for Jelena's paper. p180 (small) had 740 subjects, p150 (large) had 556
# Note: no NAs and no zeros in these data
source("CRC_prep_data.R")
rm(list = ls(pattern = "crc|dist|prox|rect"))

mat0 <- colon %>% select(contains("Aminoacid_")) %>% 
  select_if(~ sum(is.na(.)) < (nrow(colon)*0.31))

mat <- colon %>% select(colnames(mat0)) #738, 698 w/o Greece
scalemat <- scale(mat)

### Continuous models per SD increase
# Define function to apply across quartiles (already matched by lab)

library(survival)
# Original co-variates. Breslow needed for >10 years
multiclr <- function(AAconc, dat) { 
  clogit(Cncr_Caco_Clrt ~ AAconc + Bmi_Cat + Smoke_Stat + Alc_Drinker + Pa_Index + #Pa_Total + 
           strata(Match_Caseset), method = "breslow", 
           data = dat) 
}

library(broom)
mods <- apply(scalemat, 2, multiclr, dat = colon) %>% 
  map_df( ~ tidy(., exponentiate = T, conf.int = T)) %>% filter(term == "AAconc") %>% 
  mutate(p.adj = p.adjust(p.value, method = "fdr"), compound = colnames(mat), Phase = "Phase I")

results <- mods %>% select(Phase, compound, estimate, conf.low, conf.high, p.adj) 

# Format OR and 95% CI
pooled <- results %>%
  mutate_at(vars(estimate:conf.high), ~ format(round(., digits = 2), nsmall = 2)) %>%
  #mutate_at(vars(p.adj:phet), ~ round(., digits = 3)) %>% 
  mutate(B1 = " (", hyph = "-", B2 = ")") %>%
  unite("CI", estimate, B1, conf.low, hyph, conf.high, B2, sep = "")

# Using map() (not used)
contD <- bind_rows(cont1, cont2) %>% group_by(compound) %>% nest() %>% 
  mutate(mods = map(data, ~ rma(estimate, sei = std.error, data = ., method="REML")),
         tidied = map(mods, tidy, conf.int = T)) %>% unnest(tidied)

# Format sensitivity analyses
followup <- MA.cont %>%
  mutate_at(vars(estimate:conf.high), ~ format(round(., digits = 2), nsmall = 2)) %>%
  mutate_at(vars(phet), ~ round(., digits = 3)) %>% 
  mutate(B1 = " (", hyph = "-", B2 = ")") %>%
  unite("CI", estimate, B1, conf.low, hyph, conf.high, B2, sep = "")


### Categorical analysis. Need to alter function to get categories with cutpoints based 
# on controls
ctrl <- colon$Cncr_Caco_Clrt == 0

# Using control cutpoints only
#cutct <- function(x) cut(x, breaks = quantile(x[ctrl]), include.lowest = T, labels = 1:4)

# Using control inner cutpoints and full range outer cutpoints
cutct <- function(x, ...) { 
  inner <- quantile(x[ctrl])[2:4]
  outer <- quantile(x)[c(1,5)]
  cut(x, breaks = sort(c(outer, inner)), include.lowest = T, ...)
}

# Apply across continuous matrix to get categories as factor or integer (for p-trend)
mat1 <- apply(mat, 2, cutct, labels = 1:4)
mat1a <- apply(mat, 2, cutct, labels = F)

library(broom)
fits1 <- apply(mat1, 2, multiclr, dat = colon) %>% 
  map_df( ~tidy(., exponentiate = T, conf.int = T)) %>%
  filter(grepl("AAconc", term)) %>%
  mutate(Phase = "Phase I", compound = rep(colnames(mat), each = 3)) %>%
  mutate(p.adj = p.adjust(p.value, method = "fdr"))

# Get p-trend by entering quartile numbers into the model as continuous variables (1,2,3,4)
trends <- apply(mat3a, 2, multiclr, dat = colon1) %>% map_df( ~ tidy(.)) %>% 
  filter(term == "AAconc") %>% select(p.value)

# Get control subset for phase II and apply cutpoints across mat2
ctrl <- colon2$Cncr_Caco_Clrt == 0
mat4 <- apply(mat2, 2, cutct, labels = 1:4)
mat4a <- apply(mat2, 2, cutct, labels = F)

fits2 <- apply(mat4, 2, multiclr, dat = colon2) %>% 
  map_df( ~tidy(., exponentiate = T, conf.int = T)) %>% 
  filter(grepl("AAconc", term)) %>%
  mutate(Phase = "Phase II", compound = rep(colnames(mat2), each = 3)) %>%
  mutate(p.adj = p.adjust(p.value, method = "fdr"))
  
results3 <- fits1 %>% select(Phase, term, compound, estimate, conf.low, conf.high, p.adj) 
results4 <- fits2 %>% select(Phase, term, compound, estimate, conf.low, conf.high, p.adj)

# P-trends
trends <- apply(mat4a, 2, multiclr, dat = colon2) %>% map_df( ~ tidy(.)) %>% 
  filter(term == "AAconc") %>% select(p.value)

# Meta-analysis
library(metafor)
cats <- bind_rows(fits1, fits2) %>% group_by(compound, term) %>% nest() %>% 
  mutate(mods = lapply(data, function(df) rma(estimate, sei = std.error, data=df, method="REML")))

# Extract data from models and make data frame
ma.cat <- lapply(cats$mods, "[", c("b", "ci.lb", "ci.ub", "se", "I2", "QEp"))
MA.cat <- map_df(ma.cat, bind_rows) %>%
  bind_cols(compound = fits1$compound, term = cats$term) %>%
  select(compound, term, estimate = "b", conf.low = "ci.lb", conf.high = "ci.ub", I2, phet = "QEp")


tabcat <- bind_rows(results3, results4, MA.cat) %>%
  mutate_at(vars(estimate:conf.high), ~ format(round(., digits = 2), nsmall = 2)) %>% 
  mutate_at(vars(p.adj:phet), ~ round(., digits = 3)) %>% 
  mutate(B1 = " (", hyph = "-", B2 = ")") %>%
  unite("OR", estimate, B1, conf.low, hyph, conf.high, B2, sep = "") %>%
  arrange(compound, term)

# Get a table for each quartile with analyses and meta-analysis
Q2 <- filter(tabcat, term == "AAconc2") %>% select(OR:phet)
Q3 <- filter(tabcat, term == "AAconc3") %>% select(OR:phet)
Q4 <- filter(tabcat, term == "AAconc4") %>% select(OR:phet)

# Now copy and paste tabcont, Q2, Q3, Q4 into Excel

# Extra covariates for sensitivity analysis
library(survival)
multiclr <- function(aa.conc, dat) { 
  clogit(Cncr_Caco_Clrt ~ aa.conc + Bmi_Cat + Smoke_Stat + Alc_Drinker + Pa_Total + Qe_Energy + 
           L_School + Smoke_Int + Height_C + Qe_Alc + Qge0701 + Qge0704 + Qge05 +
           strata(Match_Caseset), data = dat) 
}

