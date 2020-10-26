# CRC lipids study
# 112 lipids in 3223 participants, no missing values
library(zoo)

# Continuous models (needs imputation and log transformation)
# Remove non-matched subjects
#crc1 <- crc %>% group_by(Match_Caseset) %>% filter(n() == 2) %>% ungroup()

scalemat <- crc %>%
  #group_by(Match_Caseset) %>% filter(n() == 2) %>% ungroup() %>%
  select(Acylcarn_C10:Acylcarn_C8, Glyceroph_Lysopc_A_C16_0:Sphingo_Sm_C26_1) %>%
  na_if(0) %>% na.aggregate(FUN = function(x) min(x)/2) %>%
  log2 %>% scale

missmap(as_tibble(scalemat))
heatmap.2(scalemat, trace = "none", col = colpalette1)

library(corrplot)
cormat <- cor(scalemat)
corrplot(cormat, method = "color", order = "hclust", tl.col = "black", tl.cex = 0.7)

# Continuous models per SD increase
# Define function to apply across quartiles (already matched by lab)
library(survival)
multiclr <- function(x) { 
  clogit(Cncr_Caco_Clrt ~ x + Bmi_C + Qe_Energy + L_School + Smoke_Stat + Smoke_Int + Height_C +
    Qe_Alc + Qge0701 + Qge0704 + strata(Match_Caseset), data = crc) 
  }

mods <- apply(scalemat, 2, multiclr)
mods1 <- map_df(mods, ~ tidy(., exponentiate = T)) %>% filter(grepl("x", term)) %>% 
  mutate_if(is.numeric, ~round(., 4)) %>%
  mutate(p.adj = p.adjust(p.value, method = "fdr")) %>%
  unite(OR.CI, estimate, conf.low, conf.high, sep = "-") %>%
  cbind(cmpd = colnames(scalemat))


# Categorical associations (doesn't need imputation or scaling)
# First need to run CRC_data_prep to get crc and crc3 objects
library(ggplot)
library(broom)

# Biocrates compounds. Split into quartiles with cut_number
df3 <- crc %>% select(Acylcarn_C10:Acylcarn_C8, Glyceroph_Lysopc_A_C16_0:Sphingo_Sm_C26_1) %>% 
  mutate_all(~cut_number(., n = 4, labels = 1:4))

fits2 <- apply(df3, 2, multiclr)
mods2 <- map_df(fits2, ~tidy(., exponentiate = T)) %>% filter(grepl("x", term)) %>% 
  mutate_if(is.numeric, ~round(., 3)) %>%
  mutate(p.adj = p.adjust(p.value, method = "fdr")) %>%
  unite(OR.CI, estimate, conf.low, conf.high, sep = "-") %>%
  cbind(cmpd = colnames(df3)) %>% group_by(cmpd) %>%
  filter(min(p.value) < 0.05)
  arrange(cmpd)
  
  
# Polygenic risk scores
snps <- read_dta("clrt_gwas_gecco_snps_GRS.dta")
csnp <- inner_join(crc, snps, by = "Idepic")
table(CT = csnp$Cncr_Caco_Clrt, lab = csnp$lab)
table(CT = csnp$Cncr_Caco_Clrt)
plot(csnp$GRS)
csnp$GRScat <- cut_number(csnp$GRS, n = 4, labels = 1:4)

table(CT = csnp$Cncr_Caco_Clrt, Q = csnp$GRScat)
chisq.test(csnp$Cncr_Caco_Clrt, csnp$GRScat)

ggplot(csnp, aes(x = GRS, group = as.factor(Cncr_Caco_Clrt))) + 
  geom_density(aes(fill = as.factor(Cncr_Caco_Clrt)), alpha = 0.4)
