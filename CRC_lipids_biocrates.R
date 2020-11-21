# CRC lipids study
# 112 lipids in 3223 participants, no missing values

# Get data from CRC_prep_data_rev and remove unneeded objects
source("CRC_prep_data_rev.R")
rm(list = ls(pattern = "1|2|bmi"))

# Join GRS data
snps <- read_dta("clrt_gwas_gecco_snps_GRS.dta")
csnp <- inner_join(crc, snps, by = "Idepic")
table(CT = csnp$Cncr_Caco_Clrt, lab = csnp$lab)
table(CT = csnp$Cncr_Caco_Clrt)

library(zoo)

# Continuous models (needs imputation and log transformation)
# Remove non-matched subjects
#crc1 <- crc %>% group_by(Match_Caseset) %>% filter(n() == 2) %>% ungroup()


# Select the lipids only (112) and convert zeros to NA
mat1 <- crc %>% select(Acylcarn_C10:Acylcarn_C8, Glyceroph_Lysopc_A_C16_0:Sphingo_Sm_C26_1) %>% na_if(0)
mat2 <- colon %>% select(Acylcarn_C10:Acylcarn_C8, Glyceroph_Lysopc_A_C16_0:Sphingo_Sm_C26_1) %>% na_if(0)
mat3 <- prox %>% select(Acylcarn_C10:Acylcarn_C8, Glyceroph_Lysopc_A_C16_0:Sphingo_Sm_C26_1) %>% na_if(0)
mat4 <- dist %>% select(Acylcarn_C10:Acylcarn_C8, Glyceroph_Lysopc_A_C16_0:Sphingo_Sm_C26_1) %>% na_if(0)
mat5 <- rectal %>% select(Acylcarn_C10:Acylcarn_C8, Glyceroph_Lysopc_A_C16_0:Sphingo_Sm_C26_1) %>% na_if(0)

# Check for missings
hist(colSums(is.na(mat)), breaks = 30)
scalemat1 <- mat1 %>% na.aggregate(FUN = function(x) min(x)/2) %>% log2 %>% scale
scalemat2 <- mat2 %>% na.aggregate(FUN = function(x) min(x)/2) %>% log2 %>% scale
scalemat3 <- mat3 %>% na.aggregate(FUN = function(x) min(x)/2) %>% log2 %>% scale
scalemat4 <- mat4 %>% na.aggregate(FUN = function(x) min(x)/2) %>% log2 %>% scale
scalemat5 <- mat5 %>% na.aggregate(FUN = function(x) min(x)/2) %>% log2 %>% scale

library(Amelia)
missmap(as_tibble(scalemat))
heatmap.2(scalemat1, trace = "none", col = colpalette1)

library(corrplot)
cormat <- cor(scalemat1)
corrplot(cormat, method = "color", order = "hclust", tl.col = "black", tl.cex = 0.7)

### Continuous models per SD increase
# Define function to apply across quartiles (already matched by lab)
library(survival)
multiclr <- function(x, dat) { 
  clogit(Cncr_Caco_Clrt ~ x + Bmi_C + Qe_Energy + L_School + Smoke_Stat + Smoke_Int + Height_C +
    Qe_Alc + Qge0701 + Qge0704 + strata(Match_Caseset), data = dat) 
  }

# Apply models by subsite (use 2nd fn parameter as an option)
mods1 <- apply(scalemat1, 2, multiclr, dat = crc) %>% map_df( ~ tidy(., exponentiate = T)) %>% 
  filter(grepl("x", term)) %>% mutate(p.adj = p.adjust(p.value, method = "fdr"))

mods2 <- apply(scalemat2, 2, multiclr, dat = colon) %>% map_df( ~ tidy(., exponentiate = T)) %>% 
  filter(grepl("x", term)) %>% mutate(p.adj = p.adjust(p.value, method = "fdr"))

mods3 <- apply(scalemat3, 2, multiclr, dat = prox) %>% map_df( ~ tidy(., exponentiate = T)) %>% 
  filter(grepl("x", term)) %>% mutate(p.adj = p.adjust(p.value, method = "fdr"))

mods4 <- apply(scalemat4, 2, multiclr, dat = dist) %>% map_df( ~ tidy(., exponentiate = T)) %>% 
  filter(grepl("x", term)) %>% mutate(p.adj = p.adjust(p.value, method = "fdr"))

mods5 <- apply(scalemat5, 2, multiclr, dat = rectal) %>% map_df( ~ tidy(., exponentiate = T)) %>% 
  filter(grepl("x", term)) %>% mutate(p.adj = p.adjust(p.value, method = "fdr"))


#crccon <- map_df(mods, ~ tidy(., exponentiate = T)) %>% filter(grepl("x", term)) %>% 
#  mutate_if(is.numeric, ~round(., 4)) %>% mutate(p.adj = p.adjust(p.value, method = "fdr")) %>%
#  unite(OR.CI, estimate, conf.low, conf.high, sep = "-") %>% cbind(cmpd = colnames(scalemat))


### Categorical associations (doesn't need imputation or scaling)
# First need to run CRC_data_prep to get crc and crc3 objects
library(ggplot)
library(broom)

# Biocrates compounds. Split into quartiles with cut_number
mat6 <- mat1 %>% mutate_all(~cut_number(., n = 4, labels = 1:4))
mat7 <- mat2 %>% mutate_all(~cut_number(., n = 4, labels = 1:4))
mat8 <- mat3 %>% mutate_all(~cut_number(., n = 4, labels = 1:4))
mat9 <- mat4 %>% mutate_all(~cut_number(., n = 4, labels = 1:4))
mat10 <- mat5 %>% mutate_all(~cut_number(., n = 4, labels = 1:4))


fits1 <- apply(mat6, 2, multiclr, dat = crc) %>% map_df( ~tidy(., exponentiate = T)) %>% 
  filter(grepl("x", term)) %>% mutate_if(is.numeric, ~round(., 3)) %>%
  mutate(p.adj = p.adjust(p.value, method = "fdr")) %>% 
  #unite(OR.CI, estimate, conf.low, conf.high, sep = "-") %>%
  cbind(cmpd = colnames(df3)) %>% group_by(cmpd) %>% filter(min(p.value) < 0.05)
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
