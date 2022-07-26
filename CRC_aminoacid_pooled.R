# EPIC P180 and P150 studies now pooled instead of meta-analysed. 
# p180 (small) had 740 subjects, p150 (large) had 556. Note: no NAs and no zeros in these data
source("CRC_prep_data.R")
rm(list = ls(pattern = "colon|dist|prox|rect|crc1|crc2"))
# save.image(file = "amino_acids_models.Rdata")

library(survival)
# Define function to apply across quartiles (already matched by lab)
# Original co-variates. Breslow needed for >10 years
multiclr <- function(AAconc, dat) { 
  clogit(Cncr_Caco_Clrt ~ AAconc + Bmi_Cat + Smoke_Stat + Alc_Drinker + Pa_Index + #Pa_Total + 
           strata(Match_Caseset), method = "breslow", 
         data = dat) 
}

# Subsets for sensitivity analysis
# Follow up time
crc <- crc %>% filter(Tfollowup >= 2)
crc <- crc %>% filter(Tfollowup >= 5)
crc <- crc %>% filter(Tfollowup >= 10)

# Tumour stage, EPIC classification
crc <- crc %>% filter(Stagclrt == 2)
crc <- crc %>% filter(Stagclrt %in% 3:5)

# By sex
crc <- crc %>% filter(Sex == 1)
crc <- crc %>% filter(Sex == 2)

# Count amino acids
crc %>% select(contains("Aminoacid_")) %>% ncol() #22

# Filter amino acids with over 31% missings
# Update: keep all amino acids and give missings instead (to keep p180 AAs)
mat0 <- crc %>% select(contains("Aminoacid_")) #%>% 
  #select_if(~ sum(is.na(.)) < (nrow(crc) * 0.31))
mat <- crc %>% select(colnames(mat0)) #1308 w/o Greece
scalemat <- scale(mat)

### Continuous models per SD increase
library(broom)
mods <- apply(scalemat, 2, multiclr, dat = crc) %>% 
  map_df( ~ tidy(., exponentiate = T, conf.int = T)) %>% filter(term == "AAconc") %>% 
  mutate(p.adj = p.adjust(p.value, method = "fdr"), compound = colnames(mat))

# Count non-missings for each amino acid
samp <- tibble(measured = colSums(!is.na(mat)))

results <- mods %>% select(compound, estimate, conf.low, conf.high, p.value, p.adj) %>% bind_cols(samp)

# Format OR (to convert to function?)
pooled <- results %>%
  mutate_at(vars(estimate:conf.high), ~ format(round(., digits = 2), nsmall = 2)) %>%
  mutate(B1 = " (", hyph = "-", B2 = ")") %>%
  unite("CI", estimate, B1, conf.low, hyph, conf.high, B2, sep = "")

### Categorical analysis. Need to alter function to get categories with cutpoints based on controls
ctrl <- crc$Cncr_Caco_Clrt == 0

# Using control cutpoints only
#cutct <- function(x) cut(x, breaks = quantile(x[ctrl]), include.lowest = T, labels = 1:4)

# Using control inner cutpoints and full range outer cutpoints
cutct <- function(x, ...) { 
  inner <- quantile(x[ctrl])[2:4]
  outer <- quantile(x)[c(1, 5)]
  cut(x, breaks = sort(c(outer, inner)), include.lowest = T, ...)
}

# Apply across continuous matrix to get categories as factor or integer (for p-trend)
mat1 <- apply(mat, 2, cutct, labels = 1:4)
mat1a <- apply(mat, 2, cutct, labels = F)

library(broom)
fits1 <- apply(mat1, 2, multiclr, dat = crc) %>% 
  map_df( ~tidy(., exponentiate = T, conf.int = T)) %>%
  filter(grepl("AAconc", term)) %>%
  mutate(compound = rep(colnames(mat1), each = 3)) %>%
  mutate(p.adj = p.adjust(p.value, method = "fdr"))

# Get p-trend by entering quartile numbers into the model as continuous variables (1,2,3,4)
trends <- apply(mat1a, 2, multiclr, dat = crc) %>% map_df( ~ tidy(.)) %>% 
  filter(term == "AAconc") %>% select(p.value)

# Extract results
results1 <- fits1 %>% select(term, compound, estimate, conf.low, conf.high, p.adj) 

tabcat <- results1 %>%
  mutate_at(vars(estimate:conf.high), ~ format(round(., digits = 2), nsmall = 2)) %>% 
  mutate_at(vars(p.adj), ~ round(., digits = 3)) %>% 
  mutate(B1 = " (", hyph = "-", B2 = ")") %>%
  unite("OR", estimate, B1, conf.low, hyph, conf.high, B2, sep = "") %>%
  arrange(compound, term)

# Get a table for each quartile with analyses and meta-analysis
Q2 <- filter(tabcat, term == "AAconc2") %>% select(OR:p.adj)
Q3 <- filter(tabcat, term == "AAconc3") %>% select(OR:p.adj)
Q4 <- filter(tabcat, term == "AAconc4") %>% select(OR:p.adj)

# Categorical analysis for P180 compounds only (run cutct first)
crc180 <- crc %>% filter(lab == 1)
ctrl <- crc180$Cncr_Caco_Clrt == 0
# Remove aspartate (too many of the same value) and Xleu (not included)
mat180 <- crc180 %>% select(colnames(mat0[-c(4,22)])) #708 w/o Greece

# Apply across continuous matrix to get categories as factor or integer labels = F (for p-trend)
mat1 <- apply(mat180, 2, cutct, labels = 1:4)
mat1a <- apply(mat180, 2, cutct, labels = F)

library(broom)
fits180 <- apply(mat1, 2, multiclr, dat = crc180) %>% 
  map_df( ~tidy(., exponentiate = T, conf.int = T)) %>% filter(grepl("AAconc", term)) %>%
  mutate(compound = rep(colnames(mat180), each = 3)) %>%
  mutate(p.adj = p.adjust(p.value, method = "fdr"))

# Get p-trend by entering quartile numbers into the model as continuous variables (1,2,3,4)
trends <- apply(mat1a, 2, multiclr, dat = crc180) %>% map_df( ~ tidy(.)) %>% 
  filter(term == "AAconc") %>% select(p.value) %>% bind_cols(compound = colnames(mat180))

# Extract results
results1 <- fits180 %>% select(term, compound, estimate, conf.low, conf.high, p.adj) 

tabcat <- results1 %>%
  mutate_at(vars(estimate:conf.high), ~ format(round(., digits = 2), nsmall = 2)) %>% 
  mutate_at(vars(p.adj), ~ round(., digits = 3)) %>% 
  mutate(B1 = " (", hyph = "-", B2 = ")") %>%
  unite("OR", estimate, B1, conf.low, hyph, conf.high, B2, sep = "") %>%
  arrange(compound, term)

# Correlations
mat <- colon %>% filter(Cncr_Caco_Clrt == 0) %>% select(contains("Aminoacid_")) %>% 
  select_if(~ sum(is.na(.)) < (nrow(colon)*0.31)) %>% select(-Aminoacid_Xleu)
scalemat <- scale(mat)
library(corrplot)
cormat <- cor(scalemat, use = "pairwise.complete.obs")
corrplot(cormat, method = "square", order = "hclust")


