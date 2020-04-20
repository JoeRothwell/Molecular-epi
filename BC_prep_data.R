# Preparation of data for BC study (after exploration of data)
# 2 datasets for continuous and categorical analysis

library(tidyverse)
library(readxl)
library(survival)

# First get metadata. Add follow up time and subset variables needed
meta <- read_csv("Lifepath_meta.csv", na = "9999") %>%
  group_by(MATCH) %>% mutate(Tfollowup = max(DIAGSAMPLING, na.rm = T)) %>% ungroup %>%
  select(CT, BMI, SMK, DIABETE, RTH, DURTHSBMB, CENTTIME, STOCKTIME, MATCH, ALCOHOL, MENOPAUSE,
         DIAGSAMPLING, Tfollowup) %>%
  mutate_at(vars(SMK, DIABETE), as.factor) %>%
  mutate(DURTHSBMBCat = ifelse(DURTHSBMB > 0, 1, 0))

# For subsetting
#pre <- meta$MENOPAUSE == 0
#post <- meta$MENOPAUSE == 1 

# For removal of problem racks for Ethanol
#meta1 <- meta %>% filter(!RACK %in% c(10, 29, 33, 34))

# For removal of case-control pairs not matched by menopausal status
unmatch_pairs <- meta %>% group_by(MATCH) %>% summarise(sum.men = sum(MENOPAUSE)) %>% 
  filter(sum.men == 1) %>% select(MATCH) %>% pull() %>% as.numeric

log.vec <- !(meta$MATCH %in% unmatch_pairs)
meta <- meta[log.vec, ]

# For subsetting
pre <- meta$MENOPAUSE == 0
post <- meta$MENOPAUSE == 1 
pre0 <- meta$MENOPAUSE == 0 & meta$Tfollowup > 2

# Compound data for forest plots
cmpd.meta <- read.csv("NMR_cmpd_metadata_new.csv")
cmpds.ordered <- cmpd.meta %>% arrange(description) %>% mutate(row = 1:n() + (as.numeric(description)-1))
rowvec <- cmpds.ordered$row

# Unscaled data
ints0 <- read_tsv("1510_XMetaboliteE3N_cpmg_unscaled.txt")

#ints.all <- read.delim("1507_XMetabolite_std_cpmg_E3N.txt")
#samples <- ints.all$CODBMB %in% meta$CODBMB
#ints <- ints.all[samples, -1]

# Remove problem samples
#filt <- !(meta$RACK %in% c(10, 29, 33, 34))
ints0 <- ints0[log.vec, ]

# Remove top and bottom 1% of each compound and replace with NA
outliers <- function(x) x > quantile(x, probs = 0.99) | x < quantile(x, probs = 0.01)
logicalmat <- apply(ints0, 2, outliers)
ints0[logicalmat] <- NA

# Scale to unit variance
ints <- scale(ints0)


# Replace negative values with half the minimum positive value
#rm.neg.values <- function(x) ifelse(x < 0, min(x[x > 0])/2, x)
#ints0 <- apply(ints0, 2, rm.neg.values)
# Check
#which(apply(as.matrix(ints0), 2, min) < 0)

# Get a subset of the most discriminating compounds
# Histidine, NAC, glycerol, ornithine, ethanol, pyruvate, albumin, glutamate, glutamine, 3 fatty acids
ss <- ints[, c(5,14,15,16,17,19,20,27,28,33,41,42)]

# Follow up time
boxplot2(meta$DIAGSAMPLING ~ meta$MENOPAUSE, varwidth = T, col = "dodgerblue")






