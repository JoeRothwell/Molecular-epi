# Preparation of data for BC study (after exploration of data)
# 2 datasets for continuous and categorical analysis

library(tidyverse)
library(readxl)

# First get metadata. Subset variables needed
meta <- read_csv("Lifepath_meta.csv", na = "9999") %>%
  select(CT, BMI, SMK, DIABETE, RTH, ALCOHOL, DURTHSDIAG, CENTTIME, STOCKTIME, RACK, MATCH, MENOPAUSE) %>%
  mutate_at(vars(SMK, DIABETE, RACK, MATCH), as.factor)

cmpd_meta <- read.csv("NMR_cmpd_metadata_new.csv")

# For removal of problem racks for Ethanol
meta1 <- meta %>% filter(!RACK %in% c(10, 29, 33, 34))

# Continuous. Read scaled data and subset to get subjects included in CC

ints0 <- read_tsv("1510_XMetaboliteE3N_cpmg_unscaled.txt")

#ints.all <- read.delim("1507_XMetabolite_std_cpmg_E3N.txt")
#samples <- ints.all$CODBMB %in% meta$CODBMB
#ints <- ints.all[samples, -1]

# Remove problem racks
filt <- !(meta$RACK %in% c(10, 29, 33, 34))
ints.filt <- ints0[filt, ]

# Remove top and bottom 1% of each compound and replace with NA
outliers <- function(x) x > quantile(x, probs = 0.99) | x < quantile(x, probs = 0.01)
logicalmat <- apply(ints.con, 2, outliers)
ints.filt[logicalmat] <- NA

# Scale to unit variance
ints <- scale(ints0)

# Replace negative values with half the minimum positive value
#rm.neg.values <- function(x) ifelse(x < 0, min(x[x > 0])/2, x)
#ints0 <- apply(ints0, 2, rm.neg.values)
# Check
#which(apply(as.matrix(ints0), 2, min) < 0)

# Categorical

# Update 11/10/19: Updated with unscaled data (see BC_compounds_scaling for old data)
#ints0 <- read_tsv("1510_XMetaboliteE3N_cpmg_unscaled.txt")
quartiles <- ints.filt %>% mutate_all(funs(cut_number(., n = 4, labels = 1:4))) 
dat1 <- cbind(meta[filt, ], quartiles)
