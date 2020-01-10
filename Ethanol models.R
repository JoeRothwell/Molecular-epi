# Ethanol models
library(tidyverse)
library(readxl)

ints0 <- read_tsv("1510_XMetaboliteE3N_cpmg_unscaled.txt")
cmpd_meta <- read.csv("NMR_cmpd_metadata_new.csv")

ints <- scale(ints0)
eth <- ints[, 14]

# Lifestyle data. Subset variables needed
meta <- read_csv("Lifepath_meta.csv", na = "9999") %>% mutate_at(vars(SMK, DIABETE, RACK, MATCH, SAMPYEAR), as.factor)
meta <- cbind(meta, eth)
filt <- !(meta$RACK %in% c(10, 29, 33, 34))
meta1 <- meta[filt, ]

library(survival)
fit1 <- clogit(CT ~ eth + BMI + SMK + DIABETE + RTH + ALCOHOL + DURTHSDIAG + CENTTIME + STOCKTIME + strata(MATCH), data = meta)
fit2 <- clogit(CT ~ eth + BMI + SMK + DIABETE + RTH + ALCOHOL + DURTHSDIAG + CENTTIME + STOCKTIME + strata(MATCH), data = meta1)

# categorical

quartiles <- ints0 %>% mutate_all(funs(cut_number(., n = 4, labels = 1:4))) 
eth_q <- cut_number(pull(ints0[, 14]), n = 4, labels = 1:4)
Q1Q4 <- eth_q %in% c(1, 4)
meta <- cbind(meta, eth_q)
meta2 <- meta[Q1Q4, ]
# Need to drop levels 2 and 3 or else model does not converge
meta2$eth_q <- droplevels(meta2$eth_q)

fit3 <- clogit(CT ~ eth_q + BMI + SMK + DIABETE + RTH + ALCOHOL + DURTHSDIAG + CENTTIME + STOCKTIME +
         strata(MATCH), data = meta2)
