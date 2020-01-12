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
base <- CT ~ eth + BMI + SMK + DIABETE + RTH + ALCOHOL + DURTHSDIAG + CENTTIME + STOCKTIME + strata(MATCH)
fit1 <- clogit(update(base, . ~ .), data = meta)
fit2 <- clogit(update(base, . ~ .), data = meta1)

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

library(broom)
t <- map_df(list(fit1, fit2), tidy) %>% filter(str_detect(term, "eth"))

par(mar=c(5,4,2,2))
library(metafor)
forest(t$estimate, ci.lb = t$conf.low, ci.ub = t$conf.high, refline = 1,
       xlab = "OR per SD increase in concentration", 
       transf = exp, efac = 0.5,
       pch = 18, cex = 0.8, psize = 1.5, 
       annosym = c("  (", " to ", ")"))

t1 <- map_df(fit1, tidy) %>% filter(str_detect(term, "eth"))

forest(t1$estimate, ci.lb = t1$conf.low, ci.ub = t1$conf.high, refline = 1,
       xlab = "OR per SD increase in concentration", 
       transf = exp, efac = 0.5,
       pch = 18, cex = 0.8, psize = 1.5, 
       annosym = c("  (", " to ", ")"))