# BC risk models for metabolites
library(tidyverse)
library(readxl)

ints0 <- read_tsv("1510_XMetaboliteE3N_cpmg_unscaled.txt")
cmpd_meta <- read.csv("NMR_cmpd_metadata_new.csv")

# Compounds as continuous variables ----
#ints <- scale(ints0)

# Replace negative values with half the minimum positive value
rm.neg.values <- function(x) ifelse(x < 0, min(x[x > 0])/2, x)
ints0 <- apply(ints0, 2, rm.neg.values)

ints <- scale(ints0)

# Lifestyle data. Subset variables needed
meta <- read_csv("Lifepath_meta.csv", na = "9999") %>%
  select(CT, BMI, SMK, DIABETE, RTH, ALCOHOL, DURTHSDIAG, CENTTIME, STOCKTIME, RACK, MATCH, MENOPAUSE, SAMPYEAR) %>%
  mutate_at(vars(SMK, DIABETE, RACK, MATCH, SAMPYEAR), as.factor)

dat <- cbind(meta, ints)

library(survival)
fits <-  apply(ints, 2, function(x) clogit(CT ~ BMI + SMK + DIABETE + RTH + ALCOHOL + 
                                                SAMPYEAR +
                                                DURTHSDIAG + 
                                                CENTTIME +  STOCKTIME + # only have effect on Lactate
                                                #RACK + 
                                                strata(MATCH) + x, data = dat))

library(broom)
t2  <- map_df(fits, tidy) %>% filter(str_detect(term, "x")) %>% bind_cols(cmpd_meta) %>% arrange(description)

# Plot data with Metafor
# Get vectors for row spacings using groups (may add compound classes later)
cmpds_ordered <- cmpd_meta %>% arrange(description) %>% mutate(row = 1:n() + (as.numeric(description) - 1))
rowvec <- cmpds_ordered$row

par(mar=c(5,4,2,2))
library(metafor)
forest(t2$estimate, ci.lb = t2$conf.low, ci.ub = t2$conf.high, refline = 1,
       ylim = c(1, max(rowvec) + 3), #at = 0:5,
       xlab = "OR per SD increase in concentration", 
       transf = exp, rows = rowvec, efac = 0.5,
       pch = 18, cex = 0.8, psize = 1.5, 
       annosym = c("  (", " to ", ")"),
       slab = t2$display_name, 
       main = fits[[1]]$formula[[3]], cex.main = 0.8)

hh <- par("usr")
text(hh[1], max(rowvec) + 2, "Metabolite", pos = 4, cex = 0.8)
text(hh[2], max(rowvec) + 2, "OR [95% CI]", pos = 2, cex = 0.8)

# -----------------------------------------------------------------

fits <-  apply(ints, 2, function(x) clogit(CT ~ BMI + SMK + DIABETE + RTH + ALCOHOL + 
                                             SAMPYEAR +
                                             DURTHSDIAG + 
                                             CENTTIME +  STOCKTIME + # only have effect on Lactate
                                             #RACK + 
                                             strata(MATCH) + x, data = dat, subset = MENOPAUSE == 1))

t2  <- map_df(fits, tidy) %>% filter(str_detect(term, "x")) %>% bind_cols(cmpd_meta) %>% arrange(description)

# Plot data with Metafor
# Get vectors for row spacings using groups (may add compound classes later)
cmpds_ordered <- cmpd_meta %>% arrange(description) %>% mutate(row = 1:n() + (as.numeric(description) - 1))
rowvec <- cmpds_ordered$row

par(mar=c(5,4,2,2))
library(metafor)
forest(t2$estimate, ci.lb = t2$conf.low, ci.ub = t2$conf.high, refline = 1,
       ylim = c(1, max(rowvec) + 3), #at = 0:5,
       xlab = "OR per SD increase in concentration", 
       transf = exp, rows = rowvec, efac = 0.5,
       pch = 18, cex = 0.8, psize = 1.5, 
       annosym = c("  (", " to ", ")"),
       slab = t2$display_name, 
       main = fits[[1]]$formula[[3]], cex.main = 0.8)

hh <- par("usr")
text(hh[1], max(rowvec) + 2, "Metabolite", pos = 4, cex = 0.8)
text(hh[2], max(rowvec) + 2, "OR [95% CI]", pos = 2, cex = 0.8)

# -------------------------------------------------------------------

fits <-  apply(ints, 2, function(x) clogit(CT ~ BMI + SMK + DIABETE + RTH + ALCOHOL + 
                                             SAMPYEAR +
                                             DURTHSDIAG + 
                                             CENTTIME +  STOCKTIME + # only have effect on Lactate
                                             #RACK + 
                                             strata(MATCH) + x, data = dat, subset = MENOPAUSE == 0))

t2  <- map_df(fits, tidy) %>% filter(str_detect(term, "x")) %>% bind_cols(cmpd_meta) %>% arrange(description)

# Plot data with Metafor
# Get vectors for row spacings using groups (may add compound classes later)
cmpds_ordered <- cmpd_meta %>% arrange(description) %>% mutate(row = 1:n() + (as.numeric(description) - 1))
rowvec <- cmpds_ordered$row

par(mar=c(5,4,2,2))
library(metafor)
forest(t2$estimate, ci.lb = t2$conf.low, ci.ub = t2$conf.high, refline = 1,
       ylim = c(1, max(rowvec) + 3), #at = 0:5,
       xlab = "OR per SD increase in concentration", 
       transf = exp, rows = rowvec, efac = 0.5,
       pch = 18, cex = 0.8, psize = 1.5, 
       annosym = c("  (", " to ", ")"),
       slab = t2$display_name, 
       main = fits[[1]]$formula[[3]], cex.main = 0.8)

hh <- par("usr")
text(hh[1], max(rowvec) + 2, "Metabolite", pos = 4, cex = 0.8)
text(hh[2], max(rowvec) + 2, "OR [95% CI]", pos = 2, cex = 0.8)
