# BC risk models for metabolites
library(tidyverse)
library(readxl)

# Update 11/10/19: Updated with unscaled data (see BC_compounds_scaling for old data)
ints0 <- read_tsv("1510_XMetaboliteE3N_cpmg_unscaled.txt")
cmpd_meta <- read.csv("NMR_cmpd_metadata_new.csv")

# Lifestyle data. Subset variables needed
meta <- read_csv("Lifepath_meta.csv", na = "9999") %>%
  select(CT, BMI, SMK, DIABETE, RTH, ALCOHOL, DURTHSDIAG, CENTTIME, STOCKTIME, RACK, MATCH, MENOPAUSE, SAMPYEAR) %>%
  mutate_at(vars(SMK, DIABETE, RACK, MATCH, SAMPYEAR), as.factor)

# Remove problem racks
#meta1 <- meta %>% filter(!RACK %in% c(10, 29, 33, 34))
#filt <- !(meta$RACK %in% c(10, 29, 33, 34))

# Analysis by quartiles of metabolite concentration (resist outliers). Get quartiles for each compound from unscaled data
quartiles <- ints0 %>% mutate_all(funs(cut_number(., n = 4, labels = 1:4))) 

library(survival)

# All subjects, pre or post menopausal

# Apply across all compounds
fits0 <- apply(quartiles, 2, function(x) {
  Q1Q4 <- x == 1 | x == 4
  clogit(CT ~  x[Q1Q4] + BMI + SMK + DIABETE + RTH + ALCOHOL + 
                         DURTHSDIAG +
                         CENTTIME + STOCKTIME + #RACK +
                         strata(MATCH), data = meta[Q1Q4, ])
} )

library(broom)
t2 <- map_df(fits0, tidy) %>% filter(str_detect(term, "x")) %>% bind_cols(cmpd_meta) %>% arrange(description)

# Plot data with Metafor
# Get vectors for row spacings using groups (may add compound classes later)
cmpds_ordered <- cmpd_meta %>% arrange(description) %>% mutate(row = 1:n() + (as.numeric(description) - 1))
rowvec <- cmpds_ordered$row

par(mar=c(5,4,2,2))
library(metafor)
forest(t2$estimate, ci.lb = t2$conf.low, ci.ub = t2$conf.high, refline = 1,
       ylim = c(1, max(rowvec) + 3), #at = 0:5,
       xlab = "OR 4th vs 1st quartile", 
       transf = exp, rows = rowvec, efac = 0.5,
       pch = 18, cex = 0.8, psize = 1.5, 
       annosym = c("  (", " to ", ")"),
       slab = t2$display_name, 
       main = fits0[[1]]$formula[[3]], cex.main = 0.8)

hh <- par("usr")
text(hh[1], max(rowvec) + 2, "Metabolite", pos = 4, cex = 0.8)
text(hh[2], max(rowvec) + 2, "OR [95% CI]", pos = 2, cex = 0.8)

# ------------------------------------------

# Post-menopausal only

post <- meta$MENOPAUSE == 1
meta1 <- meta[post, ]
ints.post <- ints0[post, ]
quartiles1 <- ints.post %>% mutate_all(funs(cut_number(., n = 4, labels = 1:4))) 
dat1 <- cbind(meta1, quartiles1)

# Apply across all compounds
fits0 <- apply(quartiles1, 2, function(x) {
  Q1Q4 <- x == 1 | x == 4
  clogit(CT ~  x[Q1Q4] + BMI + SMK + DIABETE + RTH + ALCOHOL + 
                         DURTHSDIAG + #SAMPYEAR +
                         CENTTIME + STOCKTIME + #RACK +
                         strata(MATCH), data = meta1[Q1Q4, ])
} )

library(broom)
t3 <- map_df(fits0, tidy) %>% filter(str_detect(term, "x")) %>% bind_cols(cmpd_meta) %>% arrange(description)

cmpds_ordered <- cmpd_meta %>% arrange(description) %>% mutate(row = 1:n() + (as.numeric(description) - 1))
rowvec <- cmpds_ordered$row

par(mar=c(5,4,2,2))
library(metafor)
forest(t3$estimate, ci.lb = t3$conf.low, ci.ub = t3$conf.high, refline = 1,
       ylim = c(1, max(rowvec) + 3), #at = 0:5,
       xlab = "OR 4th vs 1st quartile", 
       transf = exp, rows = rowvec, efac = 0.5,
       pch = 18, cex = 0.8, psize = 1.5, 
       annosym = c("  (", " to ", ")"),
       slab = t3$display_name, 
       main = fits0[[1]]$formula[[3]], cex.main = 0.8)

hh <- par("usr")
text(hh[1], max(rowvec) + 2, "Metabolite", pos = 4, cex = 0.8)
text(hh[2], max(rowvec) + 2, "OR [95% CI]", pos = 2, cex = 0.8)

# ------------------------------------------

# Pre-menopausal only

pre <- meta$MENOPAUSE == 0
meta1 <- meta[pre, ]
ints.pre <- ints0[pre, ]
quartiles1 <- ints.pre %>% mutate_all(funs(cut_number(., n = 4, labels = 1:4))) 
dat1 <- cbind(meta1, quartiles1)

# Apply across all compounds
fits0 <- apply(quartiles1, 2, function(x) {
  Q1Q4 <- x == 1 | x == 4
  clogit(CT ~  x[Q1Q4] + BMI + SMK + DIABETE + RTH + ALCOHOL + 
           DURTHSDIAG + #SAMPYEAR +
           CENTTIME + STOCKTIME + #RACK +
           strata(MATCH), data = meta1[Q1Q4, ])
} )

library(broom)
t3 <- map_df(fits0, tidy) %>% filter(str_detect(term, "x")) %>% bind_cols(cmpd_meta) %>% arrange(description)

cmpds_ordered <- cmpd_meta %>% arrange(description) %>% mutate(row = 1:n() + (as.numeric(description) - 1))
rowvec <- cmpds_ordered$row

par(mar=c(5,4,2,2))
library(metafor)
forest(t3$estimate, ci.lb = t3$conf.low, ci.ub = t3$conf.high, refline = 1,
       ylim = c(1, max(rowvec) + 3), #at = 0:5,
       xlab = "OR 4th vs 1st quartile", 
       transf = exp, rows = rowvec, efac = 0.5,
       pch = 18, cex = 0.8, psize = 1.5, 
       annosym = c("  (", " to ", ")"),
       slab = t3$display_name, 
       main = fits0[[1]]$formula[[3]], cex.main = 0.8)

hh <- par("usr")
text(hh[1], max(rowvec) + 2, "Metabolite", pos = 4, cex = 0.8)
text(hh[2], max(rowvec) + 2, "OR [95% CI]", pos = 2, cex = 0.8)
