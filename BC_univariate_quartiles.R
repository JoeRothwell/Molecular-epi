# BC risk models for metabolites
library(tidyverse)
library(readxl)

# Update 11/10/19: Updated with unscaled data (see BC_compounds_scaling for old data)
ints0 <- read_tsv("1510_XMetaboliteE3N_cpmg_unscaled.txt")

# Look at intensity ranges for each compound
walk2(ints0, colnames(ints0), ~ plot(.x, main = .y, col = ifelse(.x < 0, "red", "grey")))
which(apply(as.matrix(ints0), 2, min) < 0)
# 3 compounds have values < 0: formate, hypoxanthine, inosine

# Lifestyle data. Subset variables needed
meta <- read_csv("Lifepath_meta.csv", na = "9999") %>%
  select(CT, BMI, SMK, DIABETE, RTH, ALCOHOL, DURTHSDIAG, CENTTIME, STOCKTIME, RACK, MATCH, MENOPAUSE) %>%
  mutate_at(vars(SMK, DIABETE, RACK, MATCH), as.factor)

# Analysis by quartiles of metabolite concentration (resist outliers). Get quartiles for each compound from unscaled data
quartiles <- ints0 %>% mutate_all(funs(cut_number(., n = 4, labels = 1:4))) 
dat <- cbind(meta, quartiles)

library(survival)

# CLR models to get odds ratios for metabolites
# All subjects
fits0 <- apply(quartiles, 2, function(x) {
  Q1Q4 <- x == 1 | x == 4
  clogit(CT ~  x[Q1Q4] + BMI + SMK + DIABETE + RTH + ALCOHOL + DURTHSDIAG + CENTTIME + 
           STOCKTIME + strata(MATCH), data = meta[Q1Q4, ])
} )

library(lme4)
fitsX <- apply(quartiles, 2, function(x) {
  Q1Q4 <- x == 1 | x == 4
  glmer(CT ~  x[Q1Q4] + BMI + SMK + DIABETE + RTH + ALCOHOL + DURTHSDIAG + CENTTIME + 
          STOCKTIME + (1|MATCH), family = "binomial", data = meta[Q1Q4, ])
} )


# Subgroup analysis

# Pre-menopausal only
dat1 <- dat[dat$MENOPAUSE == 0, ]



fits1 <- apply(quartiles, 2, function(x) {
  Q1Q4 <- x == 1 | x == 4
  clogit(CT ~  x[Q1Q4] + BMI + SMK + DIABETE + RTH + ALCOHOL + DURTHSDIAG + CENTTIME + 
           STOCKTIME + RACK + strata(MATCH), data = meta[Q1Q4, ])
} )

# Post-menopausal only
dat2 <- dat[dat$MENOPAUSE == 1, ]

fits2 <- apply(quartiles, 2, function(x) {
  Q1Q4 <- x == 1 | x == 4
  clogit(CT ~  x[Q1Q4] + BMI + SMK + DIABETE + RTH + ALCOHOL + DURTHSDIAG + CENTTIME + 
           STOCKTIME + strata(MATCH), data = meta[Q1Q4, ])
} )

cmpd_meta <- read.csv("NMR_cmpd_metadata_new.csv")

library(broom)
t2 <- map_df(fits0, tidy) %>% filter(str_detect(term, "x")) %>% bind_cols(cmpd_meta) %>% arrange(description)
t2a <- map_df(fits1b, tidy) %>% filter(term == "x") %>% bind_cols(cmpd_meta) %>% arrange(description)

# Plot data with Metafor (pre-menopausal for manuscript)
# Get vectors for row spacings using groups (may add compound classes later)
cmpds_ordered <- cmpd_meta %>% arrange(description) %>% mutate(row = 1:n() + (as.numeric(description)-1))
rowvec <- cmpds_ordered$row

par(mar=c(5,4,1,2))
library(metafor)
forest(t2$estimate, ci.lb = t2$conf.low, ci.ub = t2$conf.high, refline = 1,
       ylim = c(1, max(rowvec) + 3), #at = 0:5,
       xlab = "Multivariable-adjusted odds ratio", 
       transf = exp, rows = rowvec, efac = 0.5,
       pch = 18, cex = 0.8, psize = 1.5, 
       annosym = c("  (", " to ",")"),
       slab = t2$display_name)

hh <- par("usr")
text(hh[1], max(rowvec) + 2, "Metabolite", pos = 4, cex = 0.8)
text(hh[2], max(rowvec) + 2, "OR [95% CI]", pos = 2, cex = 0.8)
