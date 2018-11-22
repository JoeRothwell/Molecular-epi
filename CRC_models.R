# Models for colorectal cancer status and WCRF score. First source CRC_data_prep.R which preps the two CCs.
source(CRC_data_prep.R)
source(Metabolic_signature_WCRF.R)

library(haven)
library(tidyverse)
meta <- read_dta("D:/clrt_caco.dta")

# Convert categorical variables to factors
var.list  <- c("Country", "Center", "Sex", "Match_Caseset", "Smoke_Stat", "L_School")
meta      <- meta %>% mutate_at(vars(var.list), as.factor)

# Get names of individual contributor variables to WCRF score
wcrf_vars <- c("Wcrf_Bmi", "Wcrf_Pa", "Wcrf_Fwg_Cal", "Wcrf_Pf_Cal", "Wcrf_Fv_Cal", "Wcrf_Fibt_Cal",
               "Wcrf_Meat_Cal", "Wcrf_C_Cal")

# Get variable descriptions (for forest plot)
scorecomp <- c("1. Body fatness", "2. Physical activity", "3. Energy density/sugary drinks", "4. FV intake", 
               "5. Foods of plant origin", "6. Fibre intake", "7. Meat intake", 
               "       Overall WCRF score (cal.)", "       Signature metabolites")

library(survival)
library(broom)
# Function to run model case-control status from main WCRF score components and get tidy output
clr_crc <- function(dat, sig = "none") {
  
  # model with covariates
  clr_fit <- function(x) clogit(Cncr_Caco_Clrt ~ x + Qe_Energy + L_School + Smoke_Stat + strata(Match_Caseset), data = dat)
  mat <- dat %>% select(wcrf_vars)
  
  # bind predicted scores from Metabolic_signature_WCRF.R
    mat <- if(sig == "large") mat %>% bind_cols(large.nofast) else if(sig == "small") mat %>% bind_cols(small) else mat
  print(paste("dimensions:", dim(mat)))
  multifit <- apply(mat, 2, clr_fit)
  
  # subset risk estimate for score only
  output <- map_df(multifit, tidy) %>% filter(term == "x")
}
t1 <- clr_crc(meta)
t1 <- clr_crc(crc1, sig = "small")
t1 <- clr_crc(crc2, sig = "large")

# Plot output of each study as forest
library(metafor)
par(mar=c(5,4,1,2))
forest(t1$estimate, ci.lb = t1$conf.low, ci.ub = t1$conf.high, refline = 1, 
       xlab = "Odds ratio (per unit increase in score)", pch = 18, transf = exp, psize = 1.5, #xlim = c(-1.8, 4),
       slab = scorecomp) 

hh <- par("usr")
text(hh[1], nrow(t1) + 2, "Component of score", pos = 4)
text(hh[2], nrow(t1) + 2, "OR [95% CI]", pos = 2)

# Plot only overall score and signature from each study
t1 <- clr_crc(meta) %>% slice(n())
t2 <- clr_crc(crc1, sig = "small") %>% slice(n())
t3 <- clr_crc(crc2, sig = "large") %>% slice(n())
t4 <- bind_rows(t1, t2, t3)

par(mar=c(5,4,1,2))
forest(t4$estimate, ci.lb = t4$conf.low, ci.ub = t4$conf.high, refline = 1, 
       xlab = "Odds ratio (per unit increase in score)", pch = 18, transf = exp, psize = 1.5, #xlim = c(-1.8, 4),
       slab = c("Overall score, full study", "Signature, small metabolomics", "Signature, large metabolomics"))

hh <- par("usr")
text(hh[1], nrow(t4) + 2, "Parameter", pos = 4)
text(hh[2], nrow(t4) + 2, "OR [95% CI]", pos = 2)





