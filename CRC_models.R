# Models for colorectal cancer status and WCRF score. First source CRC_data_prep.R which preps the two CCs.
source("CRC_data_prep.R")
source(Metabolic_signature_WCRF.R)

library(haven)
library(tidyverse)
meta <- read_dta("clrt_caco.dta")

# Convert categorical variables to factors
var.list  <- c("Country", "Center", "Sex", "Match_Caseset", "Smoke_Stat", "L_School")
meta      <- meta %>% mutate_at(vars(var.list), as.factor)

# Get names of individual contributor variables to WCRF score
wcrf_vars <- c("Wcrf_Bmi", "Wcrf_Pa", "Wcrf_Fwg_Cal", "Wcrf_Pf_Cal", "Wcrf_Fv_Cal", "Wcrf_Fibt_Cal",
               "Wcrf_Meat_Cal", "Wcrf_C_Cal")

# Get variable descriptions (for forest plot)
scorecomp <- c("1. Body fatness", "2. Physical activity", "3. Energy density/sugary drinks", "4. FV intake", 
               "5. Foods of plant origin", "6. Fibre intake", "7. Meat intake", 
           "       Overall WCRF score (cal.)")
scorecomp2 <- c(scorecomp, "       Signature metabolites")
xtitle     <- "Odds ratio (per unit increase in score)"

library(survival)
library(broom)
# Function to run model case-control status from main WCRF score components and get tidy output
clr_crc <- function(dat, sig = "none") {
  
  # model with covariates
  mat <- dat %>% select(wcrf_vars)
  
  # bind predicted scores from Metabolic_signature_WCRF.R, print dimensions as check
  mat <- if(sig == "large") cbind(mat, large$score.2.comps) else if(sig == "small") cbind(mat, small$score.2.comps) else mat
  print(paste("dimensions:", dim(mat)))
  
  # Run model and apply across matrix
  fit <- function(x) clogit(Cncr_Caco_Clrt ~ x + Qe_Energy + L_School + Smoke_Stat + strata(Match_Caseset), data = dat)
  multifit <- apply(mat, 2, fit)
  
  # subset risk estimate for score only
  output <- map_df(multifit, tidy) %>% filter(term == "x")
}
t1 <- clr_crc(meta, sig = "none") %>% mutate(Study = scorecomp)
t2 <- clr_crc(crc1, sig = "small") %>% mutate(Study = scorecomp2)
t3 <- clr_crc(crc2, sig = "large") %>% mutate(Study = scorecomp2)

# Plot output of each study as forest
library(metafor)
par(mar=c(5,4,1,2))
# Full study ~ 14,000 subjects
forest(t1$estimate, ci.lb = t1$conf.low, ci.ub = t1$conf.high, refline = 1, xlab = xtitle, 
       transf = exp, pch = 18, psize = 1.5, slab = t1$Study) 
# Large study ~ 2370 subjects
forest(t2$estimate, ci.lb = t2$conf.low, ci.ub = t2$conf.high, refline = 1, xlab = xtitle, 
       transf = exp, pch = 18, psize = 1.5, slab = t2$Study) 
# SMall study  ~ 900 subjects
forest(t3$estimate, ci.lb = t3$conf.low, ci.ub = t3$conf.high, refline = 1, xlab = xtitle, 
       transf = exp, pch = 18, psize = 1.5, slab = t3$Study) 

hh <- par("usr")
text(hh[1], nrow(t2) + 2, "Component of score", pos = 4)
text(hh[2], nrow(t2) + 2, "OR [95% CI]", pos = 2)

# Plot only overall score and signature from each study
t4 <- bind_rows(t1[nrow(t1), ], t2[(nrow(t2)-1):nrow(t2), ], t3[(nrow(t3)-1):nrow(t3), ])

par(mar=c(5,4,1,2))
forest(t4$estimate, ci.lb = t4$conf.low, ci.ub = t4$conf.high, refline = 1, 
       xlab = "Odds ratio (per unit increase in score)", pch = 18, transf = exp, psize = 1.5,
       slab = c("WCRF score, full study", "WCRF score, small", "Metabolic signature, small", 
                "WCRF score, large", "Metabolic signature, large"))

hh <- par("usr")
text(hh[1], nrow(t4) + 2, "Parameter", pos = 4)
text(hh[2], nrow(t4) + 2, "OR [95% CI]", pos = 2)

# Fixed-effects meta-analysis of signatures
ma1 <- rma(estimate, sei = std.error, data=t4, method="FE", subset = 2:3)
# Random-effects meta-analysis of signatures
ma2 <- rma(estimate, sei = std.error, data=t4, method="REML", subset = 2:3)

# Of scores
ma3 <- rma(estimate, sei = std.error, data=t4, method="FE", subset = c(2,4))
ma4 <- rma(estimate, sei = std.error, data=t4, method="REML", subset = c(2,4))



par(mfrow = c(1,2))
slab <- c("WCRF small \nmetabolomics", "WCRF large \nmetabolomics")
forest(ma1, transf = exp, refline = 1, slab = slab, xlab = "OR")
forest(ma2, transf = exp, refline = 1, slab = slab, xlab = "OR")
forest(ma3, transf = exp, refline = 1, slab = slab, xlab = "OR")
forest(ma4, transf = exp, refline = 1, slab = slab, xlab = "OR")





