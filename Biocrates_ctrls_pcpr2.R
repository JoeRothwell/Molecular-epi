# Adjustment of biocrates controls to remove confounding
source("CRC_data_prep.R")

# Biocrates dataset
# Prepare controls matrix. Replace zero, impute with half mins, scale
concs <- as.matrix(controls)
concs[concs == 0] <- NA
concs1 <- na.aggregate(concs, FUN = function(x) min(x)/2)
logconcs <- log2(concs1) %>% scale

# PCPR2 of raw data
library(tidyverse)
Y_meta <- ctrl %>% select(Center, batch_no, Sex, Bmi_C, Study)

# Convert categoricals to factors
varlist <- c("Center", "Sex", "Study")
Y_meta <- Y_meta %>% mutate_at(vars(varlist), as.factor)

library(pcpr2)
props1 <- runPCPR2(logconcs, Y_meta)
plotProp(props1)

# adjust matrix for study, centre, sex, batch, BMI
adj <- function(x) residuals(lmer(x ~ Center + batch_no + 
                                    #Sex + Bmi_C + 
                                    (1|Study), data = Y_meta))
adjmat <- apply(logconcs, 2, adj)

props2 <- runPCPR2(adjmat, Y_meta)
plotProp(props2)

# Compare before and after adjustment
par(mfrow = c(1,2))
plotProp(props1, main = "A", adj = 0)
plotProp(props2)

#-------------------------------------------------------------------

# Fatty acids dataset
# Subset FA concentrations
library(tidyverse)
dat <- output[[1]]
common_cols <- output[[2]]
concs <- dat %>% select(one_of(common_cols))

# Processing steps
concs <- as.matrix(concs)
concs[concs == 0] <- NA
library(zoo)
concs <- na.aggregate(concs, FUN = function(x) min(x)/2)
logconcs <- log2(concs) %>% scale

#Y_meta <- dat %>% select(LABO, Country, STUDY)
#Y_meta <- dat %>% select(LABO, Country, STUDY, Center)

# Country and center are singular. Center explains more variability
Y_meta <- dat %>% select(LABO, Center, STUDY)

library(pcpr2)
# Warning: takes a long time due to many observations and factor levels
props1 <- runPCPR2(logconcs, Y_meta)
plotProp(props1)

library(lme4)
adj <- function(x) residuals(lmer(x ~ LABO + STUDY + (1|Center), data = Y_meta))
resmat <- apply(logconcs, 2, adj)
props2 <- runPCPR2(resmat, Y_meta)

par(mfrow = c(1,2))
plotProp(props1, main = "B", adj = 0)
plotProp(props2)

# Final adjustment to be used before PLS in Metabolic_signatures.R
