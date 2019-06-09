# Adjustment of biocrates controls to remove confounding
source("CRC_data_prep.R")

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
