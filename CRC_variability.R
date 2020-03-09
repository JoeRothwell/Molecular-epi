# EPIC controls. First dataset, 3771 obs; updated November 2018 7191 obs
# Rename factor levels, split Batch_MetBio into 2 cols, extract numeric variable,
# Remove 1694 CRC controls
source("CRC_prep_data.R")

# 1741 fasted controls without Greece
# Subset compounds (X data)
cmpds <- ctrl %>% 
 select(matches("Acylcarn_|Aminoacid_|Biogenic_|Glyceroph_|Sphingo_|Sugars_"), -starts_with("Outdq"))

#cmpds <- ctrlB

zerocols <- apply(cmpds, 2, function(x) sum(x, na.rm = T)) != 0
concs <- cmpds[, zerocols]

# Subset metadata (Z data)
meta <- ctrl %>% select(Center, batch_no, Sex, Bmi_C, Study)
meta$Center <- as.factor(meta$Center)
meta$Sex <- as.factor(meta$Sex)
meta$Study <- droplevels(meta$Study)

# Prepare controls matrix. Replace zero, impute with half mins, scale
concs[concs == 0] <- NA
library(zoo)
concs1 <- na.aggregate(concs, FUN = function(x) min(x)/2)
logconcs <- log2(concs1) %>% scale

# Run pcpr2
library(pcpr2)
obj <- runPCPR2(logconcs, meta)
par(mfrow = c(1,2))
plot(obj, col = "red", main = "A", adj = 0.1)

# Biocrates
library(lme4)
adj   <- function(x) residuals(lmer(x ~ Center + batch_no + Sex + (1|Study), data = meta))
adjmat <- apply(logconcs, 2, adj)

obj2 <- runPCPR2(adjmat, meta)
plot(obj2, col = "red")

# Fatty acids
# Adjust matrix for study, centre, batch, sex for Biocrates, subset calibrated scores   
concs <- ctrldat %>% select(one_of(common.cols))
adj   <- function(x) residuals(lmer(x ~ LABO + STUDY + (1|Center), data = ctrldat))



adjmat <- apply(logconcs, 2, adj)

