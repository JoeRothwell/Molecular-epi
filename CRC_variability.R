# Find main sources of variability in EPIC controls. 
# Biocrates: First dataset, 3771 obs; updated November 2018 7191 obs
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
plot(obj, main = "Raw metabolite matrix")

# Biocrates
library(lme4)
adj   <- function(x) residuals(lmer(x ~ Center + batch_no + Sex + (1|Study), data = meta))
adjmat <- apply(logconcs, 2, adj)

obj2 <- runPCPR2(adjmat, meta)
plot(obj2, main = "After transformation")

# -------------------------------------

# Fatty acids
library(tidyverse)
# Exclude compounds with many missings. New version from Carine received 18/11/2018 with technical covariates
fa.ctrl <- readRDS("FA_WCRF_scores1.rds") %>% filter(STUDY != "Colorectum") %>% filter(Country != 6)
fa.ctrl$N_Serie <- as.numeric(fa.ctrl$N_Serie)

epic.vars <- read.csv("full_epic_vars.csv")
fa.ctrl <- left_join(fa.ctrl, epic.vars, by = "Idepic")

concs <- fa.ctrl %>% select(P14_0 : PCLA_9t_11c, -P24_0, -P20_0) %>% as.matrix
concs[concs == 0] <- NA
library(zoo)
concs1 <- na.aggregate(concs, FUN = function(x) min(x)/2)

# categorical variables to factors
meta <- fa.ctrl %>% select(Center, LABO, Sex, Bmi_C, STUDY)
var.list <- c("Center", "STUDY", "LABO")
meta <- meta %>% mutate_at(vars(var.list), as.factor)

# Run pcpr2
output1 <- runPCPR2(concs1, meta)
par(mfrow = c(1,2))
plot(output1, main = "Raw metabolite matrix")

# Adjust matrix for study, centre, batch, sex for Biocrates, subset calibrated scores 
library(lme4)
adj  <- function(x) residuals(lmer(x ~ LABO + STUDY + (1|Center), data = meta))
adjmat <- apply(concs1, 2, adj)

output2 <- runPCPR2(adjmat, meta)
plot(output2, main = "After transformation")

