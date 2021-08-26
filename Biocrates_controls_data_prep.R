# Exploratory analysis of biocrates controls
library(haven)
ctrl <- read_dta("obes_metabo.dta")

# Outlier removal and adjustment for confounders
# Subset metabolite matrix, replace zero with NA
concs <- ctrl %>% select(Acylcarn_C0 : Sugars_H1) %>% as.matrix
concs[concs == 0] <- NA
  
# Impute missing concentrations with half minimum value (otherwise gives contrasts error for sex and study)
library(zoo)
concs1 <- na.aggregate(concs, FUN = function(x) min(x)/2)
logconcs <- log2(concs1)
pca <- prcomp(logconcs, scale. = T)
  
# Plot with group labels
library(pca3d)
pca2d(pcaconcs, group = ctrl$Study, show.group.labels = F, legend = "bottomright")
box(which = "plot", lty = "solid")
#PCA reveals two outliers on PC1, 3600 and 3722
outliers <- which(pcaconcs$x[, 1] > 40)

# Export data for PCPR2
  
#subset data and replot
ctrl1 <- ctrl[ -outliers, ]
ctrl1 <- ctrl %>% mutate(batch_no = as.numeric(flatten(str_extract_all(Batch_MetBio, "[0-9]+")))) %>% 
    slice(-outliers)
logconcs1 <- logconcs[ -outliers, ]
  #ctrl1 <- cbind(ctrl1, logconcs1)
return(list(ctrl1, logconcs1))

output <- prepdata()

# first is whole dataset outliers removed, second is matrix of log concs
ctrl1 <- output[[1]]
concs <- output[[2]]

# Adjust for different covariates by taking the residuals of a linear model
library(lme4)
adj1 <- function(x) residuals(lm(x ~ Sex, data = ctrl1))
adj2 <- function(x) residuals(lm(x ~ Sex + Study, data = ctrl1))
adj3 <- function(x) residuals(lm(x ~ Sex + Center, data = ctrl1))
adj4 <- function(x) residuals(lm(x ~ Sex + Study + Center, data = ctrl1))
adj5 <- function(x) residuals(lmer(x ~ Center + batch_no + Sex + (1|Study), data = ctrl1))

# Generate the matrices of residuals and rerun and plot PCA
pcares <- function(adjust = T, covar = "adj1", pca = F) {
  resmat <- if(adjust == F) concs else apply(concs, 2, covar)
  if(pca == T) {
    res.pca <- prcomp(resmat, scale. = T)
    pca2d(res.pca, group = ctrl1$Study, show.group.labels = F, legend = "bottomleft")
    box(which = "plot", lty = "solid")
  }
  return(resmat)
}

# Plot adjusted PCAs
resmat <- pcares(adjust = F)
# resmat <- pcares(covar = "adj4")
# resmat <- pcares(covar = "adj5")

# Extracts code names and label names and converts to a dataframe
# (previously Biocrates_cmpd_metadat.R script)

library(tidyverse)
library(haven)

# Biocrates  
ctrl <- read_dta("obes_metabo.dta") %>%
  select(matches("Acylcarn_|Aminoacid_|Biogenic_|Glyceroph_|Sphingo_|Sugars_"), 
         -starts_with("Outdq"))
labels <- sapply(controls, attr, "label")
tibble(name = names(labels), label = labels)


# Fatty acids
CRCfa1 <- read_dta("Database_Fatty acids.dta")

CRCfa <- CRCfa1 %>% select(P14_0 : PCLA_9t_11c) 
df <- tibble(Compound = colnames(CRCfa))