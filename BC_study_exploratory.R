# Breast cancer NMR metabolomics data exploratory analysis
# Scaled and unscaled data (already done in SIMCA). Updated after receiving scaled data.
library(tidyverse)
library(readxl)

# Read scaled and unscaled data (for comparison) and metadata
ints0 <- read_tsv("1510_XMetaboliteE3N_cpmg_unscaled.txt")
ints <- read_tsv("1507_XMetabolite_std_cpmg_E3N.txt")
meta <- read.csv("Lifepath_meta.csv")

# Plot metabolites one by one (need walk2 because names are an attribute of df)
walk2(ints, colnames(ints), ~ plot(.x, main = .y, col = meta$RACK))

# subset IDs to get subjects included in CC. Get positions of final CC samples in metadata
samples <- ints$CODBMB %in% meta$CODBMB

# For baseline characteristics table: see BC_baseline_char.R

# Distributions all data
par(mfrow = c(2,1))
hist(as.matrix(ints0[, -1]), breaks = 50, col = "dodgerblue")
hist(as.matrix(ints[ , -1 ]), breaks = 50, col = "dodgerblue")

# Median intensities before scaling
dev.off()
points0 <- apply(ints0, 2, median)
plot(points0, col = "white", main = "Median compound intensities")
text(points0, labels = colnames(ints0))

# After scaling
points <- apply(ints[, -1], 2, median)
plot(points, col = "white", main = "Median compound intensities")
text(points, labels = colnames(ints[, -1]))

# Time to centrifugation vs fasting status
boxplot(CENTTIME ~ FASTING, data = meta, varwidth = T)
meta1 <- meta[meta$CENTTIME < 100, ]
library(gplots)
boxplot2(CENTTIME ~ FASTING, data = meta1, varwidth = T, col = "dodgerblue",
         xlab = "Fasting status", ylab = "TBC")

# Check correlations
library(corrplot)
cormat <- cor(ints[, -1])
colnames(cormat) <- NULL
corrplot(cormat, method = "square", tl.col = "black", tl.cex = 0.8)

# The three fatty acids are highly correlated, valine and leucine, NAC1 and 2
# The fatty acids are inversely correlated with many compounds

# Run PCA of all samples
pca0 <- prcomp(ints0, scale. = F, center = F)
biplot(pca0)

pca <- prcomp(ints[, -1], scale.=F)
biplot(pca)

# Plot 
library(pca3d)
par(mfrow = c(1, 2))
pca2d(pca)
title("A", font.main = 1, adj = 0)
box(which = "plot", lty = "solid")

# Run PCPR2 to compare sources of variability
library(pcpr2)
alldata <- inner_join(meta, ints, by = "CODBMB")
X_DataMatrixScaled <- select(alldata, `3Hydroxybutyrate`:Succinate) %>% as.matrix
Z_Meta <- alldata %>%
  select(WEEKS, PLACE, AGE, BMI, MENOPAUSE, FASTING, SMK, DIABETE, CENTTIMECat1, SAMPYEAR, STOCKTIME) %>%
  mutate_at(vars(-AGE, -BMI, -WEEKS, -STOCKTIME), as.factor)

par(mar=c(6,5,4,2))
props <- runPCPR2(X_DataMatrixScaled, Z_Meta)
plotProp(props, main = "B", font.main = 1, adj = 0)

# PCA of subject IDs only
ints1 <- alldata %>% select(`3Hydroxybutyrate`:Succinate)
pca1 <- prcomp(ints1, scale.=F)

par(mfrow = c(1,2))
plt <- pca2d(pca1, group = alldata$CT)
title("Metabolite profiles of 1582 samples", font.main = 1)
box(which = "plot", lty = "solid")
legend("topleft", legend = plt$groups, col=plt$colors, pch=plt$pch)

# Adjust using residuals method
adj <- function(x) residuals(lm(x ~ BMI + SMK + DIABETE, data = alldata))
adjmat <- apply(ints1, 2, adj)

# Repeat PCA
pca2 <- prcomp(adjmat, scale.=F)
pca2d(pca2, group = alldata$CT)
title("Residuals-adjusted metabolite\nprofiles of 1582 samples", font.main = 1)
box(which = "plot", lty = "solid")
legend("topleft", legend = plt$groups, col=plt$colors, pch=plt$pch)



# Other files ----

# Rawest feature data
raw <- read_delim("D:/J_ROTHWELL/X_AlignedCohorteE3NData_cpmg_ssCitPEG_0612.txt", delim = ";")

# Excel files with intermediate steps
# List of 54 compounds and IDs
#xl1 <- read_xlsx("C:/J_ROTHWELL/1505_E3N_Identification.xlsx")

# Previous version of intensities above without scaling
#dat2 <- read_tsv("C:/J_ROTHWELL/1507_XMetaboliteE3N_cpmg.txt")

# List of metabolites and corresponding clusters
#xl2 <- read_xlsx("C:/J_ROTHWELL/1507_ClusterAnnotation_E3N.xlsx")

# Metadata and intensities
#xl4 <- read_xlsx("C:/J_ROTHWELL/1603_MatriceY_CohorteE3N_Appar.xlsx", sheet = 6)

#ints <- read_tsv("C:/J_ROTHWELL/1507_XMetabolite_std_cpmg_E3N.txt")
#ints <- ints[, -1]
#ints1 <- read_tsv("C:/J_ROTHWELL/1507_XMetaboliteE3N_cpmg.txt", skip = 1)

# Food intake data
#path <- "Y:/RepPerso/Fabienne WILM/02_Demandes_ponctuelles/10_LIFEPATH/TABLES"
#list.files(path)

# Food intake and other data
#meta1 <- read_csv("D01_20171031_LIFEPATH.csv")
#meta2 <- read_csv("D01_20161018_LIFEPATH.csv")
#meta3 <- read_csv("D01_20150917_LIFEPATH.csv")



