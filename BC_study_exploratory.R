# LIFEPATH data exploratory
library(tidyverse)
library(readxl)

# Read 1623 observations of 44 intensity variables (appears to be final scaled data) and metadata
ints <- read_tsv("1507_XMetabolite_std_cpmg_E3N.txt")
meta <- read.csv("Lifepath_meta.csv")

# subset IDs to get subjects included in CC. Get positions of final CC samples in metadata
samples <- ints$CODBMB %in% meta$CODBMB

# baseline characteristics table: see baseline_lifepath.R

# Exploratory analysis --------------------------------------------
# Check total intensities for each metabolite
ints <- ints[samples, ]
plot(colSums(ints[ , -1]), xlab = "Compound number", ylab = "Scaled intensity",
     pch = 19, col = "dodgerblue", font.main = 1,
     main = "Summed intensities of 44 metabolites for 1582 subjects")

# Check correlations
library(corrplot)
cormat <- cor(ints[, -1])
colnames(cormat) <- NULL
corrplot(cormat, method = "square", tl.col = "black", tl.cex = 0.8)

# The three fatty acids are highly correlated, valine and leucine, NAC1 and 2
# The fatty acids are inversely correlated with many compounds

# Run PCA of all samples
pca <- prcomp(ints[, -1], scale.=F)
plot(pca)

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
#raw <- read_delim("C:/J_ROTHWELL/X_AlignedCohorteE3NData_cpmg_ssCitPEG_0612.txt", delim = ";")

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



