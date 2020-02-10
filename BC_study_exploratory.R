# Breast cancer NMR metabolomics data exploratory analysis
# Scaled and unscaled data (already done in SIMCA). Updated after receiving scaled data.
library(tidyverse)
library(readxl)
library(gplots)

# Read unscaled data and metadata (each 1882 obs)
ints0 <- read_tsv("1510_XMetaboliteE3N_cpmg_unscaled.txt")
ints <- scale(ints0)
meta <- read.csv("Lifepath_meta.csv")

# Read scaled data and subset samples from meta
ints.all <- read.delim("1507_XMetabolite_std_cpmg_E3N.txt")

which(apply(as.matrix(ints0), 2, min) < 0)
# 3 compounds have values < 0: formate, hypoxanthine, inosine

# subset IDs to get subjects included in CC. Get positions of final CC samples in metadata
samples <- ints.all$CODBMB %in% meta$CODBMB
ints <- ints.all[samples, ]

# Distributions all data
par(mfrow = c(2,1))
hist(as.matrix(ints0), breaks = 50, col = "dodgerblue")
hist(as.matrix(ints),  breaks = 50, col = "dodgerblue")

# Plot metabolites individually
library(RColorBrewer)

#palette(rainbow(5))
palette(brewer.pal(5, "Dark2"))
plot.ts(ints[, 1:10], type = "p", col = meta$SAMPYEAR, main = "Hydroxybutyrate - Creatine")
plot.ts(ints[, 11:20], type = "p", col = meta$SAMPYEAR, main = "Creatinine - Glutamate")
plot.ts(ints[, 21:30], type = "p", col = meta$SAMPYEAR, main = "Glycerol - Isoleucine")
plot.ts(ints[, 31:40], type = "p", col = meta$SAMPYEAR, main = "Lactate - Malonate")
plot.ts(ints[, 41:ncol(ints)], type = "p", col = meta$SAMPYEAR, main = "Malonate - Succinate")

# Problem compounds only
plot.ts(ints[, c(13, 14, 23, 24, 39, 43)], type = "p", col = meta$SAMPYEAR)

# Or use walk2 to go through plotting all columns
par(mfrow = c(5, 1), mai = c(0.3, 0.5, 0.2, 0.1))
walk2(ints, colnames(ints), ~ plot(.x, main = .y, col = meta$SAMPYEAR))
walk2(ints, colnames(ints), ~ boxplot2(.x ~ meta$MENOPAUSE, main = .y, top = T))

# Plot distributions by compound and menopause
ints1 <- cbind(meno = meta$MENOPAUSE, ints) %>% data.frame
ints.melt <- gather(ints1, compound, value, -meno)
ggplot(ints.melt, aes(as.factor(meno), value)) + geom_boxplot(outlier.shape = NA) + theme_bw() +
  geom_jitter(width=0.3, alpha = 0.3, size = 0.1) + facet_wrap( ~ compound, ncol = 8, scales = "free_y")

# For baseline characteristics table: see BC_baseline_char.R

# Median intensities before scaling
dev.off()
Intensity <- apply(ints0, 2, median)
plot(Intensity, col = "white", main = "Median compound intensities")
text(Intensity, labels = colnames(ints0))
# Fatty acids and glucose strongest signals, then lactate, glycerophosphocholine.

# Visualise case-control differences
ggplot(ints, aes(x= as.factor(meta$CT), y=log(Hypoxanthine))) + 
  geom_line(aes(group = meta$MATCH), alpha = 0.5) + 
  geom_point(aes(group = meta$MATCH), alpha = 0.5) + 
  theme_minimal() + facet_grid( ~ meta$MENOPAUSE)

# Check correlations: fatty acids are highly correlated
library(corrplot)
#cormat <- cor(ints, use = "pairwise.complete.obs")
cormat <- cor(ints[, -1])
colnames(cormat) <- NULL
corrplot(cormat, method = "square", tl.col = "black", tl.cex = 0.8, order = "hclust")

# Dendrogram
library(dendextend)
dend <- hclust(dist(cormat)) %>% as.dendrogram
par(mar=c(1,1,1,8))
dend %>% 
  set("labels_col", value = c("skyblue", "orange", "grey"), k = 3) %>%
  set("branches_k_color", value = c("skyblue", "orange", "grey"), k = 3) %>%
  plot(horiz = T, axes = F)


# Run PCA of all samples
pca0 <- prcomp(ints0, scale. = F, center = F)
biplot(pca0)

pca <- prcomp(ints[, -1], scale.=F)
biplot(pca, cex = 0.7)
scores <- data.frame(pca$x)
ggplot(scores, aes(x = PC1, y = PC2, col = as.factor(meta$MENOPAUSE))) + 
  geom_point() + scale_colour_manual(values = c("black", "grey"))



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
pca2 <- prcomp(adjmat, scale. = F)
pca2d(pca2, group = alldata$CT)
title("Residuals-adjusted metabolite\nprofiles of 1582 samples", font.main = 1)
box(which = "plot", lty = "solid")
legend("topleft", legend = plt$groups, col=plt$colors, pch=plt$pch)

# Test datasets for scaling -----------------------

library(tidyverse)
library(readxl)

# Original data, 1623 observations of 44 compounds - appears to be Pareto scaled
ints <- read_tsv("1507_XMetabolite_std_cpmg_E3N.txt")

# Update 11/10/19: Updated with unscaled data (will scale to unit variance). Added "_unscaled" to filename.
# Note: NAC1 and NAC2 are merged and only 1582 observations
ints0 <- read_tsv("1510_XMetaboliteE3N_cpmg_unscaled.txt") #%>% as.matrix
ints1 <- scale(ints0)

# Lifestyle data. Subset variables needed
meta <- read_csv("Lifepath_meta.csv", na = "9999")
samples <- ints$CODBMB %in% meta$CODBMB
ints <- ints[samples, -1]

# Exploratory of different datasets. Find differences between received datasets (formerly BC_compounds_scaling)
library(MetabolAnalyze)

# 1. Original data (1507); 2. New data (no scaling); 3. New data (unit scaled); New data (Pareto scaled)
new     <- ints0 %>% as.matrix
unit    <- scale(ints0)
pareto1 <- scaling(ints0, type = "pareto")

ll <- lapply(list(ints, new, unit, pareto1), function(x) prcomp(x, scale. = F, center = F))
par(mfrow = c(2,2))
lapply(ll, function(x,y) { pca2d(x)
  title(main = paste("Scores plot"))
  box(which = "plot")
})

#how to make hotpink, limegreen, orange, dodgerblue and give different titles?

# Conclusion: the old dataset seems to be a unit variance scaled version of the new dataset
# (although they are not exactly the same; the new datasets has merged NAC1 and NAC2)
# Decision: use unscaled dataset and apply unit scaling, discard old scaled dataset

# Old ------------------------------

# Time to centrifugation vs fasting status
boxplot(CENTTIME ~ FASTING, data = meta, varwidth = T)
meta1 <- meta[meta$CENTTIME < 100, ]
library(gplots)
boxplot2(CENTTIME ~ FASTING, data = meta1, varwidth = T, col = "dodgerblue",
         xlab = "Fasting status", ylab = "TBC")

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



