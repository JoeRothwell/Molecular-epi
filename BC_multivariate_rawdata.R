# Breast cancer study multivariate models from unbinned NMR data

# Data from Elodie Jobard 27-6-2019

library(tidyverse)
library(readxl)

# 1694 obs. of 8501 NMR variables (outliers removed, one NA variable in 8501st col)
raw <- read_tsv("1510_XAlignedE3NcpmgssCitPEGfinal.txt")

# 1739 obs. of 8500 NMR variables (all samples)
# raw1 <- read_delim("X_AlignedE3NData_cpmg_ssCitPEG_1112.txt", delim = ";")

# Metadata are stored in the following Excel file. The third sheet has the metadata w/o outliers
# (hope it's in the same order as the NMR data!)
metaQC <- read_xlsx("1510_MatriceY_CohorteE3N_Appar.xlsx", sheet = 4)

# Remove QCs from metadata and select columns
meta <- metaQC[samp, ] %>% select(ID, RACK, WEEKS, CT, MATCH, PLACE, CENTTIME, CENTTIMECat1, FASTING, SMK, BMI, BMICat1, SAMPYEAR,
                                  AGE, BP, RTH, ALCOHOL, Trait_Horm, MENOPAUSE, DIABETE, STOCKTIME)

#remove QCs and erroneous column from raw data
samp <- metaQC$TYPE_ECH == 1
samples <- raw[samp, -8501] %>% as.matrix

# pca0 <- prcomp(raw, scale. = F, center = T)
pca <- prcomp(samples, scale. = F, center = T)
library(pca3d)
pca2d(pca, group = as.factor(meta$PLACE), show.labels = T)
box(which = "plot", lty = "solid")

scores <- data.frame(pca$x) %>% bind_cols(meta)

library(ggplot2)
ggplot(scores, aes(PC1, PC2)) + #geom_point() + 
  theme_bw() +
  geom_text(aes(label = ID))

# Remove QCs from metadata


library(pcpr2)

Y_meta <- meta %>% select(WEEKS, PLACE, AGE, BMI, MENOPAUSE, FASTING, SMK, DIABETE, CENTTIMECat1, SAMPYEAR, STOCKTIME) %>%
  mutate_at(vars(-AGE, -BMI, -WEEKS, -STOCKTIME), as.factor)

Y_meta <- meta %>% select(BMI, MENOPAUSE, FASTING, SMK, DIABETE) %>% mutate_at(vars(-BMI), as.factor)
output <- runPCPR2(samples, Y_meta)
plotProp(output)
