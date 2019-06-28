# Breast cancer study multivariate models from unbinned NMR data

# Data from Elodie Jobard 27-6-2019

library(tidyverse)
library(readxl)

# 1694 obs. of 8501 NMR variables (outliers removed, one NA variable in 8501st col)
raw <- read_tsv("1510_XAlignedE3NcpmgssCitPEGfinal.txt")

# 1739 obs. of 8500 NMR variables (all samples)
raw1 <- read_delim("X_AlignedE3NData_cpmg_ssCitPEG_1112.txt", delim = ";")

# Metadata are stored in the following Excel file. The third sheet has the metadata w/o outliers
# (hope it's in the same order as the NMR data!)
metaQC <- read_xlsx("1510_MatriceY_CohorteE3N_Appar.xlsx", sheet = 4)


#remove QCs and erroneous column from raw data
samp <- metaQC$TYPE_ECH == 1
samples <- raw[samp, -8501]

# Remove QCs from metadata
meta <- metaQC[samp, ] %>% select(WEEKS, CT, MATCH, PLACE, CENTTIME, CENTTIMECat1, FASTING, SMK, BMI, BMICat1, SAMPYEAR,
                                  AGE, BP, RTH, ALCOHOL, Trait_Horm, MENOPAUSE, DIABETE)

library(pcpr2)

