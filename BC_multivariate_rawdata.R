# Breast cancer study multivariate models from unbinned NMR data

# Data from Elodie Jobard 27-6-2019
# Read in data with outliers and QCs (n=1739) or with QCs only


# Subset samples only and plot PCA with pareto scaling
explore.data <- function(dataset = "QCs", data.only = F) {
  
  library(tidyverse)
  library(readxl)
  library(MetabolAnalyze)
  
  # Read in NMR data and metadata if necessary 
  
  if(dataset == "All") {
    mat <- read_delim("X_AlignedE3NData_cpmg_ssCitPEG_1112.txt", delim = ";", n_max = 1738) 
  } else {
    
    dat <- read_tsv("1510_XAlignedE3NcpmgssCitPEGfinal.txt") %>% select(-8501)
    metaQC <- read_xlsx("1510_MatriceY_CohorteE3N_Appar.xlsx", sheet = 4, na = ".")
    
    # Subset samples only
    samp <- metaQC$TYPE_ECH == 1
    mat <- dat[samp, ] %>% as.matrix
    
    # Subset samples only from metadata   
    meta <- metaQC[samp, ]
  }
  
  
  anyNA(mat)
  sum(apply(mat, 2, anyNA))
  
  # remove zero variance columns
  zerovar <- sum(apply(mat, 2, var) == 0)
  print(paste("There are", zerovar, "zero variance columns"))
  
  # There are 1116. Place in logical vector
  nonzerovar <- apply(mat, 2, var) != 0
  which(nonzerovar)
  
  mat0 <- mat[ , nonzerovar]
  print(paste("Dimensions", dim(mat0)))

  # Scale and run PCA  
  scalemat <- scaling(mat0, type = "pareto")
  
  if(data.only == T) return(scalemat)
  
  pca <- prcomp(scalemat, scale. = F, center = T, rank. = 10)
  
  output <- data.frame(pca$x) %>% bind_cols(meta)
}

scores <- explore.data()
#scores1 <- explore.data(dataset = "All")

# Plot data
library(ggplot2)
ggplot(scores, aes(PC1, PC2)) + geom_text(aes(label = WEEKS)) + theme_bw()

library(pca3d)
pca2d(as.matrix(scores[, 1:10]), group = scores$WEEKS)
box(which = "plot", lty = "solid")


#--------------------------------------------------------------------------------------------------


library(tidyverse)
library(readxl)

# 1694 obs. of 8501 NMR variables (outliers removed, one NA variable in 8501st col)
raw <- read_tsv("1510_XAlignedE3NcpmgssCitPEGfinal.txt") %>% select(-8501)

# 1739 obs. of 8500 NMR variables (all samples)
raw1 <- read_delim("X_AlignedE3NData_cpmg_ssCitPEG_1112.txt", delim = ";")

# Metadata are stored in the following Excel file. The third sheet has the metadata w/o outliers
# (hope it's in the same order as the NMR data!)
metaQC <- read_xlsx("1510_MatriceY_CohorteE3N_Appar.xlsx", sheet = 4, na = ".")


# Remove QCs from metadata and select columns
meta <- metaQC[samp, ] %>% 
  select(ID, RACK, WEEKS, CT, MATCH, PLACE, CENTTIME, CENTTIMECat1, FASTING, SMK, BMI, BMICat1, SAMPYEAR,
        AGE, BP, RTH, ALCOHOL, Trait_Horm, MENOPAUSE, DIABETE, STOCKTIME)


# pca0 <- prcomp(raw, scale. = F, center = T)
pca <- prcomp(samples0, scale. = T, center = T)
pca <- prcomp(raw1, scale. = T, center = T)

library(pca3d)
pca2d(pca, group = as.factor(meta$WEEKS), show.labels = T)
box(which = "plot", lty = "solid")

scores <- data.frame(pca$x) %>% bind_cols(meta)

library(ggplot2)
ggplot(scores, aes(PC1, PC2)) + geom_text(aes(label = MATCH, col = "blue"))

# Remove QCs from metadata

library(pcpr2)

Y_meta <- meta %>% select(WEEKS, PLACE, AGE, BMI, MENOPAUSE, FASTING, SMK, DIABETE, CENTTIMECat1, SAMPYEAR, STOCKTIME) %>%
  mutate_at(vars(-AGE, -BMI, -WEEKS, -STOCKTIME), as.factor)

Y_meta <- meta %>% select(BMI, MENOPAUSE, FASTING, SMK, DIABETE) %>% mutate_at(vars(-BMI), as.factor)
output <- runPCPR2(samples, Y_meta)
plotProp(output)
