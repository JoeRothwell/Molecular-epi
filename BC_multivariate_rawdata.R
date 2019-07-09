# Breast cancer study multivariate models from unbinned NMR data (raw features)

# Data from Elodie Jobard 27-6-2019
# Read in data with outliers and QCs (n=1739) or with QCs only


# Subset samples only and plot PCA with pareto scaling
explore.data <- function(dataset = "QCs", data.only = F, exclude.tbc = F) {
  
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
    samp <- if(exclude.tbc == F) metaQC$TYPE_ECH == 1 else metaQC$TYPE_ECH == 1 & metaQC$CENTTIME < 100
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
  
  if(data.only == T) return(list(scalemat, meta))
  
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

# ----------------------------------------------------------------------------------------------------------

# Output Pareto-scaled data and perform PCPR2

dat1 <- explore.data(data.only = T, exclude.tbc = T)

concs <- dat[[1]]
concs1 <- dat1[[1]]

meta <- dat[[2]] %>%
  select(CT, MATCH, WEEKS, PLACE, AGE, BMI, MENOPAUSE, FASTING, SMK, DIABETE, CENTTIMECat1, CENTTIME, SAMPYEAR, STOCKTIME) %>%
  mutate_at(vars(-AGE, -BMI, -STOCKTIME, -CENTTIME), as.factor)

Z_Meta <- meta %>% select(-MATCH, -CT, -CENTTIME)

library(pcpr2)
props.raw <- runPCPR2(concs, Z_Meta)
plotProp(props.raw)

# Greatest sources of variability are BMI > DIABETE > PLACE
# Adjust for fixed effects only. FASTING covers CENTTIME
adj <- function(x) residuals(lm(x ~ PLACE + WEEKS + AGE + BMI + DIABETE + FASTING + SAMPYEAR, data = meta))


# Random effects model does not work, either boundary fit or does not converge
#library(lme4)
#adj <- function(x) residuals(lmer(x ~ BMI + DIABETE + 1|PLACE, data = meta))

# Samples were matched on age and menopausal status (at blood collection), collection centre, fasting status,
# blood collection data. Need to avoid including linearly dependent variables
adjmat <- apply(concs, 2, adj)

props.adj <- runPCPR2(adjmat, Z_Meta)

par(mfrow = c(2,1))
plotProp(props.raw, main = "Raw feature intensities", font.main = 1)
plotProp(props.adj, main = "Transformed to residuals of linear model of 
  intensity on BMI, place, diabetes status, fasting status", font.main = 1)

# Check PCA of transformed matrix
library(pca3d)
scores.adj <- prcomp(adjmat, scale. = F)
pca2d(scores.adj, group = scores$WEEKS)

#--------------------------------------------------------------------------------------------------

# Descriptions of files

# 1694 obs. of 8501 NMR variables (outliers removed, one NA variable in 8501st col)
#raw <- read_tsv("1510_XAlignedE3NcpmgssCitPEGfinal.txt") %>% select(-8501)

# 1739 obs. of 8500 NMR variables (all samples)
#raw1 <- read_delim("X_AlignedE3NData_cpmg_ssCitPEG_1112.txt", delim = ";")

# Metadata are stored in the following Excel file. The third sheet has the metadata w/o outliers
# (hope it's in the same order as the NMR data!)
#metaQC <- read_xlsx("1510_MatriceY_CohorteE3N_Appar.xlsx", sheet = 4, na = ".")






