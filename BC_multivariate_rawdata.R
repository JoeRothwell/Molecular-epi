# Breast cancer study multivariate models from unbinned NMR data (raw features)

# Data from Elodie Jobard 27-6-2019
# Read in data with outliers and QCs (n=1739) or with QCs only

# Subset samples only and plot PCA with pareto scaling
prep.data <- function(dataset = "QCs", data.only = F, exclude.tbc = F) {
  
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
    samp <- if(exclude.tbc == F) metaQC$TYPE_ECH == 1 else metaQC$TYPE_ECH == 1 & !is.na(metaQC$CENTTIME)
    mat <- dat[samp, ] %>% as.matrix
    #mat1 <- dat[samp1, ] %>% as.matrix
    
    # Subset samples only from metadata   
    meta <- metaQC[samp, ]
    #meta1 <- metaQC[samp1, ]
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

scores <- prep.data()
#scores1 <- explore.data(dataset = "All")

# Plot data
library(ggplot2)
ggplot(scores, aes(PC1, PC2)) + geom_text(aes(label = WEEKS)) + theme_bw()

library(pca3d)
pca2d(as.matrix(scores[, 1:10]), group = scores$WEEKS)
box(which = "plot", lty = "solid")

# ----------------------------------------------------------------------------------------------------------

# Output Pareto-scaled data and perform PCPR2

dat <- prep.data(data.only = T, exclude.tbc = T)

concs <- dat[[1]]
meta <- dat[[2]] %>%
  select(CT, MATCH, WEEKS, PLACE, AGE, BMI, MENOPAUSE, FASTING, SMK, DIABETE, CENTTIMECat1, CENTTIME, SAMPYEAR, 
         DIAGSAMPLINGCat1, STOCKTIME) %>%
  mutate_at(vars(-AGE, -BMI, -CENTTIME), as.factor)

Z_Meta <- meta %>% select(-MATCH, -CT, -CENTTIME, -DIAGSAMPLINGCat1)

library(pcpr2)
props.raw <- runPCPR2(concs, Z_Meta)
plotProp(props.raw)

# Greatest sources of variability are BMI > DIABETE > PLACE
# Adjust for fixed effects only. Random effects model with lme4 did not work, boundary fit or didn't converge
adj <- function(x) residuals(lm(x ~ PLACE + WEEKS + AGE + BMI + DIABETE + FASTING + SAMPYEAR + 
                                  CENTTIMECat1 + STOCKTIME, data = meta))
adjmat <- apply(concs, 2, adj)

props.adj <- runPCPR2(adjmat, Z_Meta)

par(mfrow = c(2,1))
plotProp(props.raw, main = "Raw feature intensities", font.main = 1)
plotProp(props.adj, main = "Transformed to residuals of linear model of 
  intensity on confounders*", font.main = 1)

# Check PCA of transformed matrix
library(pca3d)
scores.adj <- prcomp(adjmat, scale. = F)
pca2d(scores.adj, group = scores$WEEKS)

# Final data matrix is adjmat, n = 1572

#--------------------------------------------------------------------------------------------------

# Multivariate analysis. Subsets to be made:
# 1. All samples; 2. Pre-menopausal only; 3. Post-menopausal only; 4. Diagnosed < 5 years only; 5. Diagnosed > 5 years only

# First give each control the time to diagnosis time of the corresponding case
meta <- meta %>% group_by(MATCH) %>% mutate(tdiag = max(as.numeric(DIAGSAMPLINGCat1), na.rm = T))
all <- data.frame(class = as.factor(meta$CT), adjmat)

# Logical vectors (need to group for diagnosed > and < 5 years)
pre   <- meta$MENOPAUSE == 0
post  <- meta$MENOPAUSE == 1 
early <- meta$tdiag == 1
late  <- meta$tdiag == 2

# Get logical vectors for pre and post menopaual, < 5 and > 5 years

library(caret)
# Split into training and test sets on class, 75% to training set
inTrain <- createDataPartition(y = all$class, p = 0.75, list = F)
training <- all[inTrain, ]
testing <- all[-inTrain, ]

# Cross validation. The sample size is quite large so can use a large number of folds (10)
set.seed(111)
folds <- createMultiFolds(y = training$class, k = 10, times = 5)
control <- trainControl("repeatedcv", index = folds, selectionFunction = "oneSE")
sapply(folds, length)

# Train PLS model
set.seed(111)
mod <- train(class ~ ., data = all, method = "pls", metric = "Accuracy", 
             trControl = control, tuneLength = 20) 

plot(mod)
confusionMatrix(mod)
plot(varImp(mod), 10)

# Predict test set
predictions <- predict(mod, newdata = testing)
confusionMatrix(predictions, reference = testing$class)



# Description of files

# 1694 obs. of 8501 NMR variables (outliers removed, one NA variable in 8501st col)
#raw <- read_tsv("1510_XAlignedE3NcpmgssCitPEGfinal.txt") %>% select(-8501)

# 1739 obs. of 8500 NMR variables (all samples)
#raw1 <- read_delim("X_AlignedE3NData_cpmg_ssCitPEG_1112.txt", delim = ";")

# Metadata are stored in the following Excel file. The third sheet has the metadata w/o outliers
# (hope it's in the same order as the NMR data!)
#metaQC <- read_xlsx("1510_MatriceY_CohorteE3N_Appar.xlsx", sheet = 4, na = ".")






