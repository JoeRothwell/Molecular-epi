# Breast cancer study multivariate models from unbinned NMR data (raw features)

# Data from Elodie Jobard 27-6-2019
# Read in data with outliers and QCs (n=1739) or with QCs only
# Read in NMR data and metadata. Dataset containing samples and QCs only 

library(tidyverse)
library(readxl)
dat <- read_tsv("1510_XAlignedE3NcpmgssCitPEGfinal.txt") %>% select(-8501)
meta <- read_xlsx("1510_MatriceY_CohorteE3N_Appar.xlsx", sheet = 4, na = ".")

# Dataset containing all samples, QCs and blanks
#mat <- read_delim("X_AlignedE3NData_cpmg_ssCitPEG_1112.txt", delim = ";", n_max = 1738) 

# Subset samples only and plot PCA with pareto scaling
prep.data <- function(incl.qc = F, pc.scores = T) {

  samp <- meta$TYPE_ECH == 1 & !is.na(meta$CENTTIME)
    
  # Subset samples only (removing 112 QCs) from data and metadata
  mat  <- if(incl.qc == F) dat[samp, ] %>% as.matrix else as.matrix(dat)
  meta <- if(incl.qc == F) meta[samp, ]

  #anyNA(mat)
  #sum(apply(mat, 2, anyNA))
  
  # remove zero variance columns
  zerovar <- sum(apply(mat, 2, var) == 0)
  print(paste("There are", zerovar, "zero variance columns"))
  
  # There are 1116. Place in logical vector
  nonzerovar <- apply(mat, 2, var) != 0
  which(nonzerovar)
  
  mat0 <- mat[ , nonzerovar]
  print(paste("Dimensions", dim(mat0)))

  # Scale and run PCA  
  library(MetabolAnalyze)
  scalemat <- scaling(mat0, type = "pareto")
  
  if(pc.scores == F) return(list(scalemat, meta))
  pca <- prcomp(scalemat, scale. = F, center = T, rank. = 10)
  output <- data.frame(pca$x) %>% bind_cols(meta)
}

scores <- prep.data()
scores1 <- prep.data(incl.qc = T)

# Scores plot for manuscript
library(ggplot2)
ggplot(scores, aes(PC1, PC2, colour = as.factor(TYPE_ECH))) + geom_point() + theme_bw() +
  xlab("Score on PC1") + ylab("Score on PC2") +
  scale_color_discrete(labels = c("Experimental samples", "QCs")) +
  theme(legend.position = "bottom", #legend.justification = c(0, 0), 
        legend.title = element_blank()) +
  geom_hline(yintercept = 0, linetype = "dashed") +
  geom_vline(xintercept = 0, linetype = "dashed")

library(pca3d)
pca2d(as.matrix(scores[, 1:10]), group = scores$WEEKS)
box(which = "plot", lty = "solid")

# ----------------------------------------------------------------------------------------------------------

# Output Pareto-scaled data and perform PCPR2

dat <- prep.data(pc.scores = F)
#dat <- prep.data(data.only = T, exclude.tbc = T)

concs <- dat[[1]]
meta <- dat[[2]] %>%
  select(CT, MATCH, WEEKS, PLACE, AGE, BMI, MENOPAUSE, FASTING, SMK, DIABETE, CENTTIMECat1, CENTTIME, SAMPYEAR, 
         DIAGSAMPLINGCat1, STOCKTIME) %>%
  mutate_at(vars(-AGE, -BMI, -CENTTIME), as.factor)

Z_Meta <- meta %>% select(-MATCH, -CT, -CENTTIME, -DIAGSAMPLINGCat1)

library(pcpr2)
props.raw <- runPCPR2(concs, Z_Meta)
plot(props.raw)

# Greatest sources of variability are BMI > DIABETE > PLACE
# Adjust for fixed effects only. Random effects model with lme4 did not work, boundary fit or didn't converge
adj <- function(x) residuals(lm(x ~ PLACE + WEEKS + AGE + BMI + DIABETE + FASTING + SAMPYEAR + 
                                  CENTTIMECat1 + STOCKTIME, data = meta))
adjmat <- apply(concs, 2, adj)
#saveRDS(adjmat, "adjusted_NMR_features.rds")

props.adj <- runPCPR2(adjmat, Z_Meta)

par(mfrow = c(1, 2))
plot(props.raw, main = "Raw feature intensities", font.main = 1)
plot(props.adj, main = "Transformed to residuals of linear model of 
  intensity on confounders*", font.main = 1)

# Check PCA of transformed matrix
library(pca3d)
scores.adj <- prcomp(adjmat, scale. = F)
pca2d(scores.adj, group = scores$WEEKS)

# Final data matrix is adjmat, n = 1572

#--------------------------------------------------------------------------------------------------

# Multivariate analysis. Subsets to be made:
# 1. All samples; 2. Pre-menopausal only; 3. Post-menopausal only; 4. Diagnosed < 5 years only; 5. Diagnosed > 5 years only

adjmat <- load("adjusted_NMR_features.rds")

# First give each control the time to diagnosis time of the corresponding case
meta <- meta %>% group_by(MATCH) %>% mutate(tdiag = max(as.numeric(DIAGSAMPLINGCat1), na.rm = T))
all <- data.frame(class = as.factor(meta$CT), adjmat)

# Logical vectors (need to group for diagnosed > and < 5 years)
pre   <- meta$MENOPAUSE == 0
post  <- meta$MENOPAUSE == 1 
early <- meta$tdiag == 1
late  <- meta$tdiag == 2

library(caret)
library(pROC)
# Function to train and fit the model and do a ROC analysis
bc.roc <- function(dat, ...) {
  
  # Split into training and test sets on class, 75% to training set
  inTrain <- createDataPartition(y = dat$class, p = 0.75, list = F)
  training <- dat[inTrain, ]
  testing <- dat[-inTrain, ]
  
  # Cross validation. If the sample size is large can use a large number of folds (10)
  set.seed(111)
  folds <- createMultiFolds(y = training$class, ...)
  control <- trainControl("repeatedcv", index = folds, selectionFunction = "oneSE")
  print(sapply(folds, length))
  
  # Train PLS model
  mod0 <- train(class ~ ., data = training, method = "pls", metric = "Accuracy", 
                trControl = control, tuneLength = 20) 
  
  #plot(mod0)
  confusionMatrix(mod0)
  
  # Predict test set
  predictions0 <- predict(mod0, newdata = testing)
  confusionMatrix(predictions0, reference = testing$class)
  
  # Get AUC
  predictions0_1 <- predict(mod0, newdata = testing, type = "prob")
  result0 <- roc(testing$class, predictions0_1$`0`, ci = T)
  
}

p0 <- bc.roc(all, k = 10)
p1 <- bc.roc(all[post, ], k = 10)
p2 <- bc.roc(all[pre, ], k = 5, times = 5)
p3 <- bc.roc(all[early, ], k = 10)
p4 <- bc.roc(all[late, ], k = 10)

plot.roc(p0, ci = T, grid = T, print.auc = T)
plot.roc(p1, grid = T, print.auc = T)
plot.roc(p2, grid = T, print.auc = T)
plot.roc(p3, grid = T, print.auc = T)
plot.roc(p4, grid = T, print.auc = T)

library(ggROC)
ggroc(p4, colour = "darkblue", size = 1) + theme_bw() +
  geom_segment(aes(x = 1, xend = 0, y = 0, yend = 1), color="grey", linetype="dashed") +
  ggsave("roc_late.png")

# Compare models
models <- resamples(list("All" = mod0, "Post-menopausal" = mod1,
                         "Early diagnosis" = mod3, "Late diagnosis" = mod4))
bwplot(models, metric = "Accuracy")

# Description of files

# 1694 obs. of 8501 NMR variables (outliers removed, one NA variable in 8501st col)
#raw <- read_tsv("1510_XAlignedE3NcpmgssCitPEGfinal.txt") %>% select(-8501)

# 1739 obs. of 8500 NMR variables (all samples)
#raw1 <- read_delim("X_AlignedE3NData_cpmg_ssCitPEG_1112.txt", delim = ";")

# Metadata are stored in the following Excel file. The third sheet has the metadata w/o outliers
# (hope it's in the same order as the NMR data!)
#metaQC <- read_xlsx("1510_MatriceY_CohorteE3N_Appar.xlsx", sheet = 4, na = ".")






