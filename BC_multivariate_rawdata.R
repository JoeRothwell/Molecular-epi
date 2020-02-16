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

# First give each control the time to diagnosis time of the corresponding case
meta <- meta %>% group_by(MATCH) %>% mutate(tdiag = max(as.numeric(DIAGSAMPLINGCat1), na.rm = T))
all <- data.frame(class = as.factor(meta$CT), adjmat)

# Logical vectors (need to group for diagnosed > and < 5 years)
pre   <- meta$MENOPAUSE == 0
post  <- meta$MENOPAUSE == 1 
early <- meta$tdiag == 1
late  <- meta$tdiag == 2


library(caret)

# 0. All subjects
# Split into training and test sets on class, 75% to training set
inTrain <- createDataPartition(y = all$class, p = 0.75, list = F)
training <- all[inTrain, ]
testing <- all[-inTrain, ]

# Cross validation. The sample size is quite large so can use a large number of folds (10)
set.seed(111)
#folds <- createMultiFolds(y = training$class, k = 10, times = 5)
folds <- createFolds(y = training$class, k = 10)
control <- trainControl("repeatedcv", index = folds, selectionFunction = "oneSE")
sapply(folds, length)

# Train PLS model
mod0 <- train(class ~ ., data = training, method = "pls", metric = "Accuracy", 
             trControl = control, tuneLength = 20) 

plot(mod0)
confusionMatrix(mod0)
#plot(varImp(mod0), 10)

# Predict test set
predictions0 <- predict(mod0, newdata = testing)
confusionMatrix(predictions0, reference = testing$class)

# Get AUC
library(pROC)
predictions0_1 <- predict(mod0, newdata = testing, type = "prob")
result0 <- roc(testing$class, predictions0_1$`0`, ci = T)
ci.auc(result0)

# Plot ROC curve
plot(result0, main = "All subjects")
plot.roc(result0, main = "All subjects", ci = T, auc.polygon.col = "lightblue")


#-----------------------------------------------------------------------------------------

# 1. Post-menopausal only
# Split into training and test sets on class, 75% to training set
post <- all[post, ]
inTrain <- createDataPartition(y = post$class, p = 0.75, list = F)
training <- post[inTrain, ]
testing <- post[-inTrain, ]

set.seed(222)
folds <- createFolds(y = training$class, k = 10)
control <- trainControl("repeatedcv", index = folds, selectionFunction = "oneSE")

# Train PLS model
mod1 <- train(class ~ ., data = training, method = "pls", metric = "Accuracy", 
             trControl = control, tuneLength = 20) 

plot(mod1)
#confusionMatrix(mod1)
#plot(varImp(mod1), 10)

# Predict test set
predictions1 <- predict(mod1, newdata = testing)
confusionMatrix(predictions1, reference = testing$class)

# Get AUC
library(pROC)
predictions1_1 <- predict(mod1, newdata = testing, type = "prob")
result1 <- roc(testing$class, predictions1_1$`0`, ci = T)
ci.auc(result1)

# Plot ROC curve
plot(result1, main = "Post-menopausal")


#---------------------------------------------------------------------------------------

# 2. Pre-menopausal only
# Split into training and test sets on class, 75% to training set
pre <- all[pre, ]
inTrain <- createDataPartition(y = pre$class, p = 0.75, list = F)
training <- pre[inTrain, ]
testing <- pre[-inTrain, ]

# Use a 10 times repeated 5-fold cross validation
set.seed(333)
folds <- createMultiFolds(y = training$class, k = 5, times = 5)
control <- trainControl("repeatedcv", index = folds, selectionFunction = "oneSE")

# Train PLS model
mod2 <- train(class ~ ., data = training, method = "pls", metric = "Accuracy", 
             trControl = control, tuneLength = 20) 

plot(mod2)
#confusionMatrix(mod2)
#plot(varImp(mod2), 10)

# Predict test set
predictions2 <- predict(mod2, newdata = testing)
confusionMatrix(predictions2, reference = testing$class)

# Get AUC
library(pROC)
predictions2_1 <- predict(mod2, newdata = testing, type = "prob")
result2 <- roc(testing$class, predictions2_1$`0`, ci = T)
ci.auc(result2)

# Plot ROC curve
plot(result2, main = "Pre-menopasual")

#---------------------------------------------------------------------------------------

# 3. Early diagnosis only
# Split into training and test sets on class, 75% to training set
early <- all[early, ]
inTrain <- createDataPartition(y = early$class, p = 0.75, list = F)
training <- early[inTrain, ]
testing <- early[-inTrain, ]

# Use a 10 times repeated 5-fold cross validation
set.seed(444)
folds <- createFolds(y = training$class, k = 10)
control <- trainControl("repeatedcv", index = folds, selectionFunction = "oneSE")

# Train PLS model
mod3 <- train(class ~ ., data = training, method = "pls", metric = "Accuracy", 
             trControl = control, tuneLength = 20) 

plot(mod3)
#confusionMatrix(mod3)
#plot(varImp(mod3), 10)

# Predict test set
predictions3 <- predict(mod3, newdata = testing)
confusionMatrix(predictions3, reference = testing$class)

# Get AUC
library(pROC)
predictions3_1 <- predict(mod3, newdata = testing, type = "prob")
result3 <- roc(testing$class, predictions3_1$`0`, ci = T)
ci.auc(result3)

# Plot ROC curve
plot(result3, main = "Diagnosis < 5y")

#---------------------------------------------------------------------------------------

# 4. Late diagnosis only
# Split into training and test sets on class, 75% to training set
late <- all[late, ]
set.seed(555)
inTrain <- createDataPartition(y = late$class, p = 0.75, list = F)
training <- late[inTrain, ]
testing <- late[-inTrain, ]

# Use a 10 times repeated 5-fold cross validation
set.seed(555)
folds <- createFolds(y = training$class, k = 10)
control <- trainControl("repeatedcv", index = folds, selectionFunction = "oneSE")

# Train PLS model
set.seed(555)
mod4 <- train(class ~ ., data = training, method = "pls", metric = "Accuracy", 
             trControl = control, tuneLength = 20) 

plot(mod4)
#confusionMatrix(mod4)
#plot(varImp(mod), 10)

# Predict test set
predictions4 <- predict(mod4, newdata = testing)
confusionMatrix(predictions4, reference = testing$class)

# Get AUC
library(pROC)
predictions4_1 <- predict(mod4, newdata = testing, type = "prob")
result4 <- roc(testing$class, predictions4_1$`0`, ci = T)
result4
ci.auc(result4)

# Plot ROC curve
plot(result4, main = "Diagnosis > 5y")


#----------------------------------------------------------------------------------------

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






