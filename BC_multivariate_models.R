# OPLS multivariate models (on annotated compounds)
library(tidyverse)

# Read 1623 observations of 44 intensity variables (appears to be final scaled data) and metadata
ints <- read_tsv("1507_XMetabolite_std_cpmg_E3N.txt")
meta <- read.csv("Lifepath_meta.csv")
dat <- inner_join(meta, ints, by = "CODBMB")
ints1 <- dat %>% select(`3Hydroxybutyrate`:Succinate) %>% as.matrix

# Adjust using residuals method (mixed effects works better)
library(lme4)
#adj <- function(x) residuals(lm(x ~ PLACE + AGE + BMI + DIABETE + CENTTIMECat1, data = dat))
adj1 <- function(x) residuals(lmer(x ~ AGE + BMI + DIABETE + SAMPYEAR + (1|PLACE), data = dat)) 

adjmat <- apply(ints1, 2, adj1)

Z_Meta <- dat %>%
  select(WEEKS, PLACE, AGE, BMI, MENOPAUSE, FASTING, SMK, DIABETE, CENTTIMECat1, SAMPYEAR, STOCKTIME) %>%
  mutate_at(vars(-AGE, -BMI, -WEEKS, -STOCKTIME), as.factor)

# Check with PCPR2 before and after adjustment
library(pcpr2)
props.raw <- runPCPR2(ints1, Z_Meta)
props.adj <- runPCPR2(adjmat, Z_Meta)
par(mfrow = c(2,1))
plotProp(props.raw, main = "Sources of variation in 1582 cases and controls")
plotProp(props.raw, main = "No adjustment to metabolite concentrations", font.main = 1)
plotProp(props.adj, main = "Residual-adjusted metabolite matrix", font.main = 1)


# Multivariate PLS with caret

library(caret)

# Compile cross-validation settings
set.seed(100)

# Outcome needs to be factor
dat$CT <- as.factor(dat$CT)
myfolds <- createMultiFolds(dat$CT, k = 5, times = 10)
control <- trainControl("repeatedcv", index = myfolds, selectionFunction = "oneSE")

# Train PLS model
dat1 <- data.frame(CT = as.factor(dat$CT), adjmat)

# Split into training and test sets
inTrain <- createDataPartition(dat1$CT, p=0.75, list = F)
training <- dat1[inTrain, ]
testing <- dat1[-inTrain, ]

# Run model on training set
# Note: works with tcControl but not trControl (gives strange error)
mod1 <- train(CT ~ ., data = training, method = "pls", metric = "Accuracy",
              tuneLength = 20, tcControl = control)  

plot(mod1) 
plot(varImp(mod1, scale = F), 10, main = "PLS-DA")  

predictions <- predict(mod1, newdata = testing)
confusionMatrix(predictions, testing$CT)
# Model is not better than random chance.



# With mixOmics
library(mixOmics)
pca.lp <- pca(adjmat, ncomp = 10, center = F, scale = F)
plot(pca.lp)
plotIndiv(pca.lp, group = as.factor(dat$CT), ind.names = T, legend = T, title = "Subject metabolic profiles")

# Run PLS-DA specifying number of components
plsda.res <- plsda(adjmat, as.factor(dat$CT), ncomp = 5)

set.seed(2543) # for reproducibility here, only when the `cpus' argument is not used
perf.plsda <- perf(plsda.res, validation = "Mfold", folds = 5, progressBar = T, auc = T, nrepeat = 10) 

plot(perf.plsda, col = color.mixo(1:3), sd = T, legend.position = "horizontal")
plot(perf.plsda)

#coeff <- plsda.res$X
#plot(coeff[1, ])
perf.plsda$choice.ncomp # 2 components appear to be best

# Rerun the PLS-DA with 3 components
plsda.res1 <- plsda(adjmat, as.factor(dat$CT), ncomp = 2)
plotVar(plsda.res1)
plotIndiv(plsda.res1, ind.names = F, legend = T, ellipse = T, title = 'PLS-DA')
plotLoadings(plsda.res1, contrib = "max", ndisplay = 50)
auroc(plsda.res1)


