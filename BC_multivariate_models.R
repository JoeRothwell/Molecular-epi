# OPLS multivariate models
library(tidyverse)

# Read 1623 observations of 44 intensity variables (appears to be final scaled data) and metadata
ints <- read_tsv("1507_XMetabolite_std_cpmg_E3N.txt")
meta <- read.csv("Lifepath_meta.csv")
alldata <- inner_join(meta, ints, by = "CODBMB")
ints1 <- alldata %>% select(`3Hydroxybutyrate`:Succinate)

# Adjust using residuals method
adj <- function(x) residuals(lm(x ~ BMI + SMK + DIABETE, data = alldata))
adjmat <- apply(ints1, 2, adj)

# OPLS-DA on residual-adjusted concentrations
library(mixOmics)
pca.lp <- pca(adjmat, ncomp = 10, center = F, scale = F)
plot(pca.lp)
plotIndiv(pca.lp, group = as.factor(alldata$CT), ind.names = T, legend = T, title = "Subject metabolic profiles")

# Run PLS-DA specifying number of components
plsda.res <- plsda(adjmat, as.factor(alldata$CT), ncomp = 5)

set.seed(2543) # for reproducibility here, only when the `cpus' argument is not used
perf.plsda <- perf(plsda.res, validation = "Mfold", folds = 5, progressBar = T, auc = T, nrepeat = 10) 

plot(perf.plsda, col = color.mixo(1:3), sd = T, legend.position = "horizontal")
plot(perf.plsda)

#coeff <- plsda.res$X
#plot(coeff[1, ])
perf.plsda$choice.ncomp # 2 components appear to be best

# Rerun the PLS-DA with 3 components
plsda.res1 <- plsda(adjmat, as.factor(alldata$CT), ncomp = 2)
plotVar(plsda.res1)
plotIndiv(plsda.res1, ind.names = F, legend = T, ellipse = T, title = 'PLS-DA')
plotLoadings(plsda.res1, contrib = "max", ndisplay = 50)
auroc(plsda.res1)