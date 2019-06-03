# LIFEPATH data exploratory
library(tidyverse)
library(readxl)

# Read 1623 observations of 44 intensity variables (appears to be final scaled data) and metadata
ints <- read_tsv("1507_XMetabolite_std_cpmg_E3N.txt")
meta <- read.csv("Lifepath_meta.csv")

# subset IDs to get subjects included in CC. Get positions of final CC samples in metadata
samples <- ints$CODBMB %in% meta$CODBMB

# subset for baseline characteristics table (same order as manuscript)
meta0 <- meta %>% 
  select(CT,                    # case-control status
         AGE,                   
         BMICat1, 
         RTHCat1,               # Waist-hip ratio categorical
         MENOPAUSE,             # menopausal status at blood collection
         SMK, 
         DIABETE, 
         Life_Alcohol_Pattern_1, 
         BP, 
         CO,                    # taking oral contraceptives
         Trait_Horm,            # menopausal treatment therapy taken 24h before blood collection
         DURTHSDIAG,            # Duration of use of therapy at date of diagnosis
         CENTTIMECat1,          # time before centrifugation (?)
         FASTING, 
         STOCKTIME,             # Storage time (years)
         BEHAVIOUR,             # Tumour behaviour
         SUBTYPE, 
         CERB2,                 # HER2 receptor
         ER,                    # Estrogen receptor
         PR,                    # Progesterone receptor
         SBR, 
         GRADE, 
         STADE, 
         DIAGSAMPLINGCat1) %>% 
  mutate_at(vars(-AGE, -STOCKTIME), as.factor)

# Exploratory analysis ----
# Check total intensities for each metabolite
ints <- ints[samples, ]
plot(colSums(ints[ , -1]), xlab = "Compound number", ylab = "Scaled intensity",
     pch = 19, col = "dodgerblue", font.main = 1,
     main = "Summed intensities of 44 metabolites for 1582 subjects")

# Check correlations
library(corrplot)
cormat <- cor(ints[, -1])
colnames(cormat) <- NULL
corrplot(cormat, method = "square", tl.col = "black", tl.cex = 0.8)

# The three fatty acids are highly correlated, valine and leucine, NAC1 and 2
# The fatty acids are inversely correlated with many compounds

# Run PCA of all samples
pca <- prcomp(ints[, -1], scale.=F)
plot(pca)

# Plot 
library(pca3d)
par(mfrow = c(1, 2))
pca2d(pca)
title("A", font.main = 1, adj = 0)
box(which = "plot", lty = "solid")

# Run PCPR2
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
pca2 <- prcomp(adjmat, scale.=F)
pca2d(pca2, group = alldata$CT)
title("Residuals-adjusted metabolite\nprofiles of 1582 samples", font.main = 1)
box(which = "plot", lty = "solid")
legend("topleft", legend = plt$groups, col=plt$colors, pch=plt$pch)

# Breast cancer risk model. Subset variables needed
meta1 <- meta %>%
  select(CODBMB, CT, MATCH, PLACE, AGE, BMI, BP, RTH, ALCOHOL, MENOPAUSE, FASTING, SMK, DIABETE, CENTTIMECat1, SAMPYEAR, STOCKTIME) %>%
  mutate_at(vars(-CODBMB, -CT, -AGE, -BMI, -STOCKTIME, -RTH, -ALCOHOL), as.factor)

# Conditional logistic regression to get odds ratios for lifestyle factors
library(survival)
fit <- clogit(CT ~ BMI + SMK + DIABETE + BP + RTH + CENTTIMECat1 + STOCKTIME + strata(MATCH), data = meta1) 
# output <- cbind(exp(coef(fit)), exp(confint(fit)))
library(broom)
t1 <- tidy(fit) %>% # mutate_if(is.numeric, exp) %>% 
  select(-(std.error:p.value))

library(metafor)
dev.off()
par(mar=c(5,4,1,2))
forest(t1$estimate, ci.lb = t1$conf.low, ci.ub = t1$conf.high, refline = 1, 
       xlab = "Multivariable adjusted odds ratio",
       transf = exp, pch = 18, psize = 1.5, slab = t1$term, alim = c(0,2), xlim = c(-1, 3))
text(-1, nrow(t1) + 2, "Variable", pos = 4)
text(3, nrow(t1) + 2, "OR [95% CI]", pos = 2)
# matching factors removed!

# CLR models to get odds ratios for metabolites
data <- left_join(meta1, ints, by = "CODBMB")
clr <- function(x) clogit(CT ~ x + BMI + SMK + DIABETE + strata(MATCH), data = meta1)
ints <- 15:ncol(data)
multifit <- apply(data[, ints], 2, clr)
t2 <- map_df(multifit, tidy) %>% filter(term == "x")

dev.off()
par(mar=c(5,4,1,2))
forest(t2$estimate, ci.lb = t2$conf.low, ci.ub = t2$conf.high, refline = 1, #xlab = xtitle, 
       xlab = "Multivariable adjusted odds ratio",
       transf = exp, pch = 18, psize = 1, slab = names(multifit))
hh <- par("usr")
text(hh[1], nrow(t2) + 2, "Compound", pos = 4)
text(hh[2], nrow(t2) + 2, "OR [95% CI]", pos = 2)
       #alim = c(0,2), xlim = c(-1, 3)) 

# Funnel plots for metabolites
funnel(x = t2$estimate, sei = t2$std.error)
funnel(x = t2$estimate, sei = t2$std.error, yaxis = "vi")

# Non-metabolite model lme4
library(lme4)
fit1 <- glmer(CT ~ BMI + SMK + DIABETE + (1|MATCH), data = meta1, family = "binomial")
summary(fit1)

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



# Other files ----

# Rawest feature data
#raw <- read_delim("C:/J_ROTHWELL/X_AlignedCohorteE3NData_cpmg_ssCitPEG_0612.txt", delim = ";")

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



