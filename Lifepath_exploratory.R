# LIFEPATH data exploratory
library(tidyverse)
library(readxl)

# Final data files seem to be the following:
# 1623 observations of 44 intensity variables. Looks scaled version of dat 4 and final prepared data
#ints <- read_tsv("C:/J_ROTHWELL/1507_XMetabolite_std_cpmg_E3N.txt")
ints <- read_tsv("1507_XMetabolite_std_cpmg_E3N.txt")

# metadata (from XL or csv)
#meta <- read_xlsx("C:/J_ROTHWELL/1603_MatriceY_CohorteE3N_Appar.xlsx", sheet = 6)
meta <- read.csv("Lifepath_meta.csv")

# subset for baseline characteristics table
meta0 <- meta %>% 
  select(CT, AGE, BMI, HANCHE, MENOPAUSE, SMK, DIABETE, Life_Alcohol_Pattern_1, BP, Trait_Horm, 
         CO, CENTTIMECat1, FASTING, STOCKTIME, BEHAVIOUR, SUBTYPE, HR, Estro_THM, Pg_seul, 
         SBR, GRADE, STADE, DIAGSAMPLING) %>% 
  mutate_at(vars(-AGE, -BMI, -HANCHE, -DIAGSAMPLING, -STOCKTIME), as.factor)

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

# Exploratory analysis ----
# Check total intensities for each metabolite
plot(colSums(ints[ , -1]), xlab = "Compound number", ylab = "Scaled intensity",
     pch = 19, col = "dodgerblue", main = "Summed intensities of 44 metabolites")

# Check correlations
cormat <- cor(ints)
colnames(cormat) <- NULL
library(corrplot)
corrplot(cormat, method = "square", tl.col = "black", tl.cex = 0.8)

# The three fatty acids are highly correlated, valine and leucine, NAC1 and 2
# The fatty acids are inversely correlated with many compounds

# Run PCA and get proportions of variability for components
pca <- prcomp(ints[, -1], scale.=F)
plot(pca)

# Plot 
library(pca3d)
pca2d(pca)
title("Metabolite profiles of 1622 subjects")
box(which = "plot", lty = "solid")

# Find outlier and remove
which(pca$x[, 2] < -10) #row 1409

# Rerun PCA and replot
pca1 <- prcomp(ints[-1409, -1], scale.=F)
pca2d(pca1, biplot = T)
title("Metabolite profiles of 1622 subjects")
box(which = "plot", lty = "solid")

# Breast cancer risk model. Subset variables needed
meta1 <- meta %>%
  select(CODBMB, CT, MATCH, WEEKS, PLACE, AGE, BMI, MENOPAUSE, FASTING, SMK, DIABETE, CENTTIMECat1, SAMPYEAR, STOCKTIME) %>%
  mutate_at(vars(-CODBMB, -CT, -AGE, -BMI, -WEEKS, -STOCKTIME), as.factor)

# Conditional logistic regression to get odds ratios for lifestyle factors
library(survival)
fit <- clogit(CT ~ BMI + SMK + DIABETE + strata(MATCH), data = meta1) 
# output <- cbind(exp(coef(fit)), exp(confint(fit)))
library(broom)
t1 <- tidy(fit) %>% # mutate_if(is.numeric, exp) %>% 
  select(-(std.error:p.value))

library(metafor)
par(mar=c(5,4,1,2))
forest(t1$estimate, ci.lb = t1$conf.low, ci.ub = t1$conf.high, refline = 1, 
       xlab = "Multivariable adjusted odds ratio",
       transf = exp, pch = 18, psize = 1.5, slab = t1$term, alim = c(0,2), xlim = c(-1, 3))
text(-1, 5, "Variable", pos = 4)
text(3, 5, "OR [95% CI]", pos = 2)
# matching factors removed!

# CLR models to get odds ratios for metabolites
data <- left_join(meta1, ints, by = "CODBMB")
clr <- function(x) clogit(CT ~ x + BMI + SMK + DIABETE + strata(MATCH), data = meta1)
multifit <- apply(data[, 15:ncol(data)], 2, clr)
t2 <- map_df(multifit, tidy) %>% filter(term == "x")

par(mar=c(5,4,1,2))
forest(t2$estimate, ci.lb = t2$conf.low, ci.ub = t2$conf.high, refline = 1, #xlab = xtitle, 
       transf = exp, pch = 18, psize = 1, slab = names(multifit))
hh <- par("usr")
text(hh[1], nrow(t2) + 2, "Parameter", pos = 4)
text(hh[2], nrow(t2) + 2, "OR [95% CI]", pos = 2)
       #alim = c(0,2), xlim = c(-1, 3)) 


# Non-metabolite model lme4

library(lme4)
fit1 <- lmer(CT ~ BMI + SMK + DIABETE + (1|MATCH), data = meta1)

