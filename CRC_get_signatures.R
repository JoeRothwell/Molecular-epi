# Compute Biocrates and fatty acid signatures of WCRF score by PLS
source("CRC_prep_data.R")
library(tidyverse)
library(lme4)
library(zoo)

# Get compound matrix and questionnaire scores
# Adjusts and scales the controls metabolite matrices for confounders and binds scores for model
# Define linear models for removal of confounders from matrix
adj   <- function(x) residuals(lmer(x ~ Center + batch_no + Sex + (1|Study), data = ctrl))
adjFA <- function(x) residuals(lmer(x ~ LABO + STUDY + (1|Center), data = fa.ctrl))

# Replace zero, impute with half mins, scale, transform to residuals of models above
logmat0 <- ctrls %>% na_if(0) %>% na.aggregate(FUN = function(x) min(x)/2) %>% log2 %>% scale
adjmat0 <- apply(logmat0, 2, adj) %>% data.frame # All control compounds

logmat1 <- ctrlA %>% na_if(0) %>% na.aggregate(FUN = function(x) min(x)/2) %>% log2 %>% scale
adjmat1 <- apply(logmat1, 2, adj) %>% data.frame # Overlap control/A

logmat2 <- ctrlB %>% na_if(0) %>% na.aggregate(FUN = function(x) min(x)/2) %>% log2 %>% scale
adjmat2 <- apply(logmat2, 2, adj) %>% data.frame # Overlap control/B

logmat3 <- fa.ctrl[, common.cols] %>% na_if(0) %>% na.aggregate(FUN = function(x) min(x)/2) %>% 
  log2 %>% scale
adjmat3 <- apply(logmat3, 2, adjFA) %>% data.frame # Fatty acids

# Bind WCRF scores to adjusted metabolite matrix for PLS modelling
Bioc0 <- cbind(score = ctrl$Wcrf_C_Cal, adjmat0) %>% filter(!is.na(score))
Bioc1 <- cbind(score = ctrl$Wcrf_C_Cal, adjmat1) %>% filter(!is.na(score))
Bioc2 <- cbind(score = ctrl$Wcrf_C_Cal, adjmat2) %>% filter(!is.na(score))
FAdata <- cbind(score = fa.ctrl$Wcrf_C_Cal, adjmat3) %>% filter(!is.na(score))
# Overlap control/A/B (not needed)

# Function to get optimal dimensions for each model (pls or caret)
library(pls)
get.signature <- function(plsdata, which.mod = "plsmod"){

    if(which.mod == "plsmod"){

    # Start with a sensible number of components eg 10
    set.seed(111)
    mod <- plsr(score ~ ., ncomp = 10, data = plsdata, validation = "CV")
    # Find the number of dimensions with lowest cross validation error
    cv <- RMSEP(mod)
    plot(RMSEP(mod), legendpos = "topright")
    
    # Get components with lowest RMSEP, "one SE" and "permutation" methods
    best.dims  <- which.min(cv$val[estimate = "adjCV", , ]) - 1
    onesigma   <- selectNcomp(mod, method = "onesigma", plot = F)
    permut     <- selectNcomp(mod, method = "randomization", plot = T)
    
    print(paste("Lowest RMSEP from", best.dims, "comp(s);", 
                "one SE method suggests", onesigma, "comp(s);",
                "permutation method suggests", permut, "comp(s)"))
    
  } else if(which.mod == "caretmod") {  
    
    library(caret)
    set.seed(111)
    mod <- train(score ~ ., data = plsdata, method = "pls", metric = "RMSE", tuneLength = 20) 
  }
  
}

# Fit final PLS models with optimal dimensions to get signatures (see PLS vignette p12)
set.seed(111)
mod1a <- plsr(score ~ ., data = Bioc1, ncomp = 1)
mod1b <- plsr(score ~ ., data = Bioc2, ncomp = 1)
mod2  <- plsr(score ~ ., data = FAdata, ncomp = 2)
mod0  <- plsr(score ~ ., data = Bioc0, ncomp = 1)
# explained variances, prediction, scores, loadings plots
# explvar(mod)
# plot(mod, plottype = "scores")

# Make tables of important compounds, using compound metadata to get proper names
plot.signature <- function(mod, biocrates = T, percentile = 5, all = T){
  
  library(tidyverse)
  cmpds <- if(biocrates == T) read_csv("Biocrates_cmpd_metadata.csv") else read_csv("FA_compound_data.csv")
  
  # Coefficients and variable importance. First subset one-matrix array to get matrix
  # Extract coefficients from pls or caret objects, 1st LV: (pls gives an mvr object)
  if(class(mod) == "mvr"){
    coeff <- data.frame(value = round(coef(mod)[ , 1, 1], 3))
  } else {
    coeff <- data.frame(value = mod$finalModel$coefficients[, 1, 1])
  }
  
  pc <- 100/percentile
  
  dat <- coeff %>%
    rownames_to_column(var = "Compound") %>% mutate(sm = sum(abs(value))) %>%
    left_join(cmpds, by = "Compound") %>% 
    # Subset appropriate columns and print influential compounds
    select(class, compound = displayname, Coefficient = value)
  
  # Get top and bottom 5 percentiles and subset for table
  dat <- dat %>% mutate(Perc = cut_number(Coefficient, n = pc, labels = 1:pc)) #%>%
    #arrange(class)
  dat1 <- dat %>% filter(Perc == 1 | Perc == pc) %>% 
    arrange(Perc, abs(Coefficient)) %>% map_df(rev)
  
  if(all == T) return(dat) else return(dat1)
  
  n.high = sum(dat1$Perc == pc)
  n.low <- sum(dat1$Perc == 1)
  lv <- if(class(mod) == "mvr") mod$ncomp else mod$finalModel$ncomp

}

# Data for manuscript tables (copy into ms file via Excel)
table3a <- plot.signature(mod0, percentile = 5, all = F)
table3b <- plot.signature(mod2, biocrates = F, percentile = 10, all = F)

# Data for scatter plots
pltdata <- plot.signature(mod0, all = T)
faplot  <- plot.signature(mod2, biocrates = F, percentile = 10, all = T)

# Save workspace (for .Rmd file)
#save.image(file="metabolic_signatures.Rdata")

# Finally, predict WCRF scores from Biocrates or fatty acids data for two datasets (x3)
# Small case-control (no zero intensities)
df <- crc1[, colnames(ctrlA)] %>% log2 %>% scale %>% as_tibble
crc1.ph <- cbind(crc1, comp1 = predict(mod1a, df)[,,1])

# Large case-control
df1 <- crc2[, colnames(ctrlB)] %>% na_if(0) %>% na.aggregate(FUN = function(x) min(x)/2) %>%
  log2 %>% scale %>% as_tibble
crc2.ph <- cbind(crc2, comp1 = predict(mod1b, df1)[,,1])

# Small case-control fatty acids
df2 <- crc1fa[, common.cols] %>% na_if(0) %>% na.aggregate(FUN = function(x) min(x)/2) %>%
  log2 %>% scale %>% as_tibble
crc3.ph <- cbind(crc1fa, comp2 = predict(mod2, df2)[,,1])


# Old: sex specific signatures (not used)

# PLS data for sex-specific signatures (not done for FAs because too few males)
#Bioc1m <- get.plsdata(ctrlAm, ctrl.m, bioc = T, subgroup = T)
#Bioc1f <- get.plsdata(ctrlAf, ctrl.f, bioc = T, subgroup = T)
#Bioc2m <- get.plsdata(ctrlBm, ctrl.m, bioc = T, subgroup = T)
#Bioc2f <- get.plsdata(ctrlBf, ctrl.f, bioc = T, subgroup = T)
#FAdatf <- get.plsdata(crc1faf, fa.ctrl.f, bioc = F)

# Sex-specific
#mod1m <- plsr(score ~ ., data = Bioc1m, ncomp = 1)
#mod1f <- plsr(score ~ ., data = Bioc1f, ncomp = 1)
#mod2m  <- plsr(score ~ ., data = Bioc2m, ncomp = 1)
#mod2f  <- plsr(score ~ ., data = Bioc2f, ncomp = 1)
#modFAf <- plsr(score ~ ., data = FAdatf, ncomp = 2)