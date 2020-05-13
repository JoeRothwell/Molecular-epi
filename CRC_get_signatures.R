# Compute Biocrates and fatty acid signatures of WCRF score by PLS
source("CRC_prep_data.R")
library(tidyverse)
library(lme4)
library(zoo)

# Get compound matrix and questionnaire scores
# Adjusts and scales the controls metabolite matrices for confounders and binds scores for model
# Set zeros to NA and impute with half miminum conc, log transform
# Bind WCRF scores to adjusted metabolite matrix for PLS

get.plsdata <- function(dat, ctrldat, bioc = T, subgroup = F){
  
  if(bioc == T & subgroup == F) {
    # Adjust matrix for study, centre, batch, sex for Biocrates, subset calibrated scores 
    concs <- as.matrix(dat)
    adj   <- function(x) residuals(lmer(x ~ Center + batch_no + Sex + (1|Study), data = ctrldat))
    
  } else if(bioc == T & subgroup == T) {
    concs <- as.matrix(dat)
    adj   <- function(x) residuals(lmer(x ~ Center + batch_no + (1|Study), data = ctrldat))
    
  } else {
    
    concs <- ctrldat %>% select(one_of(common.cols))
    adj   <- function(x) residuals(lmer(x ~ LABO + STUDY + (1|Center), data = ctrldat))
  }
  
  score <- tibble(score = ctrldat$Wcrf_C_Cal)  
  # Prepare controls matrix. Replace zero, impute with half mins, scale
  concs[concs == 0] <- NA
  concs1 <- na.aggregate(concs, FUN = function(x) min(x)/2)
  logconcs <- log2(concs1) %>% scale
  adjmat <- apply(logconcs, 2, adj)
  
  # Data setup. Must be a df with Bind scores to log matrix
  output <- cbind(score, adjmat) %>% filter(!is.na(score))
  
}

# Get PLS data for whole controls dataset (manuscript table, more compounds)
Bioc0  <- get.plsdata(ctrls, ctrl) # All control compounds

# Get PLS data for case-control risk models
Bioc1  <- get.plsdata(ctrlA, ctrl, bioc = T) # Overlap control/A
Bioc2  <- get.plsdata(ctrlB, ctrl, bioc = T) # Overlap control/B
FAdata <- get.plsdata(crc1fa, fa.ctrl, bioc = F) # Fatty acids
#Bioc3  <- get.plsdata(ctrls0) # Overlap control/A/B (not needed)

# PLS data for sex-specific signatures (not done for FAs because too few males)
#Bioc1m <- get.plsdata(ctrlAm, ctrl.m, bioc = T, subgroup = T)
#Bioc1f <- get.plsdata(ctrlAf, ctrl.f, bioc = T, subgroup = T)
#Bioc2m <- get.plsdata(ctrlBm, ctrl.m, bioc = T, subgroup = T)
#Bioc2f <- get.plsdata(ctrlBf, ctrl.f, bioc = T, subgroup = T)
#FAdatf <- get.plsdata(crc1faf, fa.ctrl.f, bioc = F)

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

# Sex-specific
mod1m <- plsr(score ~ ., data = Bioc1m, ncomp = 1)
mod1f <- plsr(score ~ ., data = Bioc1f, ncomp = 1)
mod2m  <- plsr(score ~ ., data = Bioc2m, ncomp = 1)
mod2f  <- plsr(score ~ ., data = Bioc2f, ncomp = 1)
modFAf <- plsr(score ~ ., data = FAdatf, ncomp = 2)

# explained variances, prediction, scores, loadings plots
# explvar(mod)
# plot(mod, plottype = "scores")


# Produce tables of important compounds, using compound metadata to get proper names
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

# Finally, predict WCRF scores from Biocrates or fatty acids data for two datasets
# (three predictions)

# Small case-control
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