# Compute Biocrates and fatty acid signatures of WCRF score by PLS
source("CRC_prep_data.R")

# Get compound matrix and questionnaire scores
get.plsdata <- function(dat, cor.data = F){
 # changes
  # Adjusts and scales the controls metabolite matrices for confounders and binds scores for model
  library(tidyverse)
  library(lme4)
  library(zoo)
  # Set zeros to NA and impute with half miminum conc, log transform
  # Bind WCRF scores to adjusted metabolite matrix
  
  if(nrow(dat) == 877) {
  
  # select common FAs in the same order. FAs is for correlation with Biocrates compounds
  FAs  <- dat %>% select(Idepic, one_of(common.cols))
  if(cor.data == T) return(FAs)
  
  CRCfa <- FAs %>% select(-Idepic)
  concs <- fa.ctrl %>% select(one_of(common.cols))
  
  # Fatty acids: adjust matrix for study, centre, sex, batch
  adj <- function(x) residuals(lmer(x ~ LABO + STUDY + (1|Center), data = fa.ctrl))
  # Subset calibrated scores
  score <- data_frame(score = fa.ctrl$Wcrf_C_Cal)
  
  } else {
    
  concs <- as.matrix(dat)
  # Biocrates: adjust matrix for study, centre, sex, batch
  adj <- function(x) residuals(lmer(x ~ Center + batch_no + Sex + (1|Study), data = ctrl))
  # Subset calibrated scores
  score <- data_frame(score = ctrl$Wcrf_C_Cal)
  
  }
  
  # Prepare controls matrix. Replace zero, impute with half mins, scale
  concs[concs == 0] <- NA
  concs1 <- na.aggregate(concs, FUN = function(x) min(x)/2)
  logconcs <- log2(concs1) %>% scale
  
  adjmat <- apply(logconcs, 2, adj)
  
  # Data setup. Must be a df with Bind scores to log matrix
  output <- cbind(score, adjmat) %>% filter(!is.na(score))
  
}

# Get PLS data for whole controls dataset (manuscript table, more compounds)
Bioc0  <- get.plsdata(ctrls) # All control compounds

# Get PLS data for case-control risk models
Bioc1  <- get.plsdata(ctrlA) # Overlap control/A
Bioc2  <- get.plsdata(ctrlB) # Overlap control/B
FAdata <- get.plsdata(CRCfa1) # Fatty acids
#Bioc3  <- get.plsdata(ctrls0) # Overlap control/A/B (not needed)

# Get metabolic signatures by PLSR for case-control risk models
library(pls)
get.signature <- function(plsdata, which.mod = "plsmod"){
  
  # Biocrates (endogenous metabolites) for metabolic signature of WCRF score on controls dataset
  # This function gives the choice of fitting the model either with the pls package or Caret
  # pls was eventually used for the analysis.
  
  if(which.mod == "plsmod"){
    
    set.seed(111)
    # Start with a sensible number of components eg 10
    mod <- plsr(score ~ ., ncomp = 10, data = plsdata, validation = "CV")
    # Find the number of dimensions with lowest cross validation error
    cv <- RMSEP(mod)
    plot(RMSEP(mod), legendpos = "topright")
    
    # Calculate optimal number of dimensions and rerun model (see PLS vignette p12)
    
    # Get components with lowest RMSEP, "one SE" and "permutation" methods
    best.dims <- which.min(cv$val[estimate = "adjCV", , ]) - 1
    onesigma <- selectNcomp(mod, method = "onesigma", plot = F)
    permut <- selectNcomp(mod, method = "randomization", plot = T)
    
    print(paste("Lowest RMSEP from", best.dims, "comp(s);", 
                "one SE method suggests", onesigma, "comp(s);",
                "permutation method suggests", permut, "comp(s)"))
    
    #mod <- plsr(score ~ ., data = plsdata, ncomp = 3)
    
    # explained variances
    # explvar(mod)
    # Prediction, scores, loadings plots
    # plot(mod)
    # plot(mod, plottype = "scores")
    # plot(mod, "loadings", legendpos = "topleft")
    
  } else if(which.mod == "caretmod") {  
    
    library(caret)
    set.seed(111)
    mod <- train(score ~ ., data = plsdata, method = "pls", metric = "RMSE", tuneLength = 20) 
    #mod.new
    # 2 LVs were used for the final model
    #plot(mod.new)
    #plot(varImp(mod.new))
  }
  
}

lapply(list(Bioc1, Bioc2, FAdata, Bioc0), get.signature)

# Fit final PLS models to get signatures
set.seed(111)
mod1a <- plsr(score ~ ., data = Bioc1, ncomp = 1)
mod1b <- plsr(score ~ ., data = Bioc2, ncomp = 1)
mod2  <- plsr(score ~ ., data = FAdata, ncomp = 2)
mod0  <- plsr(score ~ ., data = Bioc0, ncomp = 1)

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