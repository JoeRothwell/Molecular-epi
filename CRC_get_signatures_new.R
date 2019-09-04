# Compute Biocrates and fatty acid signatures of WCRF score by PLS
source("CRC_prep_data.R")

# Get compound matrix and questionnaire scores
get.plsdata <- function(dat, cor.data = F){
 
  # Adjusts and scales the controls metabolite matrices for confounders and binds scores for model
  library(tidyverse)
  library(lme4)
  library(zoo)
  # Set zeros to NA and impute with half miminum conc, log transform
  # Bind WCRF scores to adjusted metabolite matrix
  
  if(length(dat) == 2) {
  
  common.cols <- dat[[2]]
  fa.ctrl <- dat[[1]]
  
  # select common FAs in the same order. FAs is for correlation with Biocrates compounds
  FAs   <- read_dta("Database_Fatty acids.dta") %>% select(Idepic, one_of(common.cols)) 
  CRCfa <- FAs %>% select(-Idepic)
  concs <- fa.ctrl %>% select(one_of(common.cols))
  if(cor.data == T) return(FAs)
  
  # Fatty acids: adjust matrix for study, centre, sex, batch, BMI
  adj <- function(x) residuals(lmer(x ~ LABO + STUDY + (1|Center), data = fa.ctrl))
  # Subset calibrated scores
  score <- data_frame(score = fa.ctrl$Wcrf_C_Cal)
  
  } else {
    
  concs <- as.matrix(dat)
  # Biocrates: adjust matrix for study, centre, sex, batch, BMI
  adj <- function(x) residuals(lmer(x ~ Center + batch_no + Sex + Bmi_C + (1|Study), data = ctrl))
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

Bioc0  <- get.plsdata(ctrls) # All control compounds
Bioc1  <- get.plsdata(ctrlA) # Overlap control/A
Bioc2  <- get.plsdata(ctrlB) # Overlap control/B
Bioc3  <- get.plsdata(ctrls0) # Overlap control/A/B
FAdata <- get.plsdata(output) # Fatty acids
  
get.signature <- function(plsdata, which.mod = "plsmod"){
  
  # Biocrates (endogenous metabolites) for metabolic signature of WCRF score on controls dataset
  if(which.mod == "plsmod"){
    
    library(pls)
    
    # Start with a sensible number of components eg 10
    set.seed(111)
    mod <- plsr(score ~ ., ncomp = 10, data = plsdata, validation = "CV")
    # Find the number of dimensions with lowest cross validation error
    cv <- RMSEP(mod)
    plot(RMSEP(mod), legendpos = "topright")
    
    # Calculate optimal number of dimensions and rerun model
    # (See PLS vignette p12 for how to choose number of components)
    
    # Simply choosing which component has the lowest RMSEP
    best.dims <- which.min(cv$val[estimate = "adjCV", , ]) - 1
    
    # Choose using "one SE" and "permutation" methods
    ncomp.onesigma <- selectNcomp(mod, method = "onesigma", plot = F)
    ncomp.permut <- selectNcomp(mod, method = "randomization", plot = T)
    
    print(paste("Lowest RMSEP from", best.dims, "comp(s);", 
                "one SE method suggests", ncomp.onesigma, "comp(s);",
                "permutation method suggests", ncomp.permut, "comp(s)"))
    
    mod <- plsr(score ~ ., data = plsdata, ncomp = 2)
    
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

mod0  <- get.signature(Bioc0)
mod1a <- get.signature(Bioc1)
mod1b <- get.signature(Bioc2)
mod1c <- get.signature(Bioc3)
mod2  <- get.signature(FAdata)

#mod1a <- get.signature(which.mod = "caretmod")
#mod2a <- get.signature(which.mod = "caretmod")

# Produce tables of important compounds, using compound metadata to get proper names
plot.signature <- function(mod, biocrates = T, no.cmpds = 7, data.only = F){
  
  library(tidyverse)
  cmpds <- if(biocrates == T) read.csv("Biocrates_cmpd_metadata.csv") else read.csv("FA_compound_data.csv")
  
  # Coefficients and variable importance. First subset one-matrix array to get matrix
  # Extract coefficients from pls or caret objects, 1st LV: (pls gives an mvr object)
  if(class(mod) == "mvr"){
    coeff <- data.frame(value = round(coef(mod)[ , 1, 1], 3))
  } else {
    coeff <- data.frame(value = mod$finalModel$coefficients[, 1, 1])
  }
  
  dat <- coeff %>%
    rownames_to_column(var = "Compound") %>%
    mutate(sm = sum(abs(value))) %>%
    left_join(cmpds, by = "Compound") %>% 
    
    # Subset appropriate columns and print influential compounds
  select(class, compound = displayname, Coefficient = value)
  #select(class, compound = displayname2, Coefficient = value)
  
  if(data.only == T) return(dat)
  
  # choose number of influential compounds to plot
  n_top <- no.cmpds
  
  df1 <- top_n(dat, n_top)
  df2 <- top_n(dat, -n_top)
  print(df1)
  print(df2)
  
  # Vector of black and grey for plot points
  vec <- c( rep("black", n_top), rep("grey", nrow(dat) - (n_top*2)), rep("black", n_top) )
  
  lv <- if(class(mod) == "mvr") mod$ncomp else mod$finalModel$ncomp
  
  # Now plot data, adding text
  plot(sort(coeff$value), pch = 17, col=vec, xlab = "", ylab = "Coefficient",
       main = paste(nrow(mod$scores), "fasted subjects, optimal dimensions =", lv))
  # High and low labels
  text(nrow(dat) : (nrow(dat) - n_top), df1$Coefficient, df1$compound, pos=2, cex = 0.6)
  text(1:nrow(df2), df2$Coefficient, df2$compound, pos=4, cex=0.6)
  abline(a=0, b=0, lty = "dotted")
  
  output <- bind_rows(df1, df2)
}

table3a <- plot.signature(mod0)
table3b <- plot.signature(mod2, biocrates = F)

# Save workspace (for .Rmd file)
#save.image(file="metabolic_signatures.Rdata")

# Get data for coefficient plots
df1 <- plot.signature(mod1, data.only = T)
df2 <- plot.signature(mod2, biocrates = F, data.only = T)