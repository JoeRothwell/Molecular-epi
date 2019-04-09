# Preparation of CRC case-control datasets, control only datasets and 
# calculation of WCRF-derived signatures with visualisation

# Read in small CRC case control dataset (1) and subset the 146 biocrates compounds
# stored at \\inti\NME\EPIC_Projects\Epic_Colonrectum\Nested_CaCo_Study\2016
# Missings are already imputed
# In CRC1, subjects with Biocrates data must first be subset

# Prep three datasets

prepcrc1 <- function(){

  library(tidyverse)
  library(haven)
  crc <- read_sas("clrt_caco_metabo.sas7bdat")
  
  # Metadata and WCRF scores (keep on local drive because big file)
  wcrf <- read_dta("Wcrf_Score.dta") %>% select(Idepic, Wcrf_C_Cal)
  
  # Remove duplicated Idepics (with dplyr or base)
  meta <- read_dta("clrt_caco.dta") %>% select(-Match_Caseset, -Country, -Center, -Cncr_Caco_Clrt) %>%
    distinct(Idepic, .keep_all = T)
  
  # meta.dedup <- meta[!duplicated(meta["Idepic"]), ]
  
  # Converted to dta
  # write_dta(data, "clrt_caco_metabo.dta")
  
  # Subset biocrates concentrations
  # concs <- crc %>% select(Acylcarn_C0 : Sugars_H1, -Batch_MetBio)
  concs <- crc %>% select(matches("Acylcarn_|Aminoacid_|Biogenic_|Glyceroph_|Sphingo_|Sugars_"), 
                           -starts_with("Outdq"))
  
  # Subset biocrates data by taking !is.na > 0
  biocrates <- apply(concs, 1, function(x) sum(!is.na(x)) > 0)
  
  # Subset Biocrates subjects only
  crc1 <- crc[biocrates, ]
  crc1 <- crc1 %>% inner_join(meta, by = "Idepic")
  
  var.list <- c("Country", "Center", "Sex", "Match_Caseset", "Smoke_Stat", "L_School")
  crc1 <- crc1 %>% mutate_at(vars(var.list), as.factor)
  
  #class(crc1$Sex)

}
crc1 <- prepcrc1()

prepcrc2 <- function(){
  library(haven)
  wcrf <- read_dta("Wcrf_Score.dta") %>% select(ends_with("_Cal"), Idepic)
  crc2 <- read_csv("biocrates_p150.csv")
  crc2 <- crc2 %>% left_join(wcrf, by = "Idepic")
  
  var.list <- c("Country", "Center", "Sex", "Match_Caseset", "Smoke_Stat", "L_School")
  crc2 <- crc2 %>% mutate_at(vars(var.list), as.factor)
}
crc2 <- prepcrc2()

prepctrl <- function(data = c("old", "new")){
  
  if(data == "new") {
  
    library(haven)
    ctrl <- read_dta("obes_metabo.dta")
    
    # Convert Study variable to factor and rename levels
    ctrl$Study <- fct_recode(ctrl$Study, Colonrectum_S = "Colonrectum", Colonrectum_L = "Colonrectum (bbmri)")
    
    # Split Batch_MetBio into two columns, then extract numeric variable
    # Remove 1694 CRC controls
    ctrl <- ctrl %>% separate(Batch_MetBio, into = c("batch", "rest")) %>%
      mutate(batch_no = as.numeric(flatten(str_extract_all(batch, "[0-9]+")))) %>% 
      filter(!(Study %in% c("Colonrectum_L", "Colonrectum_S")), Fasting_C == 2)
    # 1799 fasted subjects left
  
  } else if(data == "old") {
  
    # Old data (without new Prostate study)
    # Note: original file is now obes_meta1.dta
    ctrl <- readRDS("Biocrates data controls.rds")
    ctrl$Fasting_C[ctrl$Fasting_C == ""] <- NA
    #outliers <- c(3600, 3722)
    #ob <- ctrl[-outliers, ]
    ob <- ctrl
    
    # make new numeric variable for batch and exclude colorectal controls and non-fasted samples
    # needs to first convert Batch$MetBio to numeric
    ob <- ob %>% mutate(batch_no = as.numeric(flatten(str_extract_all(Batch_MetBio, "[0-9]+")))) 
    
    ctrl <- ctrl %>% separate(Batch_MetBio, into = c("batch", "rest")) %>%
    mutate(batch_no = as.numeric(flatten(str_extract_all(batch, "[0-9]+")))) %>% 
    filter(Fasting_C == "Yes", Study != "Colonrectum")
  
  }
  
}
ctrl <- prepctrl(data = "new")

# Calculate Biocrates and fatty acid signatures

signature1 <- function(fasting = T) {

  library(tidyverse)
  # First dataset, 3771 obs; updated November 2018 7191 obs
  print(paste(nrow(ctrl), "Controls read"))
  
  # Subset biocrates compounds
  controls <- ctrl %>% select(matches("Acylcarn_|Aminoacid_|Biogenic_|Glyceroph_|Sphingo_|Sugars_"), 
                            -starts_with("Outdq"))
  
  zerocols <- apply(controls, 2, function(x) sum(x, na.rm = T)) != 0
  controls <- controls[, zerocols]
  
  colnames(controls) %>% length # 147 variables
  
  # Large CRC metabolomics subset ----
  
  # (from Jelena, ~ 1200 case-control pairs), call it crc2
  # Join scores and filter out non-fasted samples
  crc2 <- if(fasting == T) crc2 %>% filter(Fasting_C == 2) else crc2
  print(paste("Subjects in larger case-control: ", nrow(crc2)))
  # Variables were converted to factors in CRC_data_prep.R
  
  #setA <- crc2 %>% select(Acylcarn_C0 : Sugars_H1, -starts_with("Outdq_"))
  
  setA <- crc2 %>% select(matches("Acylcarn_|Aminoacid_|Biogenic_|Glyceroph_|Sphingo_|Sugars_"), 
                          -starts_with("Outdq"))
  zerocols1 <- apply(setA, 2, function(x) sum(x, na.rm = T)) != 0
  setA <- setA[, zerocols1]
  colnames(setA) %>% length
  # 163 variables
  
  # Get common compounds between controls and CRC CC
  common_cols <- intersect(colnames(controls), colnames(setA))   # 135 compounds in common
  
  # Small CRC metabolomics subset----
  # (from Bertrand, ~ 490 case-control pairs), call it crc1
  print(paste("Subjects in smaller case-control:", nrow(crc1)))
  
  # setB <- crc1 %>% select(Acylcarn_C0 : Sugars_H1, -Batch_MetBio)
  setB <- crc1 %>% select(matches("Acylcarn_|Aminoacid_|Biogenic_|Glyceroph_|Sphingo_|Sugars_"), 
                           -starts_with("Outdq"))
  colnames(setB) %>% length
  
  # Convert variables to factors
  var.list <- c("Country", "Center", "Sex")
  crc1 <- crc1 %>% mutate_at(vars(var.list), as.factor)
  #common_cols <- intersect(colnames(controls), colnames(setB))   # 145 compounds in common
  
  # ----
  
  # Get common cols between all three datasets (controls, small CC, big CC)
  common_cols2 <- intersect(common_cols, colnames(setB))   # Glyceroph_Lysopc_A_C24_0 is removed
  
  # ----  
  
  # subset and reorder both datasets to get the same 126 compounds in the same order
  controls <- controls %>% select(one_of(common_cols2))
  setA <- setA %>% select(one_of(common_cols2))
  setB <- setB %>% select(one_of(common_cols2))
  
  # check colnames are the same for both sets
  identical(colnames(controls), colnames(setA))
  
  # Prepare controls matrix----
  # Replace zero, impute with half mins, scale
  concs <- as.matrix(controls)
  concs[concs == 0] <- NA
  library(zoo)
  concs1 <- na.aggregate(concs, FUN = function(x) min(x)/2)
  logconcs <- log2(concs1) %>% scale
  
  # adjust matrix for study, centre, sex, batch, BMI
  library(lme4)
  adj5 <- function(x) residuals(lmer(x ~ Center + batch_no + Sex + Bmi_C + (1|Study), data = ctrl))
  adjmat <- apply(logconcs, 2, adj5)
  
  # Subset calibrated scores
  score <- data_frame(score = ctrl$Wcrf_C_Cal)
  
  # Data setup. Must be a df with Bind scores to log matrix
  plsdata <- cbind(score, adjmat) %>% filter(!is.na(score))
  
  # PLS model----
  # For metabolic signature of WCRF score on controls dataset
  
  library(pls)
  set.seed(111)
  mod <- plsr(score ~ ., data = plsdata, validation = "CV")
  #summary(mod)
  
  # Find the number of dimensions with lowest cross validation error
  cv <- RMSEP(mod)
  plot(RMSEP(mod), legendpos = "topright")
  
  # Calculate optimal number of dimensions and rerun model
  best.dims <- which.min(cv$val[estimate = "adjCV", , ]) - 1
  mod <- plsr(score ~ ., data = plsdata, ncomp = best.dims)
  
  # explained variances
  # explvar(mod)
  
  # Plots ----
  # Prediction, scores, loadings
  # plot(mod)
  # plot(mod, plottype = "scores")
  # plot(mod, "loadings", legendpos = "topleft")
  
}
mod1 <- signature1()

signature2 <- function(){
  
  # Fatty acid signatures of WCRF score
  library(tidyverse)
  library(haven)
  
  # Get WCRF scores
  wcrf <- read_dta("Wcrf_Score.dta") %>% select(Idepic, Wcrf_C_Cal)
  
  # Get CRC dataset from Elom and join WCRF scores. Convert categorical co-variates to factors
  var.list <- c("L_School", "Smoke_Stat")
  CRCfa1 <- read_dta("Database_Fatty acids.dta") %>% left_join(wcrf, by = "Idepic") %>% mutate_at(vars(var.list), as.factor)
  
  # Subset FA concentrations
  CRCfa <- CRCfa1 %>% select(P14_0 : PCLA_9t_11c) 
  
  # Get dataset for PLS modelling (all EPIC controls). Exclude compounds with many missings
  # Note: new version from Carine received 18/11/2018 with technical covariates
  fa.scores <- readRDS("FA_WCRF_scores1.rds") %>% filter(STUDY != "Colorectum")
  fa.scores$N_Serie <- as.numeric(fa.scores$N_Serie)
  
  # convert categorical variables to factors
  
  var.list <- c("Country", "Center", "STUDY", "LABO")
  fa.scores <- fa.scores %>% mutate_at(vars(var.list), as.factor)
  
  concs <- fa.scores %>% select(P14_0 : PCLA_9t_11c, -P24_0, -P20_0)
  
  #library(Amelia)
  #missmap(concs)
  #missmap(CRCfa)
  
  length(colnames(concs))
  length(colnames(CRCfa))
  common_cols <- intersect(colnames(concs), colnames(CRCfa))
  
  # select common FAs in the same order
  CRCfa <- read_dta("Database_Fatty acids.dta") %>% select(one_of(common_cols)) 
  concs <- fa.scores %>% select(one_of(common_cols))
  
  identical(colnames(CRCfa), colnames(concs))
  
  # PLS to get signature of fatty acid metabolism for high and low scorers
  # Set zeros to NA and impute with half miminum conc, log transform
  # Bind to data frame
  
  concs <- as.matrix(concs)
  concs[concs == 0] <- NA
  library(zoo)
  concs <- na.aggregate(concs, FUN = function(x) min(x)/2)
  logconcs <- log2(concs) %>% scale
  
  # adjust matrix for study, centre, sex, batch, BMI
  library(lme4)
  #getresiduals <- function(x) residuals(lm(x ~ STUDY + LABO + Country, data = fa.scores))
  getresiduals <- function(x) residuals(lmer(x ~ (1|STUDY) + LABO + Country, data = fa.scores))
  resmat <- apply(logconcs, 2, getresiduals)
  
  # Bind WCRF scores to log2 concs
  plsdata <- data.frame(score = fa.scores$Wcrf_C_Cal, resmat) %>% filter(!is.na(score))
  
  library(pls)
  set.seed(111)
  mod <- plsr(score ~ ., data = plsdata, validation = "CV")
  # summary(mod)
  
  # Find the number of dimensions with lowest cross validation error
  cv <- RMSEP(mod)
  plot(RMSEP(mod), legendpos = "topright")
  
  # Calculate optimal number of dimensions
  best.dims <- which.min(cv$val[estimate = "adjCV", , ]) - 1
  
  # Rerun the model with the optimum number of components
  mod <- plsr(score ~ ., data = plsdata, ncomp = best.dims)
  
  # explained variances
  # explvar(mod)
  
  # Plots: prediction, scores, loadings
  #plot(mod)
  #plot(mod, plottype = "scores")
  #plot(mod, "loadings", legendpos = "topleft")
  
}
mod2 <- signature2()

# Produce tables of important compounds

plot_sig1 <- function(mod){
  
  # Coefficients and variable importance. First subset one-matrix array to get matrix
  k <- 126
  coeff <- data.frame(value = coef(mod)[1:k, 1, 1])
  dat <- coeff %>% mutate(sm  = sum(abs(value)), 
           VIP = (value*100)/sm, Compound = rownames(coeff))
  
  # Subset appropriate columns and print influential compounds
  dat <- dat %>% select(Compound, Coeff = value, VIP) %>% arrange(VIP)
  
  # choose number of influential compounds to plot
  n_top <- 7
  
  df1 <- top_n(dat, n_top)  %>% arrange(-VIP)
  df2 <- top_n(dat, -n_top) %>% arrange(VIP)
  
  print(df1)
  print(df2)
  
  # Vector of black and grey for plot points
  vec <- c( rep("black", n_top), rep("grey", k-(n_top*2)), rep("black", n_top) )
  
  # Now plot data, adding text
  plot(sort(coeff$value), pch = 17, col=vec, xlab = "", ylab = "Coefficient",
       main = paste(nrow(mod$scores), "fasted subjects, optimal dimensions =", mod$ncomp))
  # High labels
  text(nrow(dat) : (nrow(dat)-n_top), df1$Coeff, df1$Compound, pos=2, cex = 0.6)
  # Low labels
  text(1:nrow(df2), df2$Coeff, df2$Compound, pos=4, cex=0.6)
  abline(a=0, b=0, lty = "dotted")
  
  output <- bind_rows(df1, df2)
  
}
table3a <- plot_sig1(mod1)

plot_sig2 <- function(mod){
  
  # Plot the variable importance
  k <- 34
  coeff <- data.frame(value = coef(mod)[1:k, 1, 1])
  dat <- coeff %>% mutate(sm  = sum(abs(value)), 
                          VIP = (value*100)/sm, Compound = rownames(coeff))
  
  # Subset appropriate columns and print influential compounds
  dat <- dat %>% select(Compound, Coeff = value, VIP) %>% arrange(VIP)
  
  # choose number of influential compounds to plot
  n_top <- 5
  
  df1 <- top_n(dat, n_top)  %>% arrange(-VIP)
  df2 <- top_n(dat, -n_top) %>% arrange(VIP)
  
  print(df1)
  print(df2)
  
  # Vector of colours for plot points
  vec <- c( rep("black", n_top), rep("grey", k-(n_top*2)), rep("black", n_top) )
  
  # Now plot data, adding text
  plot(sort(coeff$value), pch = 17, col=vec, xlab = "", ylab = "Coefficient",
       main = paste("Fatty acids:", nrow(mod$scores), "fasted subjects, optimal dimensions =", mod$ncomp))
  # High compounds
  text(nrow(dat) : (nrow(dat)-n_top), df1$Coeff, df1$Compound, pos=2, cex = 0.6)
  # Low compounds
  text(1:nrow(df2), df2$Coeff, df2$Compound, pos=4, cex=0.6)
  abline(a=0, b=0, lty = "dotted")

  output <- bind_rows(df1, df2)
}
table3b <- plot_sig2(mod2)




