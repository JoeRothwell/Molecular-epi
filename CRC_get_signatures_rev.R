# Compute Biocrates and fatty acid signatures of WCRF score by PLS and predicts scores
source("CRC_prep_data_rev.R")
library(tidyverse)
library(lme4)
library(zoo)

# Get compound matrix and questionnaire scores
# Adjusts and scales the controls metabolite matrices for confounders and binds scores for model
# Define linear models for removal of confounders from matrix

### Revised manuscript for CGH september 2020: ctrlA now represents both studies A and B merged ###

adj   <- function(x) residuals(lmer(x ~ Center + batch_no + Sex + (1|Study), data = ctrl))
adjFA <- function(x) residuals(lmer(x ~ LABO + STUDY + (1|Center), data = fa.ctrl))

# Replace zero, impute with half mins, scale, transform to residuals of models above
hm <- function(x) min(x)/2
logmat1 <- ctrlA %>% na_if(0) %>% na.aggregate(FUN = hm) %>% log2 %>% scale
adjmat1 <- apply(logmat1, 2, adj) %>% data.frame # Overlapping Biocrates compounds

logmat3 <- ctrlC %>% na_if(0) %>% na.aggregate(FUN = hm) %>% log2 %>% scale
adjmat3 <- apply(logmat3, 2, adjFA) %>% data.frame # Fatty acids

# Bind WCRF scores to adjusted metabolite matrix for PLS modelling
Bioc1 <- cbind(score = ctrl$Wcrf_C_Cal, adjmat1) %>% filter(!is.na(score))
Facid <- cbind(score = fa.ctrl$Wcrf_C_Cal, adjmat3) %>% filter(!is.na(score))

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
mod1 <- plsr(score ~ ., data = Bioc1, ncomp = 1)
mod2  <- plsr(score ~ ., data = Facid, ncomp = 2)

# Make tables of important compounds, using compound metadata to get proper names
plot.sig <- function(mod, biocrates = T, percentile = 5, all = T){
  
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
table3a <- plot.sig(mod1, percentile = 10, all = F)
table3b <- plot.sig(mod2, biocrates = F, percentile = 10, all = F)

# Data for scatter plots
pltdata <- plot.sig(mod1, all = T)
faplot  <- plot.sig(mod2, biocrates = F, percentile = 10, all = T)

# Save workspace (for .Rmd file)
#save.image(file="metabolic_signatures.Rdata")

# Finally, predict WCRF scores from Biocrates or fatty acids data for two datasets (x3)

# Revision: all case-control subjects in one. crc2, colon2, rectal2 are replaced by crc, colon, rectal, etc
# Prep data. CRC, colon, rectal, by sex and excluding 2 year follow up
hm <- function(x) min(x)/2
cols <- colnames(ctrlA)
df2 <- crc[, cols] %>% na_if(0) %>% na.aggregate(FUN = hm) %>% log2 %>% scale %>% as_tibble
df3 <- colon[, cols] %>% na_if(0) %>% na.aggregate(FUN = hm) %>% log2 %>% scale %>% as_tibble
df4 <- rectal[, cols] %>% na_if(0) %>% na.aggregate(FUN = hm) %>% log2 %>% scale %>% as_tibble
df3a <- prox[, cols] %>% na_if(0) %>% na.aggregate(FUN = hm) %>% log2 %>% scale %>% as_tibble
df3b <- dist[, cols] %>% na_if(0) %>% na.aggregate(FUN = hm) %>% log2 %>% scale %>% as_tibble

# Excluding 2 year follow up, male and female only
dff <- crcF[, cols] %>% na_if(0) %>% na.aggregate(FUN = hm) %>% log2 %>% scale %>% as_tibble
dfm <- crcM[, cols] %>% na_if(0) %>% na.aggregate(FUN = hm) %>% log2 %>% scale %>% as_tibble
dft <- crcT[, cols] %>% na_if(0) %>% na.aggregate(FUN = hm) %>% log2 %>% scale %>% as_tibble

# High and low BMI, high and low WCRF scores
dfn <- crcN[, cols] %>% na_if(0) %>% na.aggregate(FUN = hm) %>% log2 %>% scale %>% as_tibble
dfo <- crcO[, cols] %>% na_if(0) %>% na.aggregate(FUN = hm) %>% log2 %>% scale %>% as_tibble
dfl <- crcL[, cols] %>% na_if(0) %>% na.aggregate(FUN = hm) %>% log2 %>% scale %>% as_tibble
dfh <- crcH[, cols] %>% na_if(0) %>% na.aggregate(FUN = hm) %>% log2 %>% scale %>% as_tibble

# Excluding 2 year follow up, male and female only
dft1 <- colonT[, cols] %>% na_if(0) %>% na.aggregate(FUN = hm) %>% log2 %>% scale %>% as_tibble
dff1 <- colonF[, cols] %>% na_if(0) %>% na.aggregate(FUN = hm) %>% log2 %>% scale %>% as_tibble
dfm1 <- colonM[, cols] %>% na_if(0) %>% na.aggregate(FUN = hm) %>% log2 %>% scale %>% as_tibble

# Excluding 2 year follow up, male and female only
dft2 <- rectalT[, cols] %>% na_if(0) %>% na.aggregate(FUN = hm) %>% log2 %>% scale %>% as_tibble
dff2 <- rectalF[, cols] %>% na_if(0) %>% na.aggregate(FUN = hm) %>% log2 %>% scale %>% as_tibble
dfm2 <- rectalM[, cols] %>% na_if(0) %>% na.aggregate(FUN = hm) %>% log2 %>% scale %>% as_tibble

# Proximal and distal exclude 2 year follow up and by sex
dft3 <- proxT[, cols] %>% na_if(0) %>% na.aggregate(FUN = hm) %>% log2 %>% scale %>% as_tibble
dff3 <- proxF[, cols] %>% na_if(0) %>% na.aggregate(FUN = hm) %>% log2 %>% scale %>% as_tibble
dfm3 <- proxM[, cols] %>% na_if(0) %>% na.aggregate(FUN = hm) %>% log2 %>% scale %>% as_tibble

dft4 <- distT[, cols] %>% na_if(0) %>% na.aggregate(FUN = hm) %>% log2 %>% scale %>% as_tibble
dff4 <- distF[, cols] %>% na_if(0) %>% na.aggregate(FUN = hm) %>% log2 %>% scale %>% as_tibble
dfm4 <- distM[, cols] %>% na_if(0) %>% na.aggregate(FUN = hm) %>% log2 %>% scale %>% as_tibble


# Predict and bind scores
# CRC, colon and rectal overall
crc.ph <- cbind(crc, comp1 = predict(mod1b, df2)[,,1]) %>% group_by(Match_Caseset) %>% filter(n() == 2)
col.ph <- cbind(colon, comp1 = predict(mod1b, df3)[,,1]) %>% group_by(Match_Caseset) %>% filter(n() == 2)
rec.ph <- cbind(rectal, comp1 = predict(mod1b, df4)[,,1]) %>% group_by(Match_Caseset) %>% filter(n() == 2)
prox.ph <- cbind(prox, comp1 = predict(mod1b, df3a)[,,1]) %>% group_by(Match_Caseset) %>% filter(n() == 2)
dist.ph <- cbind(dist, comp1 = predict(mod1b, df3b)[,,1]) %>% group_by(Match_Caseset) %>% filter(n() == 2)

# CRC by sex, excluding 2 year follow up
crcT.ph <- cbind(crcT, comp1 = predict(mod1b, dft)[,,1]) %>% group_by(Match_Caseset) %>% filter(n() == 2)
crcF.ph <- cbind(crcF, comp1 = predict(mod1b, dff)[,,1]) %>% group_by(Match_Caseset) %>% filter(n() == 2)
crcM.ph <- cbind(crcM, comp1 = predict(mod1b, dfm)[,,1]) %>% group_by(Match_Caseset) %>% filter(n() == 2)

# CRC by BMI and WCRF score strata (for paper revisions)
crcN.ph <- cbind(crcN, comp1 = predict(mod1b, dfn)[,,1]) %>% group_by(Match_Caseset) %>% filter(n() == 2)
crcO.ph <- cbind(crcO, comp1 = predict(mod1b, dfo)[,,1]) %>% group_by(Match_Caseset) %>% filter(n() == 2)
crcL.ph <- cbind(crcL, comp1 = predict(mod1b, dfl)[,,1]) %>% group_by(Match_Caseset) %>% filter(n() == 2)
crcH.ph <- cbind(crcH, comp1 = predict(mod1b, dfh)[,,1]) %>% group_by(Match_Caseset) %>% filter(n() == 2)


# Colon by sex and excluding 2 year follow up (new for revised manuscript)
colT.ph <- cbind(colonT, comp1 = predict(mod1b, dft1)[,,1]) %>% group_by(Match_Caseset) %>% filter(n() == 2)
colF.ph <- cbind(colonF, comp1 = predict(mod1b, dff1)[,,1]) %>% group_by(Match_Caseset) %>% filter(n() == 2)
colM.ph <- cbind(colonM, comp1 = predict(mod1b, dfm1)[,,1]) %>% group_by(Match_Caseset) %>% filter(n() == 2)

# Rectal by sex and excluding 2 year follow up (new for revised manuscript)
recT.ph <- cbind(rectalT, comp1 = predict(mod1b, dft2)[,,1]) %>% group_by(Match_Caseset) %>% filter(n() == 2)
recF.ph <- cbind(rectalF, comp1 = predict(mod1b, dff2)[,,1]) %>% group_by(Match_Caseset) %>% filter(n() == 2)
recM.ph <- cbind(rectalM, comp1 = predict(mod1b, dfm2)[,,1]) %>% group_by(Match_Caseset) %>% filter(n() == 2)

# Proximal and distal by sex and excluding 2 year follow up
proxT.ph <- cbind(proxT, comp1 = predict(mod1b, dft3)[,,1]) %>% group_by(Match_Caseset) %>% filter(n() == 2)
proxF.ph <- cbind(proxF, comp1 = predict(mod1b, dff3)[,,1]) %>% group_by(Match_Caseset) %>% filter(n() == 2)
proxM.ph <- cbind(proxM, comp1 = predict(mod1b, dfm3)[,,1]) %>% group_by(Match_Caseset) %>% filter(n() == 2)

distT.ph <- cbind(distT, comp1 = predict(mod1b, dft4)[,,1]) %>% group_by(Match_Caseset) %>% filter(n() == 2)
distF.ph <- cbind(distF, comp1 = predict(mod1b, dff4)[,,1]) %>% group_by(Match_Caseset) %>% filter(n() == 2)
distM.ph <- cbind(distM, comp1 = predict(mod1b, dfm4)[,,1]) %>% group_by(Match_Caseset) %>% filter(n() == 2)


# Small case-control fatty acids (Use predictions from 2 comps)
df5 <- crc3[, colnames(ctrlC)] %>% na_if(0) %>% na.aggregate(FUN = hm) %>% log2 %>% scale %>% as_tibble
df5t <- crc3t[, colnames(ctrlC)] %>% na_if(0) %>% na.aggregate(FUN = hm) %>% log2 %>% scale %>% as_tibble
df5f <- crc3f[, colnames(ctrlC)] %>% na_if(0) %>% na.aggregate(FUN = hm) %>% log2 %>% scale %>% as_tibble
df5m <- crc3m[, colnames(ctrlC)] %>% na_if(0) %>% na.aggregate(FUN = hm) %>% log2 %>% scale %>% as_tibble

df5n <- crc3N[, colnames(ctrlC)] %>% na_if(0) %>% na.aggregate(FUN = hm) %>% log2 %>% scale %>% as_tibble
df5o <- crc3O[, colnames(ctrlC)] %>% na_if(0) %>% na.aggregate(FUN = hm) %>% log2 %>% scale %>% as_tibble
df5l <- crc3L[, colnames(ctrlC)] %>% na_if(0) %>% na.aggregate(FUN = hm) %>% log2 %>% scale %>% as_tibble
df5h <- crc3H[, colnames(ctrlC)] %>% na_if(0) %>% na.aggregate(FUN = hm) %>% log2 %>% scale %>% as_tibble

df6a <- col3[, colnames(ctrlC)] %>% na_if(0) %>% na.aggregate(FUN = hm) %>% log2 %>% scale %>% as_tibble
df6b <- col3f[, colnames(ctrlC)] %>% na_if(0) %>% na.aggregate(FUN = hm) %>% log2 %>% scale %>% as_tibble
df6c <- col3m[, colnames(ctrlC)] %>% na_if(0) %>% na.aggregate(FUN = hm) %>% log2 %>% scale %>% as_tibble

# Proximal and distal (only 12 cases of rectal)
df7 <- prox3[, colnames(ctrlC)] %>% na_if(0) %>% na.aggregate(FUN = hm) %>% log2 %>% scale %>% as_tibble
df7a <- prox3f[, colnames(ctrlC)] %>% na_if(0) %>% na.aggregate(FUN = hm) %>% log2 %>% scale %>% as_tibble
df7b <- prox3m[, colnames(ctrlC)] %>% na_if(0) %>% na.aggregate(FUN = hm) %>% log2 %>% scale %>% as_tibble

df8 <- dist3[, colnames(ctrlC)] %>% na_if(0) %>% na.aggregate(FUN = hm) %>% log2 %>% scale %>% as_tibble
df8a <- dist3f[, colnames(ctrlC)] %>% na_if(0) %>% na.aggregate(FUN = hm) %>% log2 %>% scale %>% as_tibble
df8b <- dist3m[, colnames(ctrlC)] %>% na_if(0) %>% na.aggregate(FUN = hm) %>% log2 %>% scale %>% as_tibble

# By BMI and WCRF score strata (for paper revisions)
crc3N.ph <- cbind(crc3N, comp2 = predict(mod2, df5n)[,,1]) %>% group_by(Match_Caseset) %>% filter(n() == 2)
crc3O.ph <- cbind(crc3O, comp2 = predict(mod2, df5o)[,,1]) %>% group_by(Match_Caseset) %>% filter(n() == 2)
crc3L.ph <- cbind(crc3L, comp2 = predict(mod2, df5l)[,,1]) %>% group_by(Match_Caseset) %>% filter(n() == 2)
crc3H.ph <- cbind(crc3H, comp2 = predict(mod2, df5h)[,,1]) %>% group_by(Match_Caseset) %>% filter(n() == 2)

# Sex and follow-up time
crc3.ph <- cbind(crc3, comp2 = predict(mod2, df5)[,,2]) %>% group_by(Match_Caseset) %>% filter(n() == 2)
crc3t.ph <- cbind(crc3t, comp2 = predict(mod2, df5t)[,,2]) %>% group_by(Match_Caseset) %>% filter(n() == 2)
crc3f.ph <- cbind(crc3f, comp2 = predict(mod2, df5f)[,,2]) %>% group_by(Match_Caseset) %>% filter(n() == 2)
crc3m.ph <- cbind(crc3m, comp2 = predict(mod2, df5m)[,,2]) %>% group_by(Match_Caseset) %>% filter(n() == 2)

col3.ph <- cbind(col3, comp2 = predict(mod2, df6a)[,,2]) %>% group_by(Match_Caseset) %>% filter(n() == 2)
col3f.ph <- cbind(col3f, comp2 = predict(mod2, df6b)[,,2]) %>% group_by(Match_Caseset) %>% filter(n() == 2)
col3m.ph <- cbind(col3m, comp2 = predict(mod2, df6c)[,,2]) %>% group_by(Match_Caseset) %>% filter(n() == 2)

prox3.ph <- cbind(prox3, comp2 = predict(mod2, df7)[,,2]) %>% group_by(Match_Caseset) %>% filter(n() == 2)
prox3f.ph <- cbind(prox3f, comp2 = predict(mod2, df7a)[,,2]) %>% group_by(Match_Caseset) %>% filter(n() == 2)
prox3m.ph <- cbind(prox3m, comp2 = predict(mod2, df7b)[,,2]) %>% group_by(Match_Caseset) %>% filter(n() == 2)

dist3.ph <- cbind(dist3, comp2 = predict(mod2, df8)[,,2]) %>% group_by(Match_Caseset) %>% filter(n() == 2)
dist3f.ph <- cbind(dist3f, comp2 = predict(mod2, df8a)[,,2]) %>% group_by(Match_Caseset) %>% filter(n() == 2)
dist3m.ph <- cbind(dist3m, comp2 = predict(mod2, df8b)[,,2]) %>% group_by(Match_Caseset) %>% filter(n() == 2)


# Remove unneeded variables from workspace
rm(list = ls()[!str_detect(ls(), ".ph|.both")])
# Save workspace (for .Rmd file)
#save.image(file="pred_score_tables_rev.Rdata")
#save.image(file="pred_score_tables_FA.Rdata")
