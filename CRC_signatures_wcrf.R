# For revised submission to CGH. CRC1/A and CRC2/B are merged into one
library(tidyverse)
library(haven)
library(lubridate)

# Data prep----
# Read metadata for whole CRC case-control 
# WARNING: do not use case-control status from this dataset. Use metabolomics datasets only.
# Remove duplicated Idepics (with dplyr or base). Also get follow up time and colorectal site
var.list <- c("Country", "Center", "Sex", "Match_Caseset", "L_School", #"Smoke_Int", 
              "Smoke_Stat", "Smoke_Intensity", "Fasting_C", "Menopause", "Phase_Mnscycle")

# set D_Dgclrt of controls to that of corresponding cases, calculate followup time and 
# and get colorectal subsite variables
meta <- read_dta("clrt_caco.dta") %>% 
  mutate(Tfollowup.days = D_Dgclrt - D_Bld_Coll, Tfollowup = Tfollowup.days/365.25, 
         location = case_when(
           Case_Mal_Colon_Prox == 1 ~ 1, Case_Mal_Colon_Dist == 1 ~ 2,
           Case_Mal_Colon_Nos  == 1 ~ 4, Case_Mal_Rectum     == 1 ~ 3)) %>%
  group_by(Match_Caseset) %>% fill(c(D_Dgclrt, location), .direction = "downup") %>% ungroup() %>%
  select(-Match_Caseset, -Cncr_Caco_Clrt) %>%
  distinct(Idepic, .keep_all = T)

# Small case-control subset. Use glutamate to correctly subset biocrates data, join metadata. 
# Update 3/3/2020: remove Greece. Combine smoking intensity factor levels

crc1 <- read_sas("clrt_caco_metabo.sas7bdat") %>% filter(!is.na(Aminoacid_Glu)) %>%
  left_join(meta, by = "Idepic", suffix = c("_1", "")) %>%
  mutate_at(vars(var.list), as.factor) %>% 
  mutate(Smoke_Int = fct_collapse(Smoke_Intensity, Other = c("8", "9", "10"))) %>% 
  filter(Country != 6)

# Get colon cancer only (ungroup to stop Match_Caseset from being readded later)
colon1 <- crc1 %>% filter(location == 1 | location  == 2)

#colon1 <- crc1 %>% #group_by(Match_Caseset) %>% 
#  filter(max(location, na.rm = T) == 1 | max(location, na.rm = T) == 2) #%>% ungroup(Match_Caseset)

# Subsites
rectal1 <- crc1 %>% group_by(Match_Caseset) %>% filter(max(location, na.rm = T) == 3) %>% ungroup(Match_Caseset)
prox1 <- crc1 %>% group_by(Match_Caseset) %>% filter(max(location, na.rm = T) == 1) %>% ungroup(Match_Caseset)
dist1 <- crc1 %>% group_by(Match_Caseset) %>% filter(max(location, na.rm = T) == 2) %>% ungroup(Match_Caseset)

# By BMI and WCRF score?
bmiA <- crc1 %>% group_by(Match_Caseset) %>% filter(max(location, na.rm = T) == 2) %>% ungroup(Match_Caseset)
bmiB <- crc1 %>% group_by(Match_Caseset) %>% filter(max(location, na.rm = T) == 2) %>% ungroup(Match_Caseset)
bmiC <- crc1 %>% group_by(Match_Caseset) %>% filter(max(location, na.rm = T) == 2) %>% ungroup(Match_Caseset)

# Subset male, female, cases diagnosed after 2 years only
crc1m <- crc1 %>% filter(Sex == 1)
crc1f <- crc1 %>% filter(Sex == 2)
crc1t <- crc1 %>% group_by(Match_Caseset) %>% filter(max(Tfollowup, na.rm = T) > 2) %>% ungroup()

# Large case-control subset (from Jelena)
crc2 <- read_csv("biocrates_p150.csv") %>% select(Match_Caseset, Cncr_Caco_Clrt, ends_with("Idepic"), 
         matches("(carn|oacid|genic|roph|ingo|Sugars)[_]"), -contains("tdq")) %>%
  inner_join(meta, by = "Idepic") %>% mutate_at(vars(var.list), as.factor) %>% 
  mutate(Smoke_Int = fct_collapse(Smoke_Intensity, Other = c("8", "9", "10"))) %>%
  filter(Country != 6)

# Get colon cancer or rectal cancer subsets
colon2 <- crc2 %>% group_by(Match_Caseset) %>% 
  filter(max(location, na.rm = T) == 1 | max(location, na.rm = T) == 2) %>% ungroup(Match_Caseset)

rectal2 <- crc2 %>% group_by(Match_Caseset) %>% filter(max(location, na.rm = T) == 3) %>% ungroup(Match_Caseset)
prox2 <- crc2 %>% group_by(Match_Caseset) %>% filter(max(location, na.rm = T) == 1) %>% ungroup(Match_Caseset)
#prox2a <- crc2 %>% filter(location == 1)
dist2 <- crc2 %>% group_by(Match_Caseset) %>% filter( max(location, na.rm = T) == 2) %>% ungroup(Match_Caseset)


# Get cases diagnosed after 2 years only
crc2t <- crc2 %>% group_by(Match_Caseset) %>% filter(mean(Tfollowup, na.rm = T) > 2) %>% ungroup()

# Subset male or female
crc2m <- crc2 %>% filter(Sex == 1)
crc2f <- crc2 %>% filter(Sex == 2)


### 11 september 2020, reviewers revisions for CGH ###
# Merge crc1 and crc2 to make complete dataset
crc <- bind_rows(crc1, crc2, .id = "lab")
crc$lab <- as.factor(crc$lab)

colon <- bind_rows(colon1, colon2, .id = "lab")
colon$lab <- as.factor(colon$lab)

rectal <- bind_rows(rectal1, rectal2, .id = "lab")
rectal$lab <- as.factor(rectal$lab)

prox <- bind_rows(prox1, prox2, .id = "lab")
prox$lab <- as.factor(prox$lab)

dist <- bind_rows(dist1, dist2, .id = "lab")
dist$lab <- as.factor(dist$lab)

# CRC subsets
crcM <- crc %>% filter(Sex == 1)
crcF <- crc %>% filter(Sex == 2)
crcT <- crc %>% group_by(Match_Caseset) %>% filter(max(Tfollowup, na.rm = T) > 2) %>% ungroup()

# BMI and WCRF score: underweight/normal or overweight/obese; WCRF score 1-2 or 3-5
crcN <- crc %>% filter(Bmi_C >= 18.5 & Bmi_C < 25)
crcO <- crc %>% filter(Bmi_C >= 25)
crcL <- crc %>% filter(Wcrf_C_Cal %in% 1:2)
crcH <- crc %>% filter(Wcrf_C_Cal %in% 3:5)

colonM <- colon %>% filter(Sex == 1)
colonF <- colon %>% filter(Sex == 2)
colonT <- colon %>% group_by(Match_Caseset) %>% filter(max(Tfollowup, na.rm = T) > 2) %>% ungroup()

rectalM <- rectal %>% filter(Sex == 1)
rectalF <- rectal %>% filter(Sex == 2)
rectalT <- rectal %>% group_by(Match_Caseset) %>% filter(max(Tfollowup, na.rm = T) > 2) %>% ungroup()

proxT <- prox %>% group_by(Match_Caseset) %>% filter(max(Tfollowup, na.rm = T) > 2) %>% ungroup()
proxM <- prox %>% filter(Sex == 1)
proxF <- prox %>% filter(Sex == 2)

distT <- dist %>% group_by(Match_Caseset) %>% filter(max(Tfollowup, na.rm = T) > 2) %>% ungroup()
distM <- dist %>% filter(Sex == 1)
distF <- dist %>% filter(Sex == 2)


# EPIC pooled controls. First dataset, 3771 obs; updated November 2018 7191 obs
# Rename factor levels, split Batch_MetBio into 2 cols, extract numeric variable, Remove 1694 CRC controls

ctrl <- read_dta("obes_metabo.dta") %>% mutate(Study = 
      fct_recode(Study, Colonrectum_S = "Colonrectum", Colonrectum_L = "Colonrectum (bbmri)")) %>%
  separate(Batch_MetBio, into = c("batch", "rest")) %>%
  mutate(batch_no = as.numeric(flatten(str_extract_all(batch, "[0-9]+")))) %>% 
  filter(!(Study %in% c("Colonrectum_L", "Colonrectum_S")), Fasting_C == 2) %>%
  filter(Country != 6)

# Get controls dataset for crc1 and crc2 compounds only, and put in same order
# First subset compounds only from whole data and remove zero cols
expr <- "(carn|oacid|genic|roph|ingo|Sugars)[_]"
ctrls <- ctrl %>% select(matches(expr), -contains("tdq")) %>% select_if(~ sum(., na.rm = T) != 0)

# Overlap between whole case-control and discovery controls
crcp <- crc %>% select(matches(expr), -contains("tdq")) %>% select_if(~ sum(., na.rm = T) != 0)
ctrlA <- ctrls[, intersect(colnames(ctrls), colnames(crcp))]

# Remove compounds with over 40% missings
ctrlA <- ctrlA[, colSums(is.na(ctrlA)) < 697]
ctrlA <- ctrlA %>% select(-Sphingo_Sm_C26_1, -Sphingo_Sm_C26_0)

# Fatty acids CRC dataset (from Elom)

# Gets common compounds between case-control and discovery controls and puts them in the same order.

# Get CRC dataset from Elom and join WCRF scores. Convert categorical co-variates to factors
wcrf <- meta %>% select(Idepic, Wcrf_Fwg_Cal, Wcrf_Pf_Cal, Wcrf_Meat_Cal, Wcrf_C_Cal)

crc3 <- read_dta("Database_Fatty acids.dta") %>% 
  mutate(Tfollowup.days = D_Dgclrt - D_Bld_Coll, Tfollowup = Tfollowup.days/365.25, 
         location = case_when(
           Case_Mal_Colon_Prox == 1 ~ 1, Case_Mal_Colon_Dist == 1 ~ 2,
           Case_Mal_Colon_Nos  == 1 ~ 4, Case_Mal_Rectum     == 1 ~ 3)) %>%
  group_by(Match_Caseset) %>% fill(c(Tfollowup, location), .direction = "downup") %>% ungroup() %>%
  mutate(Smoke_Int = fct_collapse(as.factor(Smoke_Intensity), Other = c("8", "9", "10"))) %>%
  left_join(wcrf, by = "Idepic") %>%
  filter(Country != 6)

# Get male and female subsets, diagnosed after 2 years
crc3m <- crc3 %>% filter(Sex == 1)
crc3f <- crc3 %>% filter(Sex == 2)
crc3t <- crc3 %>% filter(Tfollowup > 2)

# Colon cancer only
col3 <- crc3 %>% filter(location %in% 1:2)
prox3 <- crc3 %>% filter(location == 1)
dist3 <- crc3 %>% filter(location == 2)
#rect3 <- crc3 %>% filter(location == 3) # only 5 cases

# Get high and low BMI and WCRF scores (for revision)
crc3N <- crc3 %>% filter(Bmi_C >= 18.5 & Bmi_C < 25)
crc3O <- crc3 %>% filter(Bmi_C >= 25)
crc3L <- crc3 %>% filter(Wcrf_C_Cal %in% 1:2)
crc3H <- crc3 %>% filter(Wcrf_C_Cal %in% 3:5)

# Colon proximal, distal, rectal cases only
col3 <- crc3 %>% filter(location %in% 1:2)
col3m <- col3 %>% filter(Sex == 1)
col3f <- col3 %>% filter(Sex == 2)

prox3 <- crc3 %>% filter(location == 1)
prox3m <- prox3 %>% filter(Sex == 1)
prox3f <- prox3 %>% filter(Sex == 2)

dist3 <- crc3 %>% filter(location == 2)
dist3m <- dist3 %>% filter(Sex == 1)
dist3f <- dist3 %>% filter(Sex == 2)

#rect3 <- crc3 %>% filter(location == 3) # only 5 cases


# Get dataset for PLS modelling (all EPIC controls). Exclude compounds with many missings
# Note: new version from Carine received 18/11/2018 with technical covariates
fa.ctrl <- readRDS("FA_WCRF_scores1.rds") %>% filter(STUDY != "Colorectum") %>% filter(Country != 6)
fa.ctrl$N_Serie <- as.numeric(fa.ctrl$N_Serie)

# Get sex of participants (not needed, only women)
#epic.vars <- read.csv("full_epic_sex.csv")
#fa.ctrl <- left_join(fa.ctrl, epic.vars, by = "Idepic")

# Convert variables to factors
var.list <- c("Country", "Center", "STUDY", "LABO")
fa.ctrl <- fa.ctrl %>% mutate_at(vars(var.list), as.factor)

# Subset concentrations for CRC and controls
crcfa <- crc3 %>% select(P14_0 : PCLA_9t_11c) 
concs <- fa.ctrl %>% select(P14_0 : PCLA_9t_11c, -P24_0, -P20_0)

# Compounds with > 20% CV
concs <- concs %>% select(-P14_1n_5, -P18_1n_7t, -P22_0, -P18_2n_6tt)


ctrlC <- fa.ctrl[, intersect(colnames(concs), colnames(crcfa))]

# Number of control profiles for biocrates and fatty acids
nrow(ctrl)
nrow(fa.ctrl)

# Number of control subjects , biocrates and fatty acids combined
length(intersect(ctrl$Idepic, fa.ctrl$Idepic))
#intersect(ctrl$Idepic, fa.ctrl$Idepic)

rm(crc1)
rm(crc2)
rm(colon1)
rm(colon2)
rm(rectal1)
rm(rectal2)
rm(crc1f)
rm(crc1m)
rm(crc2f)
rm(crc2m)
rm(crc1t)
rm(crc2t)

# Get signatures ----
# Previously CRC_get_signatures_rev.R
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
mod2 <- plsr(score ~ ., data = Facid, ncomp = 2)

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
faplot2  <- plot.sig(mod2, biocrates = F, percentile = 10, all = T)

# Save workspace (for .Rmd file)
#save.image(file="metabolic_signatures.Rdata")

# Finally, predict WCRF scores from Biocrates or fatty acids data for two datasets (x3)

# Revision: all case-control subjects in one. crc2, colon2, rectal2 are replaced by crc, colon, rectal, etc
# Prep data. CRC, colon, rectal, by sex and excluding 2 year follow up
hm <- function(x) min(x)/2
cols <- colnames(ctrlA)

# CRC
df2 <- crc[, cols] %>% na_if(0) %>% na.aggregate(FUN = hm) %>% log2 %>% scale %>% as_tibble

# Excluding 2 year follow up, male and female only
dff <- crcF[, cols] %>% na_if(0) %>% na.aggregate(FUN = hm) %>% log2 %>% scale %>% as_tibble
dfm <- crcM[, cols] %>% na_if(0) %>% na.aggregate(FUN = hm) %>% log2 %>% scale %>% as_tibble
dft <- crcT[, cols] %>% na_if(0) %>% na.aggregate(FUN = hm) %>% log2 %>% scale %>% as_tibble

# High and low BMI, high and low WCRF scores
dfn <- crcN[, cols] %>% na_if(0) %>% na.aggregate(FUN = hm) %>% log2 %>% scale %>% as_tibble
dfo <- crcO[, cols] %>% na_if(0) %>% na.aggregate(FUN = hm) %>% log2 %>% scale %>% as_tibble
dfl <- crcL[, cols] %>% na_if(0) %>% na.aggregate(FUN = hm) %>% log2 %>% scale %>% as_tibble
dfh <- crcH[, cols] %>% na_if(0) %>% na.aggregate(FUN = hm) %>% log2 %>% scale %>% as_tibble

# Fatty acids
df5 <- crc3[, colnames(ctrlC)] %>% na_if(0) %>% na.aggregate(FUN = hm) %>% log2 %>% scale %>% as_tibble
df5t <- crc3t[, colnames(ctrlC)] %>% na_if(0) %>% na.aggregate(FUN = hm) %>% log2 %>% scale %>% as_tibble
df5f <- crc3f[, colnames(ctrlC)] %>% na_if(0) %>% na.aggregate(FUN = hm) %>% log2 %>% scale %>% as_tibble
df5m <- crc3m[, colnames(ctrlC)] %>% na_if(0) %>% na.aggregate(FUN = hm) %>% log2 %>% scale %>% as_tibble

df5n <- crc3N[, colnames(ctrlC)] %>% na_if(0) %>% na.aggregate(FUN = hm) %>% log2 %>% scale %>% as_tibble
df5o <- crc3O[, colnames(ctrlC)] %>% na_if(0) %>% na.aggregate(FUN = hm) %>% log2 %>% scale %>% as_tibble
df5l <- crc3L[, colnames(ctrlC)] %>% na_if(0) %>% na.aggregate(FUN = hm) %>% log2 %>% scale %>% as_tibble
df5h <- crc3H[, colnames(ctrlC)] %>% na_if(0) %>% na.aggregate(FUN = hm) %>% log2 %>% scale %>% as_tibble


# Predict and bind scores
# CRC, colon and rectal overall
crc.ph <- cbind(crc, comp1 = predict(mod1, df2)[,,1]) %>% group_by(Match_Caseset) %>% filter(n() == 2)
col.ph <- cbind(colon, comp1 = predict(mod1, df3)[,,1]) %>% group_by(Match_Caseset) %>% filter(n() == 2)
rec.ph <- cbind(rectal, comp1 = predict(mod1, df4)[,,1]) %>% group_by(Match_Caseset) %>% filter(n() == 2)
prox.ph <- cbind(prox, comp1 = predict(mod1, df3a)[,,1]) %>% group_by(Match_Caseset) %>% filter(n() == 2)
dist.ph <- cbind(dist, comp1 = predict(mod1, df3b)[,,1]) %>% group_by(Match_Caseset) %>% filter(n() == 2)

# CRC by sex, excluding 2 year follow up
crcT.ph <- cbind(crcT, comp1 = predict(mod1, dft)[,,1]) %>% group_by(Match_Caseset) %>% filter(n() == 2)
crcF.ph <- cbind(crcF, comp1 = predict(mod1, dff)[,,1]) %>% group_by(Match_Caseset) %>% filter(n() == 2)
crcM.ph <- cbind(crcM, comp1 = predict(mod1, dfm)[,,1]) %>% group_by(Match_Caseset) %>% filter(n() == 2)

# CRC by BMI and WCRF score strata (for paper revisions)
crcN.ph <- cbind(crcN, comp1 = predict(mod1, dfn)[,,1]) %>% group_by(Match_Caseset) %>% filter(n() == 2)
crcO.ph <- cbind(crcO, comp1 = predict(mod1, dfo)[,,1]) %>% group_by(Match_Caseset) %>% filter(n() == 2)
crcL.ph <- cbind(crcL, comp1 = predict(mod1, dfl)[,,1]) %>% group_by(Match_Caseset) %>% filter(n() == 2)
crcH.ph <- cbind(crcH, comp1 = predict(mod1, dfh)[,,1]) %>% group_by(Match_Caseset) %>% filter(n() == 2)

crc3N.ph <- cbind(crc3N, comp2 = predict(mod2, df5n)[,,1]) %>% group_by(Match_Caseset) %>% filter(n() == 2)
crc3O.ph <- cbind(crc3O, comp2 = predict(mod2, df5o)[,,1]) %>% group_by(Match_Caseset) %>% filter(n() == 2)
crc3L.ph <- cbind(crc3L, comp2 = predict(mod2, df5l)[,,1]) %>% group_by(Match_Caseset) %>% filter(n() == 2)
crc3H.ph <- cbind(crc3H, comp2 = predict(mod2, df5h)[,,1]) %>% group_by(Match_Caseset) %>% filter(n() == 2)

crc3.ph <- cbind(crc3, comp2 = predict(mod2, df5)[,,2]) %>% group_by(Match_Caseset) %>% filter(n() == 2)
crc3t.ph <- cbind(crc3t, comp2 = predict(mod2, df5t)[,,2]) %>% group_by(Match_Caseset) %>% filter(n() == 2)
crc3f.ph <- cbind(crc3f, comp2 = predict(mod2, df5f)[,,2]) %>% group_by(Match_Caseset) %>% filter(n() == 2)
crc3m.ph <- cbind(crc3m, comp2 = predict(mod2, df5m)[,,2]) %>% group_by(Match_Caseset) %>% filter(n() == 2)


# Colon
df3 <- colon[, cols] %>% na_if(0) %>% na.aggregate(FUN = hm) %>% log2 %>% scale %>% as_tibble

# Excluding 2 year follow up, male and female only
dft1 <- colonT[, cols] %>% na_if(0) %>% na.aggregate(FUN = hm) %>% log2 %>% scale %>% as_tibble
dff1 <- colonF[, cols] %>% na_if(0) %>% na.aggregate(FUN = hm) %>% log2 %>% scale %>% as_tibble
dfm1 <- colonM[, cols] %>% na_if(0) %>% na.aggregate(FUN = hm) %>% log2 %>% scale %>% as_tibble

df6a <- col3[, colnames(ctrlC)] %>% na_if(0) %>% na.aggregate(FUN = hm) %>% log2 %>% scale %>% as_tibble
df6b <- col3f[, colnames(ctrlC)] %>% na_if(0) %>% na.aggregate(FUN = hm) %>% log2 %>% scale %>% as_tibble
df6c <- col3m[, colnames(ctrlC)] %>% na_if(0) %>% na.aggregate(FUN = hm) %>% log2 %>% scale %>% as_tibble

col.ph <- cbind(colon, comp1 = predict(mod1b, df3)[,,1]) %>% group_by(Match_Caseset) %>% filter(n() == 2)

# Colon by sex and excluding 2 year follow up (new for revised manuscript)
colT.ph <- cbind(colonT, comp1 = predict(mod1, dft1)[,,1]) %>% group_by(Match_Caseset) %>% filter(n() == 2)
colF.ph <- cbind(colonF, comp1 = predict(mod1, dff1)[,,1]) %>% group_by(Match_Caseset) %>% filter(n() == 2)
colM.ph <- cbind(colonM, comp1 = predict(mod1, dfm1)[,,1]) %>% group_by(Match_Caseset) %>% filter(n() == 2)

col3.ph <- cbind(col3, comp2 = predict(mod2, df6a)[,,2]) %>% group_by(Match_Caseset) %>% filter(n() == 2)
col3f.ph <- cbind(col3f, comp2 = predict(mod2, df6b)[,,2]) %>% group_by(Match_Caseset) %>% filter(n() == 2)
col3m.ph <- cbind(col3m, comp2 = predict(mod2, df6c)[,,2]) %>% group_by(Match_Caseset) %>% filter(n() == 2)


# Colon proximal
df3a <- prox[, cols] %>% na_if(0) %>% na.aggregate(FUN = hm) %>% log2 %>% scale %>% as_tibble

# Proximal exclude 2 year follow up and by sex
dft3 <- proxT[, cols] %>% na_if(0) %>% na.aggregate(FUN = hm) %>% log2 %>% scale %>% as_tibble
dff3 <- proxF[, cols] %>% na_if(0) %>% na.aggregate(FUN = hm) %>% log2 %>% scale %>% as_tibble
dfm3 <- proxM[, cols] %>% na_if(0) %>% na.aggregate(FUN = hm) %>% log2 %>% scale %>% as_tibble

df7 <- prox3[, colnames(ctrlC)] %>% na_if(0) %>% na.aggregate(FUN = hm) %>% log2 %>% scale %>% as_tibble
df7a <- prox3f[, colnames(ctrlC)] %>% na_if(0) %>% na.aggregate(FUN = hm) %>% log2 %>% scale %>% as_tibble
df7b <- prox3m[, colnames(ctrlC)] %>% na_if(0) %>% na.aggregate(FUN = hm) %>% log2 %>% scale %>% as_tibble

# Proximal and distal by sex and excluding 2 year follow up
proxT.ph <- cbind(proxT, comp1 = predict(mod1, dft3)[,,1]) %>% group_by(Match_Caseset) %>% filter(n() == 2)
proxF.ph <- cbind(proxF, comp1 = predict(mod1, dff3)[,,1]) %>% group_by(Match_Caseset) %>% filter(n() == 2)
proxM.ph <- cbind(proxM, comp1 = predict(mod1, dfm3)[,,1]) %>% group_by(Match_Caseset) %>% filter(n() == 2)

distT.ph <- cbind(distT, comp1 = predict(mod1, dft4)[,,1]) %>% group_by(Match_Caseset) %>% filter(n() == 2)
distF.ph <- cbind(distF, comp1 = predict(mod1, dff4)[,,1]) %>% group_by(Match_Caseset) %>% filter(n() == 2)
distM.ph <- cbind(distM, comp1 = predict(mod1, dfm4)[,,1]) %>% group_by(Match_Caseset) %>% filter(n() == 2)


# Small case-control fatty acids (Use predictions from 2 comps)
df5 <- crc3[, colnames(ctrlC)] %>% na_if(0) %>% na.aggregate(FUN = hm) %>% log2 %>% scale %>% as_tibble
df5t <- crc3t[, colnames(ctrlC)] %>% na_if(0) %>% na.aggregate(FUN = hm) %>% log2 %>% scale %>% as_tibble
df5f <- crc3f[, colnames(ctrlC)] %>% na_if(0) %>% na.aggregate(FUN = hm) %>% log2 %>% scale %>% as_tibble
df5m <- crc3m[, colnames(ctrlC)] %>% na_if(0) %>% na.aggregate(FUN = hm) %>% log2 %>% scale %>% as_tibble

prox3.ph <- cbind(prox3, comp2 = predict(mod2, df7)[,,2]) %>% group_by(Match_Caseset) %>% filter(n() == 2)
prox3f.ph <- cbind(prox3f, comp2 = predict(mod2, df7a)[,,2]) %>% group_by(Match_Caseset) %>% filter(n() == 2)
prox3m.ph <- cbind(prox3m, comp2 = predict(mod2, df7b)[,,2]) %>% group_by(Match_Caseset) %>% filter(n() == 2)


# Colon distal
df3b <- dist[, cols] %>% na_if(0) %>% na.aggregate(FUN = hm) %>% log2 %>% scale %>% as_tibble

dft4 <- distT[, cols] %>% na_if(0) %>% na.aggregate(FUN = hm) %>% log2 %>% scale %>% as_tibble
dff4 <- distF[, cols] %>% na_if(0) %>% na.aggregate(FUN = hm) %>% log2 %>% scale %>% as_tibble
dfm4 <- distM[, cols] %>% na_if(0) %>% na.aggregate(FUN = hm) %>% log2 %>% scale %>% as_tibble

df8 <- dist3[, colnames(ctrlC)] %>% na_if(0) %>% na.aggregate(FUN = hm) %>% log2 %>% scale %>% as_tibble
df8a <- dist3f[, colnames(ctrlC)] %>% na_if(0) %>% na.aggregate(FUN = hm) %>% log2 %>% scale %>% as_tibble
df8b <- dist3m[, colnames(ctrlC)] %>% na_if(0) %>% na.aggregate(FUN = hm) %>% log2 %>% scale %>% as_tibble

dist.ph <- cbind(dist, comp1 = predict(mod1b, df3b)[,,1]) %>% group_by(Match_Caseset) %>% filter(n() == 2)

distT.ph <- cbind(distT, comp1 = predict(mod1b, dft4)[,,1]) %>% group_by(Match_Caseset) %>% filter(n() == 2)
distF.ph <- cbind(distF, comp1 = predict(mod1b, dff4)[,,1]) %>% group_by(Match_Caseset) %>% filter(n() == 2)
distM.ph <- cbind(distM, comp1 = predict(mod1b, dfm4)[,,1]) %>% group_by(Match_Caseset) %>% filter(n() == 2)


dist3.ph <- cbind(dist3, comp2 = predict(mod2, df8)[,,2]) %>% group_by(Match_Caseset) %>% filter(n() == 2)
dist3f.ph <- cbind(dist3f, comp2 = predict(mod2, df8a)[,,2]) %>% group_by(Match_Caseset) %>% filter(n() == 2)
dist3m.ph <- cbind(dist3m, comp2 = predict(mod2, df8b)[,,2]) %>% group_by(Match_Caseset) %>% filter(n() == 2)


# Rectal
df4 <- rectal[, cols] %>% na_if(0) %>% na.aggregate(FUN = hm) %>% log2 %>% scale %>% as_tibble

# Excluding 2 year follow up, male and female only
dft2 <- rectalT[, cols] %>% na_if(0) %>% na.aggregate(FUN = hm) %>% log2 %>% scale %>% as_tibble
dff2 <- rectalF[, cols] %>% na_if(0) %>% na.aggregate(FUN = hm) %>% log2 %>% scale %>% as_tibble
dfm2 <- rectalM[, cols] %>% na_if(0) %>% na.aggregate(FUN = hm) %>% log2 %>% scale %>% as_tibble

rec.ph <- cbind(rectal, comp1 = predict(mod1b, df4)[,,1]) %>% group_by(Match_Caseset) %>% filter(n() == 2)

# Rectal by sex and excluding 2 year follow up (new for revised manuscript)
recT.ph <- cbind(rectalT, comp1 = predict(mod1b, dft2)[,,1]) %>% group_by(Match_Caseset) %>% filter(n() == 2)
recF.ph <- cbind(rectalF, comp1 = predict(mod1b, dff2)[,,1]) %>% group_by(Match_Caseset) %>% filter(n() == 2)
recM.ph <- cbind(rectalM, comp1 = predict(mod1b, dfm2)[,,1]) %>% group_by(Match_Caseset) %>% filter(n() == 2)




# Remove unneeded variables from workspace
rm(list = ls()[!str_detect(ls(), ".ph|.both")])
# Save workspace (for .Rmd file)
#save.image(file="pred_score_tables_rev.Rdata")
#save.image(file="pred_score_tables_FA.Rdata")
save.image(file="predscore_df_subsite1.Rdata")


# Model CRC status from WCRF score or signature
# Requires crc1, crc2, controls and PLS models to be prepared from CRC_data_prep.R
#source("CRC_get_signatures_rev.R")
load("predscore_df_subsite.Rdata")
library(tidyverse)

# For smoke intensity, categories 8, 9 and 10 are collapsed into other (Smoke_Int)
# Define basic model. Note: need to remove variables for rectal subset

### Revised submission to CGH: crc1 and crc2 now merged ###
# Models for score, LH column of table (fatty acids unchanged for resubmission)
# Models for WCRF score LH column of table (corresponding subsets)

library(survival)
base <- Cncr_Caco_Clrt ~ Qe_Energy + L_School + Smoke_Stat + Smoke_Int + 
  Height_C + strata(Match_Caseset)

fit1 <- clogit(update(base, ~. + Wcrf_C_Cal), data = crc3.ph)
fit3 <- clogit(update(base, ~. + Wcrf_C_Cal), data = crc3f.ph)
fit5 <- clogit(update(base, ~. + Wcrf_C_Cal), data = crc3m.ph)
fit7 <- clogit(update(base, ~. + Wcrf_C_Cal), data = crc.ph)
fit9 <- clogit(update(base, ~. + Wcrf_C_Cal), data = crcF.ph)
fit11 <- clogit(update(base, ~. + Wcrf_C_Cal), data = crcM.ph)
fit13 <- clogit(update(base, ~. + Wcrf_C_Cal), data = col.ph)
fit15 <- clogit(update(base, ~. + Wcrf_C_Cal), data = colF.ph)
fit17 <- clogit(update(base, ~. + Wcrf_C_Cal), data = colM.ph)
fit19 <- clogit(update(base, ~. + Wcrf_C_Cal), data = rec.ph)
fit21 <- clogit(update(base, ~. - Smoke_Stat + Wcrf_C_Cal), data = recF.ph)
fit23 <- clogit(update(base, ~. - L_School + Wcrf_C_Cal), data = recM.ph)

# Models for signature, RH column of table (fatty acids unchanged for resubmission)
fit2 <- clogit(update(base, ~. + comp2), data = crc3.ph)
fit4 <- clogit(update(base, ~. + comp2), data = crc3f.ph)
fit6 <- clogit(update(base, ~. + comp2), data = crc3m.ph)
fit8 <- clogit(update(base, ~. + lab + comp1), data = crc.ph)
fit10 <- clogit(update(base, ~. + lab + comp1), data = crcF.ph)
fit12 <- clogit(update(base, ~. + lab + comp1), data = crcM.ph)
fit14 <- clogit(update(base, ~. + lab + comp1), data = col.ph)
fit16 <- clogit(update(base, ~. + lab + comp1), data = colF.ph)
fit18 <- clogit(update(base, ~. + lab + comp1), data = colM.ph)
fit20 <- clogit(update(base, ~. + lab + comp1), data = rec.ph)
fit22 <- clogit(update(base, ~. + lab - Smoke_Stat + comp1), data = recF.ph)
fit24 <- clogit(update(base, ~. + lab - Smoke_Stat - L_School + comp1), data = recM.ph)

# Signature models ----
# Formerly CRC_signature_models_rev.R

# Put score and signature models in separate lists (separate columns in table)
scomodlist <- list(fit1, fit3, fit5, fit7, fit9, fit11, fit13, fit15, fit17, fit19, fit21, fit23)
sigmodlist <- list(fit2, fit4, fit6, fit8, fit10, fit12, fit14, fit16, fit18, fit20, fit22, fit24)

scorenames <- c("CRC FA, all", "CRC FA, women", "CRC FA men", 
                "CRC endogenous all", "CRC endogenous women", "CRC endogenous men", 
                "Colon endogenous all", "Colon endogenous women", "Colon endogenous men", 
                "Rectal endogenous, all", "Rectal endogenous, women", "Rectal endogenous, men")

library(broom)
scomods <- map_df(scomodlist, ~tidy(., exponentiate = T)) %>% filter(term == "Wcrf_C_Cal") %>%
  mutate_if(is.numeric, ~round(., 2)) %>% unite(OR.CI, estimate, conf.low, conf.high, sep = "-") %>%
  add_column(model = scorenames, .before = T)

sigmods <- map_df(sigmodlist, ~tidy(., exponentiate = T)) %>% filter(term == "comp2" | term == "comp1") %>%
  mutate_if(is.numeric, ~round(., 2)) %>% unite(OR.CI, estimate, conf.low, conf.high, sep = "-")# %>%
#add_column(model = scorenames[-c(10, 13)], .before = T)


# Additional models for revision: colon proximal and distal by sex, colon and distal for fatty acids

fit1 <- clogit(update(base, ~. + comp2), data = col3.ph)
fit2 <- clogit(update(base, ~. + comp2), data = col3f.ph)
fit3 <- clogit(update(base, ~. + comp2), data = col3m.ph)
fit4 <- clogit(update(base, ~. + comp2), data = prox3.ph)
fit5 <- clogit(update(base, ~. + comp2), data = prox3f.ph)
fit6 <- clogit(update(base, ~. + comp2), data = prox3m.ph)
fit7 <- clogit(update(base, ~. + comp2), data = dist3.ph)
fit8 <- clogit(update(base, ~. + comp2), data = dist3f.ph)
fit9 <- clogit(update(base, ~. + comp2), data = dist3m.ph)

fit10 <- clogit(update(base, ~. + lab + comp1), data = prox.ph)
fit12 <- clogit(update(base, ~. + lab + comp1), data = proxF.ph)
fit13 <- clogit(update(base, ~. + lab + comp1), data = proxM.ph)
fit14 <- clogit(update(base, ~. + lab + comp1), data = dist.ph)
fit15 <- clogit(update(base, ~. + lab + comp1), data = distF.ph)
fit16 <- clogit(update(base, ~. + lab + comp1), data = distM.ph)

modlist <- list(fit1, fit2, fit3, fit4, fit5, fit6, fit7, fit8, fit9)
modlist <- list(fit1, fit2, fit3, fit4, fit5, fit6)
modlist <- list(fit10, fit12, fit13, fit14, fit15, fit16)

library(broom)
mods <- map_df(modlist, ~tidy(., exponentiate = T)) %>% filter(term == "comp2" | term == "comp1") %>%
  mutate_if(is.numeric, ~round(., 2)) %>% unite(OR.CI, estimate, conf.low, conf.high, sep = "-")


# Score, fatty acids
fit7 <- clogit(update(base, ~. + Wcrf_C_Cal), data = col3.ph)
fit8 <- clogit(update(base, ~. + Wcrf_C_Cal), data = col3f.ph)
fit9 <- clogit(update(base, ~. + Wcrf_C_Cal), data = col3m.ph)

fit10 <- clogit(update(base, ~. + Wcrf_C_Cal), data = prox3.ph)
fit11 <- clogit(update(base, ~. + Wcrf_C_Cal), data = prox3f.ph)
fit12 <- clogit(update(base, ~. + Wcrf_C_Cal), data = prox3m.ph)

fit13 <- clogit(update(base, ~. + Wcrf_C_Cal), data = dist3.ph)
fit14 <- clogit(update(base, ~. + Wcrf_C_Cal), data = dist3f.ph)
fit15 <- clogit(update(base, ~. + Wcrf_C_Cal), data = dist3m.ph)

# Score, endogenous
fit1 <- clogit(update(base, ~. + Wcrf_C_Cal), data = prox.ph)
fit2 <- clogit(update(base, ~. + Wcrf_C_Cal), data = proxF.ph)
fit3 <- clogit(update(base, ~. + Wcrf_C_Cal), data = proxM.ph)

fit4 <- clogit(update(base, ~. + Wcrf_C_Cal), data = dist.ph)
fit5 <- clogit(update(base, ~. + Wcrf_C_Cal), data = distF.ph)
fit6 <- clogit(update(base, ~. + Wcrf_C_Cal), data = distM.ph)


modlistA <- list(fit7, fit8, fit9, fit10, fit11, fit12, fit13, fit14, fit15)
modlistA <- list(fit1, fit2, fit3, fit4, fit5, fit6)
modlistA <- list(fit10, fit11, fit12, fit13, fit14, fit15)

library(broom)
modsA <- map_df(modlistA, ~tidy(., exponentiate = T)) %>% filter(term == "Wcrf_C_Cal") %>%
  mutate_if(is.numeric, ~round(., 2)) %>% unite(OR.CI, estimate, conf.low, conf.high, sep = "-")

modscoresig <- cbind(modsA, mods)


# Remove unneeded variables from workspace
rm(list = ls()[!str_detect(ls(), ".ph|.both")])
rm(list = ls()[str_detect(ls(), "vars")])

# Load the predicted score tables for modelling
load("predicted_score_tables_sex.Rdata")

# Het. tests ----
# Formerly CRC_heterogeneity_tests.R

# Heterogeneity test for sex. Biocrates A and B and Fatty acids A
# Matching factors were age, sex, study centre, follow-up time since blood collection, fasting 
# status, menopausal status and phase of menstrual cycle at blood collection.

### Revised manuscript: studies A and B merged ###

library(tidyverse)
# Normal GLM with matching factors. Follow up time is meaningless for controls
# Menopause variables are incomplete
base <- Cncr_Caco_Clrt ~ Age_Blood + #Tfollowup + #Phase_Mnscycle + #Menopause + #Qe_Energy + 
  Center + Fasting_C + L_School + Smoke_Stat + Smoke_Int + Height_C 

# Biocrates small, large, fatty acids (large fasted subset did not improve OR)

# CRC endogenous: score and signature
library(lmtest)
fit1h <- glm(update(base, ~. + lab + comp1 + Sex), data = crc.ph, family = "binomial")
fit1i <- glm(update(base, ~. + lab + comp1 * Sex), data = crc.ph, family = "binomial")
lrtest(fit1h, fit1i)
# pHET = 0.029

library(broom)
s2 <- map_df(list(fit1h, fit1i), ~ tidy(., conf.int = T)) %>% filter(str_detect(term, "score."))
forest(s2$estimate, ci.lb = s2$conf.low, ci.ub = s2$conf.high, refline = 1, transf = exp) 


fit2h <- glm(update(base, ~. + Wcrf_C_Cal + Sex), data = crc.ph, family = "binomial")
fit2i <- glm(update(base, ~. + Wcrf_C_Cal * Sex), data = crc.ph, family = "binomial")
lrtest(fit2h, fit2i)
# pHET = 0.022

s2 <- map_df(list(fit2h, fit2i), ~ tidy(., conf.int = T)) %>% filter(str_detect(term, "Wcrf_"))
forest(s2$estimate, ci.lb = s2$conf.low, ci.ub = s2$conf.high, refline = 1, transf = exp) 


# CRC A fatty acids: score and signature
fit5h <- glm(update(base, ~. + comp2 + Sex), data = crc3.ph, family = "binomial")
fit5i <- glm(update(base, ~. + comp2 * Sex), data = crc3.ph, family = "binomial")
lrtest(fit5h, fit5i)
# pHET = 0.072

s2 <- map_df(list(fit5h, fit5i), ~ tidy(., conf.int = T)) %>% filter(str_detect(term, "score."))
forest(s2$estimate, ci.lb = s2$conf.low, ci.ub = s2$conf.high, refline = 1, transf = exp) 


fit6h <- glm(update(base, ~. + Wcrf_C_Cal + Sex), data = crc3.ph, family = "binomial")
fit6i <- glm(update(base, ~. + Wcrf_C_Cal * Sex), data = crc3.ph, family = "binomial")
lrtest(fit6h, fit6i)
# p = 0.36

s2 <- map_df(list(fit6h, fit6i), ~ tidy(., conf.int = T)) %>% filter(str_detect(term, "Wcrf_"))
forest(s2$estimate, ci.lb = s2$conf.low, ci.ub = s2$conf.high, refline = 1, transf = exp) 

# Colon endogenous: score and signature
fit7h <- glm(update(base, ~. + lab + comp1 + Sex), data = col.ph, family = "binomial")
fit7i <- glm(update(base, ~. + lab + comp1 * Sex), data = col.ph, family = "binomial")
lrtest(fit7h, fit7i)
# pHET = 0.03

fit8h <- glm(update(base, ~. + Wcrf_C_Cal + Sex), data = col.ph, family = "binomial")
fit8i <- glm(update(base, ~. + Wcrf_C_Cal * Sex), data = col.ph, family = "binomial")
lrtest(fit8h, fit8i)
# pHET = 0.002

# Prox colon endogenous: score and signature
fit7h <- glm(update(base, ~. + lab + comp1 + Sex), data = prox.ph, family = "binomial")
fit7i <- glm(update(base, ~. + lab + comp1 * Sex), data = prox.ph, family = "binomial")
lrtest(fit7h, fit7i)
# pHET = 0.21

fit8h <- glm(update(base, ~. + Wcrf_C_Cal + Sex), data = prox.ph, family = "binomial")
fit8i <- glm(update(base, ~. + Wcrf_C_Cal * Sex), data = prox.ph, family = "binomial")
lrtest(fit8h, fit8i)
# pHET = 0.12

# Dist endogenous: score and signature
fit7h <- glm(update(base, ~. + lab + comp1 + Sex), data = dist.ph, family = "binomial")
fit7i <- glm(update(base, ~. + lab + comp1 * Sex), data = dist.ph, family = "binomial")
lrtest(fit7h, fit7i)
# pHET = 0.12

fit8h <- glm(update(base, ~. + Wcrf_C_Cal + Sex), data = dist.ph, family = "binomial")
fit8i <- glm(update(base, ~. + Wcrf_C_Cal * Sex), data = dist.ph, family = "binomial")
lrtest(fit8h, fit8i)
# pHET = 0.005

# Rectal endogenous: score and signature
fit3h <- glm(update(base, ~. + lab + comp1 + Sex), data = rec.ph, family = "binomial")
fit3i <- glm(update(base, ~. + lab + comp1 * Sex), data = rec.ph, family = "binomial")
lrtest(fit3h, fit3i)
# pHET = 0.46

fit4h <- glm(update(base, ~. + Wcrf_C_Cal + Sex), data = rec.ph, family = "binomial")
fit4i <- glm(update(base, ~. + Wcrf_C_Cal * Sex), data = rec.ph, family = "binomial")
lrtest(fit4h, fit4i)
# pHET = 0.8346

# Colon fatty acids: score and signature
fit9h <- glm(update(base, ~. + Wcrf_C_Cal + Sex), data = col3.ph, family = "binomial")
fit9i <- glm(update(base, ~. + Wcrf_C_Cal * Sex), data = col3.ph, family = "binomial")
lrtest(fit9h, fit9i)
# pHET = 0.2836

fit10h <- glm(update(base, ~. + comp2 + Sex), data = col3.ph, family = "binomial")
fit10i <- glm(update(base, ~. + comp2 * Sex), data = col3.ph, family = "binomial")
lrtest(fit10h, fit10i)
# pHET = 0.1148

# Prox fatty acids: score and signature
fit11h <- glm(update(base, ~. + Wcrf_C_Cal + Sex), data = prox3.ph, family = "binomial")
fit11i <- glm(update(base, ~. + Wcrf_C_Cal * Sex), data = prox3.ph, family = "binomial")
lrtest(fit11h, fit11i)
# pHET = 0.4379

fit12h <- glm(update(base, ~. + comp2 + Sex), data = prox3.ph, family = "binomial")
fit12i <- glm(update(base, ~. + comp2 * Sex), data = prox3.ph, family = "binomial")
lrtest(fit12h, fit12i)
# pHET = 0.4328

# Dist fatty acids: score and signature
fit13h <- glm(update(base, ~. + Wcrf_C_Cal + Sex), data = dist3.ph, family = "binomial")
fit13i <- glm(update(base, ~. + Wcrf_C_Cal * Sex), data = dist3.ph, family = "binomial")
lrtest(fit13h, fit13i)
# pHET = 0.4864

fit14h <- glm(update(base, ~. + comp2 + Sex), data = dist3.ph, family = "binomial")
fit14i <- glm(update(base, ~. + comp2 * Sex), data = dist3.ph, family = "binomial")
lrtest(fit14h, fit14i)
# pHET = 0.1846

# Sensitivity anal.----
# Formerly CRC_sensitivity_analysis.R

library(survival)
library(broom)
load("predscore_df_subsite.Rdata")
base <- Cncr_Caco_Clrt ~ Qe_Energy + L_School + Smoke_Stat +  Height_C + strata(Match_Caseset)

### CRC1 and 2 merged for revision ###

# For smoke intensity, categories 8, 9 and 10 are collapsed into other (Smoke_Int)
# Replace smoking duration NAs with 0
crc.ph$Dur_Smok[is.na(crc.ph$Dur_Smok)] <- 0
col.ph$Dur_Smok[is.na(col.ph$Dur_Smok)] <- 0
rec.ph$Dur_Smok[is.na(rec.ph$Dur_Smok)] <- 0
crc3.ph$Dur_Smok[is.na(crc3.ph$Dur_Smok)] <- 0

# Dairy product intake = QGE05

# CRC A Fatty acids score, signature and by sex
fit1 <- clogit(update(base, ~. + comp2), data = crc3.ph)
fit2 <- clogit(update(base, ~. + Wcrf_C_Cal), data = crc3.ph)
fit3 <- clogit(update(base, ~. + comp2 + Smoke_Int), data = crc3.ph)
fit4 <- clogit(update(base, ~. + Wcrf_C_Cal + Smoke_Int), data = crc3.ph)
fit5 <- clogit(update(base, ~. + comp2 + Dur_Smok), data = crc3.ph)
fit6 <- clogit(update(base, ~. + Wcrf_C_Cal + Dur_Smok), data = crc3.ph)
fit7 <- clogit(update(base, ~. + comp2 + Qge05), data = crc3.ph)
fit8 <- clogit(update(base, ~. + Wcrf_C_Cal + Qge05), data = crc3.ph)
fit9 <- clogit(update(base, ~. + comp2), data = crc3t.ph)
fit10 <- clogit(update(base, ~. + Wcrf_C_Cal), data = crc3t.ph)

# Extra for revision: Signature model adjusted for score
fit11 <- clogit(update(base, ~. + Wcrf_C_Cal + comp2), data = crc3.ph)

# Fatty acids
fit12 <- clogit(update(base, ~. + Smoke_Int + comp2), data = crc3N.ph)
fit13 <- clogit(update(base, ~. + Smoke_Int + comp2), data = crc3O.ph)
fit14 <- clogit(update(base, ~. + Smoke_Int + comp2), data = crc3L.ph)
fit15 <- clogit(update(base, ~. + Smoke_Int + comp2), data = crc3H.ph)

sigmodsC <- map_df(list(fit1, fit3, fit5, fit7, fit11, fit12, fit13, fit14, fit15, fit9), ~tidy(., exponentiate = T)) %>% 
  filter(term == "comp2") %>% mutate_if(is.numeric, ~round(., 2)) %>% 
  unite(OR.CI, estimate, conf.low, conf.high, sep = "-") 

scomodsC <- map_df(list(fit2, fit4, fit6, fit8, fit10), ~tidy(., exponentiate = T)) %>% filter(term == "Wcrf_C_Cal") %>% 
  mutate_if(is.numeric, ~round(., 2)) %>% unite(OR.CI, estimate, conf.low, conf.high, sep = "-")

rm(list = ls(pattern = "fit"))
# CRC Biocrates
fit1 <- clogit(update(base, ~. + lab + comp1), data = crc.ph)
fit2 <- clogit(update(base, ~. + Wcrf_C_Cal), data = crc.ph)
fit3 <- clogit(update(base, ~. + lab + comp1 + Smoke_Int), data = crc.ph)
fit4 <- clogit(update(base, ~. + Wcrf_C_Cal + Smoke_Int), data = crc.ph)
fit5 <- clogit(update(base, ~. + lab + comp1 + Dur_Smok), data = crc.ph)
fit6 <- clogit(update(base, ~. + Wcrf_C_Cal + Dur_Smok), data = crc.ph)
fit7 <- clogit(update(base, ~. + lab + comp1 + Qge05), data = crc.ph)
fit8 <- clogit(update(base, ~. + Wcrf_C_Cal + Qge05), data = crc.ph)
fit9 <- clogit(update(base, ~. + lab + comp1), data = crcT.ph)
fit10 <- clogit(update(base, ~. + Wcrf_C_Cal), data = crcT.ph)

# Signature model adjusted for score
fit11 <- clogit(update(base, ~. + Wcrf_C_Cal + comp1), data = crc.ph)

# Models for high and low BMI and WCRF score
fit12 <- clogit(update(base, ~. + lab + comp1), data = crcN.ph)
fit13 <- clogit(update(base, ~. + lab + comp1), data = crcO.ph)
fit14 <- clogit(update(base, ~. + lab + comp1), data = crcL.ph)
fit15 <- clogit(update(base, ~. + lab + comp1), data = crcH.ph)

# Summarise in tables
sigmodsA <- map_df(list(fit1, fit3, fit5, fit7, fit11, fit12, fit13, fit14, fit15, fit9), ~tidy(., exponentiate = T)) %>% 
  filter(term == "comp1") %>% mutate_if(is.numeric, ~round(., 2)) %>% 
  unite(OR.CI, estimate, conf.low, conf.high, sep = "-")  

scomodsA <- map_df(list(fit2, fit4, fit6, fit8, fit10), ~tidy(., exponentiate = T)) %>% filter(term == "Wcrf_C_Cal") %>% 
  mutate_if(is.numeric, ~round(., 2)) %>% unite(OR.CI, estimate, conf.low, conf.high, sep = "-") 

rm(list = ls(pattern = "fit"))
### Colon Biocrates
fit1 <- clogit(update(base, ~. + lab + comp1), data = col.ph)
fit2 <- clogit(update(base, ~. + Wcrf_C_Cal), data = col.ph)
fit3 <- clogit(update(base, ~. + lab + comp1 + Smoke_Int), data = col.ph)
fit4 <- clogit(update(base, ~. + Wcrf_C_Cal + Smoke_Int), data = col.ph)
fit5 <- clogit(update(base, ~. + lab + comp1 + Dur_Smok), data = col.ph)
fit6 <- clogit(update(base, ~. + Wcrf_C_Cal + Dur_Smok), data = col.ph)
fit7 <- clogit(update(base, ~. + lab + comp1 + Qge05), data = col.ph)
fit8 <- clogit(update(base, ~. + Wcrf_C_Cal + Qge05), data = col.ph)
fit9 <- clogit(update(base, ~. + lab + comp1), data = colT.ph)
fit10 <- clogit(update(base, ~. + Wcrf_C_Cal), data = colT.ph)

# Signature model adjusted for score
fit11 <- clogit(update(base, ~. + lab + Wcrf_C_Cal + comp1), data = col.ph)

sigmodsA <- map_df(list(fit1, fit3, fit5, fit7, fit11, fit10), ~tidy(., exponentiate = T)) %>% filter(term == "comp1") %>% 
  mutate_if(is.numeric, ~round(., 2)) %>% unite(OR.CI, estimate, conf.low, conf.high, sep = "-")

scomodsA <- map_df(list(fit2, fit4, fit6, fit8, fit10), ~tidy(., exponentiate = T)) %>% filter(term == "Wcrf_C_Cal") %>% 
  mutate_if(is.numeric, ~round(., 2)) %>% unite(OR.CI, estimate, conf.low, conf.high, sep = "-")

rm(list = ls(pattern = "fit"))
# Rectal Biocrates
fit1 <- clogit(update(base, ~. + lab + comp1), data = rec.ph)
fit2 <- clogit(update(base, ~. + Wcrf_C_Cal), data = rec.ph)
fit3 <- clogit(update(base, ~. + lab + comp1 + Smoke_Int), data = rec.ph)
fit4 <- clogit(update(base, ~. + Wcrf_C_Cal + Smoke_Int), data = rec.ph)
fit5 <- clogit(update(base, ~. + lab + comp1 + Dur_Smok), data = rec.ph)
fit6 <- clogit(update(base, ~. + Wcrf_C_Cal + Dur_Smok), data = rec.ph)
fit7 <- clogit(update(base, ~. + lab + comp1 + Qge05), data = rec.ph)
fit8 <- clogit(update(base, ~. + Wcrf_C_Cal + Qge05), data = rec.ph)
fit9 <- clogit(update(base, ~. + lab + comp1), data = recT.ph)
fit10 <- clogit(update(base, ~. + Wcrf_C_Cal), data = recT.ph)
fit11 <- clogit(update(base, ~. + Wcrf_C_Cal + comp1), data = rec.ph)

sigmodsD <- map_df(list(fit1, fit3, fit5, fit7, fit11, fit9), ~tidy(., exponentiate = T)) %>% filter(term == "comp1") %>% 
  mutate_if(is.numeric, ~round(., 2)) %>% unite(OR.CI, estimate, conf.low, conf.high, sep = "-")
scomodsD <- map_df(list(fit2, fit4, fit6, fit8, fit10), ~tidy(., exponentiate = T)) %>% filter(term == "Wcrf_C_Cal") %>% 
  mutate_if(is.numeric, ~round(., 2)) %>% unite(OR.CI, estimate, conf.low, conf.high, sep = "-")


# Make final table to copy and paste into manuscript
sigmods <- rbind(sigmodsC, sigmodsA, sigmodsB) %>% mutate_if(is.numeric, ~round(., 2)) %>% 
  unite(OR.CI, estimate, conf.low, conf.high, sep = "-")

scoremods <- rbind(scomodsC, scomodsA, scomodsB) %>% mutate_if(is.numeric, ~round(., 2)) %>%
  unite(OR.CI, estimate, conf.low, conf.high, sep = "-")


# By score of individual components: colorectal, colon, rectal
library(broom)
library(tidyverse)
fit1 <- clogit(update(base, ~. + Wcrf_C_Cal), data = crc3.ph)


# Variable names: Wcrf_Bmi", "Wcrf_Pa", "Wcrf_Fwg_Cal", "Wcrf_Pf_Cal", "Wcrf_Meat_Cal", "Wcrf_Alc", "Wcrf_C_Cal"
fit1 <- clogit(update(base, ~. + Wcrf_Bmi), data = crc.ph)
fit2 <- clogit(update(base, ~. + Wcrf_Pa), data = crc.ph)
fit3 <- clogit(update(base, ~. + Wcrf_Fwg_Cal), data = crc.ph)
fit4 <- clogit(update(base, ~. + Wcrf_Pf_Cal), data = crc.ph)
fit5 <- clogit(update(base, ~. + Wcrf_Meat_Cal), data = crc.ph)
fit6 <- clogit(update(base, ~. + Wcrf_Alc), data = crc.ph)

ll <- list(fit1, fit2, fit3, fit4, fit5, fit6)
mods <- map_df(ll, ~tidy(., exponentiate = T)) %>% filter(grepl("Wcrf_", term)) %>%
  mutate_if(is.numeric, ~round(., 2))

fit1 <- clogit(update(base, ~. + Wcrf_Bmi), data = rec.ph)
fit2 <- clogit(update(base, ~. + Wcrf_Pa), data = rec.ph)
fit3 <- clogit(update(base, ~. + Wcrf_Fwg_Cal), data = rec.ph)
fit4 <- clogit(update(base, ~. + Wcrf_Pf_Cal), data = rec.ph)
fit5 <- clogit(update(base, ~. + Wcrf_Meat_Cal), data = rec.ph)
fit6 <- clogit(update(base, ~. + Wcrf_Alc), data = rec.ph)

ll <- list(fit1, fit2, fit3, fit4, fit5, fit6)
mods <- map_df(ll, ~tidy(., exponentiate = T)) %>% filter(grepl("Wcrf_", term)) %>%
  mutate_if(is.numeric, ~round(., 2))





# Metabolite-CRC associations for manuscript Table 2
# First need to run CRC_data_prep to get crc and crc3 objects
library(ggplot)
# Biocrates compounds. Split into quartiles with cut_number
df2 <- crc %>% #crc[, cols] %>%
  select(Glyceroph_Lysopc_A_C17_0, Glyceroph_Pc_Ae_C40_6, Glyceroph_Pc_Ae_C36_2, Glyceroph_Pc_Ae_C38_2, Aminoacid_Ser, 
         Glyceroph_Lysopc_A_C18_2, Aminoacid_Gly, Glyceroph_Pc_Ae_C40_3, 
         # low:
         Glyceroph_Pc_Aa_C32_1, Glyceroph_Pc_Aa_C38_4, Glyceroph_Pc_Aa_C36_4, Aminoacid_Glu, Glyceroph_Pc_Aa_C34_4, 
         Glyceroph_Pc_Aa_C40_4, Glyceroph_Pc_Ae_C38_3) %>%
  mutate_all(~cut_number(., n = 4, labels = 1:4))

# Define function to apply across quartiles (already matched by lab)
clrfun <- function(x)  { clogit(Cncr_Caco_Clrt ~ x + Bmi_C + Qe_Energy + L_School + Smoke_Stat + Smoke_Int + Height_C +
                                  Qe_Alc + Qge0701 + Qge0704 + strata(Match_Caseset), data = crc) }
fits <- apply(df2, 2, clrfun)
mods1 <- map_df(fits, ~tidy(., exponentiate = T)) %>% filter(grepl("x4", term)) %>% mutate_if(is.numeric, ~round(., 2)) %>%
  unite(OR.CI, estimate, conf.low, conf.high, sep = "-")

# Fatty acids. Modelled continuously because few values.
df5 <- crc3[, colnames(ctrlC)] %>% 
  select(P17_0, P15_0, P15_1, P22_5n_6, P18_1n_9c, P16_1n_7c_n_9c, P16_0, P20_3n_9, P22_1n_9) %>% scale()

clrfun1 <- function(x)  { 
  clogit(Cncr_Caco_Clrt ~ x + Bmi_C + Qe_Energy + L_School + Smoke_Stat + Smoke_Int + Height_C + Qe_Alc +
           Qge0701 + Qge0704 + strata(Match_Caseset), data = crc3) 
}
fits <- apply(df5, 2, clrfun1)
mods <- map_df(fits, ~tidy(., exponentiate = T)) %>% filter(grepl("x", term)) %>% mutate_if(is.numeric, ~round(., 2)) %>%
  unite(OR.CI, estimate, conf.low, conf.high, sep = "-")

# MS figures----
# Formerly CRC_manuscript_figs.R

# Figures for manuscript table 2
library(ggplot2)
library(tidyverse)

# Load coefficients as generated by CRC_get_signatures.R
load("coefficient_tables.Rdata")
sig1 <- pltdata %>% filter(abs(Coefficient) > 0.019)
sig2 <- faplot2 %>% filter(abs(Coefficient) > 0.02)
sigall <- bind_rows("Fatty acid signature" = sig2, "Endogenous metabolite signature" = sig1, 
                    .id = "Signature")

# Signature plot with arrows (facetted)
ggplot(sigall, aes(y = 0, yend = Coefficient, x = reorder(compound, -abs(Coefficient)), 
                   xend = reorder(compound, Coefficient), colour = Signature)) + 
  geom_hline(yintercept = 0, linetype = "dotted") +
  geom_segment(arrow = arrow(length = unit(0.06, "inches"), type = "closed"), arrow.fill = NULL) + 
  scale_colour_manual(values = c("blue", "red")) +
  geom_vline(xintercept = 0, linetype = "dashed") +
  ylab("Coefficient on first PLSR latent variable") + theme_classic() +
  theme(axis.title.x = element_blank()) +
  #facet_grid(. ~ fct_inorder(Signature), scales = "free") +
  facet_wrap(. ~ fct_inorder(Signature), scales = "free") +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5, colour = "black"),
        axis.ticks.y = element_line(), strip.background = element_blank(),
        #axis.ticks.x = element_blank(), 
        strip.text.x = element_text(size = 11),
        panel.spacing = unit(2, "lines"),
        legend.position = "none") #+ ggtitle("B")


nlist <- c("Maintain normal\nbody weight", "Be moderately\nphysically\nactive", 
           "Limit foods\nthat promote\nweight gain", "Eat mostly\nplant foods", "Limit red and\nprocessed meat", "Avoid\nalcohol", 
           "Overall WCRF/\nAICR score")

# Correlation plot with CIs (method 1) (probably the best)
#source("WCRF_components_sig.R")
pcor.ci <- readRDS("df_wcrf_corr2.rds")


ggplot(pcor.ci, aes(x=fct_inorder(component) %>% fct_relabel(~nlist), 
                    y=estimate, colour = fct_inorder(Model), shape = fct_inorder(Model))) + 
  geom_errorbar(width=0.2, aes(ymin = conf.low, ymax = conf.high), 
                position= position_dodge(width = 0.7), colour = "grey40") +
  geom_point(size = 3, position = position_dodge(width = 0.7), stroke = 1) +
  scale_colour_manual(values = c("red","blue","blue")) +
  scale_shape_manual(values = c(18,18,20)) +
  geom_hline(yintercept = 0, linetype = "dashed", colour = "grey") +
  xlab("WCRF/AICR recommendation") +
  ylab("Partial Pearson correlation") + theme_classic() +
  scale_x_discrete(position = "bottom") +
  facet_wrap(. ~ fct_inorder(component), scales = "free_x", nrow = 1) +
  theme(legend.position = "right",
        axis.text.x = element_text(colour = "black"),
        axis.ticks.x = element_blank(),
        legend.title = element_blank(),
        legend.key.height=unit(1.2, "cm"),
        legend.box.background = element_rect(colour = "black", size = 1),
        strip.text.x = element_blank()) #+ ggtitle("C")


# Correlation plot with CIs (method 2)
ggplot(pcor.ci, aes(x = Model, y=estimate, fill = fct_inorder(Model), 
                    shape = fct_inorder(Model))) + 
  geom_point() +
  geom_errorbar(width=0.2, aes(ymin = conf.low, ymax = conf.high), 
                position= position_dodge(width = 0.5), colour = "grey60") +
  geom_point(size = 2, position = position_dodge(width = 0.5)) +
  scale_fill_manual(values = c("red","blue","blue")) +
  scale_shape_manual(values = c(21,21,24)) + theme_classic() +
  geom_hline(yintercept = 0, linetype = "dashed") +
  ylab("Partial Pearson correlation") +
  scale_x_discrete(position = "bottom") +
  facet_grid(. ~ fct_inorder(component) %>% fct_relabel(~nlist), scales = "free") +
  theme(axis.title.x = element_blank(), legend.position = "right",
        axis.text.x = element_blank(), axis.ticks.x = element_blank(),
        strip.background.x = element_blank(),
        legend.title = element_blank(),
        legend.box.background = element_rect(colour = "black"))

# Old (no facetting)
ggplot(sig1, aes(x = 0, xend = Coefficient, y = reorder(compound, Coefficient), 
                 yend = (reorder(compound, Coefficient)))) + 
  geom_segment(arrow = arrow(length = unit(0.1, "inches"), type = "closed"),
               arrow.fill = "blue") + 
  theme_bw() +
  geom_vline(xintercept = 0, linetype = "dashed") +
  xlab("Coefficient on first PLSR latent variable") +
  theme(axis.title.y = element_blank())

ggplot(sig2, aes(x = 0, xend = Coefficient, y = reorder(compound, Coefficient), 
                 yend = (reorder(compound, Coefficient)))) + 
  geom_segment(arrow = arrow(length = unit(0.1, "inches"), type = "closed"),
               arrow.fill = "red") + 
  theme_bw() +
  geom_vline(xintercept = 0, linetype = "dashed") +
  xlab("Coefficient on first PLSR latent variable") +
  theme(axis.title.y = element_blank())

# Coefficient plots

library(ggrepel)

# To change legend order, reorder factor levels
corr <- cor(Bioc0[, -1], Bioc0[, 1], method = "spearman")
pltdata1 <- pltdata %>% mutate(compound1 = ifelse(abs(Coefficient) > 0.022, compound, NA), corr)

p1 <- 
  ggplot(pltdata1, aes(x = Coefficient, y = corr, shape = str_wrap(class, 20))) + 
  geom_point() + theme_bw() +
  theme(panel.grid.major = element_blank(), legend.title = element_blank(),
        legend.position = c(0.85, 0.25),
        legend.key.size = unit(0, "lines")) +
  scale_shape_manual(values = c(15, 1, 19, 25, 22, 14, 3, 4)) +
  xlab("Coefficient on first PLSR latent variable (p1)") + 
  ylab("Spearman correlation WCRF score - metabolite") +
  geom_hline(yintercept = 0, linetype = "dashed") +
  geom_vline(xintercept = 0, linetype = "dashed") +
  geom_text_repel(aes(label = compound1), size = 3) #+
#ggtitle("A")

corr1 <- cor(Facid[, -1], Facid[, 1], method = "spearman")
faplot1 <- faplot %>% mutate(compound1 = ifelse(abs(Coefficient) > 0.045, compound, NA), corr1)

p2 <- ggplot(faplot1, aes(x = Coefficient, y = corr1, shape = class)) + geom_point() + 
  theme_bw() + xlim(-0.19, 0.19) + ylim(-0.22, 0.22) +
  theme(panel.grid.major = element_blank(), legend.title = element_blank(),
        legend.position = c(0.85, 0.2),
        legend.key.size = unit(0, "lines")) +
  scale_shape_manual(values = c(6, 16, 1, 4, 3)) +
  xlab("Coefficient on 1st PLSR latent variable (p1)") + 
  ylab("Spearman correlation WCRF score - metabolite") +
  geom_hline(yintercept = 0, linetype = "dashed") +
  geom_vline(xintercept = 0, linetype = "dashed") +
  geom_text_repel(aes(label = compound1), size = 3) #+
#ggtitle("B")


# Generate aligned plots
library(cowplot)
plot_grid(p1, p2, nrow = 2, labels = c("A", "B"), label_size = 12) %>%
  save_plot("s_plots.pdf")

both <- align_plots(p1, p2, align = "hv", axis = "tblr")
ggdraw(both[[1]])
ggsave("endogenous.png", height = 100, width = 180, units = "mm")
ggdraw(both[[2]])
ggsave("FAs.png", height = 100, width = 180, units = "mm")

# Venn diagram for compounds
library(VennDiagram)
venn.diagram(list(Controls = colnames(ctrlA), CRC1 = colnames(ctrlB), CRC2 = colnames(ctrls)), 
             imagetype = "png", filename = "compound_venn.png", 
             height=150, width=150, units="mm", cat.fontfamily = "", fontfamily = "")

# Biocrates and FAs together in facetted plot (couldn't put classes in right order)
all <- bind_rows(Endogenous = pltdata, `Fatty acids` = faplot, .id = "Metabolites")

ggplot(all, aes(y = Class, x = Coefficient, colour = Class)) + 
  geom_jitter(height = 0.05) + 
  xlab("Coefficient from first PLS latent variable") +
  theme_bw() + 
  geom_vline(xintercept = 0, linetype = "dashed") +
  theme(legend.position = "none", axis.title.y = element_blank(),
        axis.line.y = element_line()) +
  facet_grid(Metabolites ~ ., scales = "free", space = "free", switch = "y")

# Old coefficient scatter (top and bottom percentiles)

# Vector of black and grey for plot points
vec <- c( rep("black", n.low), rep("grey", nrow(dat) - nrow(dat1)), rep("black", n.high) )

# Now plot data, adding text
plot(sort(coeff$value), pch = 17, col=vec, xlab = "", ylab = "Coefficient",
     main = paste(nrow(mod$scores), "fasted subjects, optimal dimensions =", lv))
# High and low labels
text(nrow(dat) : (nrow(dat) - n_top), df1$Coefficient, df1$compound, pos=2, cex = 0.6)
text(1:nrow(df2), df2$Coefficient, df2$compound, pos=4, cex=0.6)
abline(a=0, b=0, lty = "dotted")

# Signature partial corr.----
# Formerly WCRF_components_sig.R

library(tidyverse)
#load("predicted_score_tables_sex.Rdata")
source("CRC_prep_data_rev.R")
load("predscore_df_subsite.Rdata")

# Maintain body weight within the normal range
# Be moderately physically active
# Avoid sugary drinks: Wcrf_Drinks_Cal
# Limit consumption of energy-dense foods and avoid sugary drinks
# Eat at least 5 portions of non-starchy vegetables/fruits every day: Wcrf_Fv_Cal
# Eat unprocessed cereals (grains) and/or pulses (legumes): Wcrf_Fibt_Cal
# Eat mostly foods of plant origin
# Energy dense foods: Wcrf_Ed
# Animal foods: Wcrf_Meat_Cal
# Avoid alcohol
# Breastfeed infants exclusively up to 6 months: Wcrf_Bf
# Overall score

varlist <- c("Wcrf_Bmi", "Wcrf_Pa", "Wcrf_Fwg_Cal", "Wcrf_Pf_Cal", "Wcrf_Meat_Cal", "Wcrf_Alc", "Wcrf_C_Cal")
nlist <- c("Maintain normal\nbody weight", "Be moderately\nphysically active", 
           "Limit foods that\npromote weight gain", "Eat mostly\nplant foods", "Limit red and\nprocessed meat", "Avoid\nalcohol", 
           "Overall WCRF/AICR\nscore")

crc1f <- crc3.ph %>% ungroup %>% filter(Cncr_Caco_Clrt == 0) %>% select(one_of(varlist)) %>% as.matrix
crc1b <- crc.ph %>% ungroup %>% filter(Cncr_Caco_Clrt == 0) %>% select(one_of(varlist)) %>% as.matrix

comp2.ct1 <- crc3.ph$comp2[crc3.ph$Cncr_Caco_Clrt == 0]
comp1.ct1 <- crc.ph$comp1[crc.ph$Cncr_Caco_Clrt == 0]

# Partial correlations. Subset controls only from prediction tables
dat1 <- crc3.ph[crc3.ph$Cncr_Caco_Clrt == 0, ]
dat2 <- crc.ph[crc.ph$Cncr_Caco_Clrt == 0, ]

# Function to get partial correlations, omitting NAs
get.pcor <- function(x, dat, predscore) {
  mod1 <- lm(x[!is.na(x)] ~         L_School + Qe_Energy + Smoke_Stat + Smoke_Int + Height_C, data = dat[!is.na(x), ] )
  mod2 <- lm(predscore[!is.na(x)] ~ L_School + Qe_Energy + Smoke_Stat + Smoke_Int + Height_C, data = dat[!is.na(x), ] )
  output <- cor.test(residuals(mod1), residuals(mod2), method = "pearson")
}

# Biocrates data also adjusted for laboratory of analysis
get.pcor.lab <- function(x, dat, predscore) {
  mod1 <- lm(x[!is.na(x)] ~         L_School + lab + Qe_Energy + Smoke_Stat + Smoke_Int + Height_C, data = dat[!is.na(x), ] )
  mod2 <- lm(predscore[!is.na(x)] ~ L_School + lab + Qe_Energy + Smoke_Stat + Smoke_Int + Height_C, data = dat[!is.na(x), ] )
  output <- cor.test(residuals(mod1), residuals(mod2), method = "pearson")
}

# Run function and extract data 
pcor1 <- apply(crc1f, 2, function(x, dat, predscore) get.pcor(x, dat1, comp2.ct1)) 
pcor2 <- apply(crc1b, 2, function(x, dat, predscore) get.pcor.lab(x, dat2, comp1.ct1))

library(broom)
pcor.ci <- bind_rows("Fatty acid\nsignature" = map_df(pcor1, tidy), "Endogenous\nmetabolite\nsignature" = map_df(pcor2, tidy), 
                     .id = "Model") %>% bind_cols(tibble(component = rep(varlist, 2)))

# Plot data: partial correlation
saveRDS(pcor.ci, "wcrf_corr_df.rds")
#load("df_wcrf_correlations.rds")


# Raw correlations
library(psych)
fa1 <- corr.test(comp2.ct1, crc1f, use = "pairwise.complete.obs")
em1 <- corr.test(comp1.ct1, crc1b, use = "pairwise.complete.obs")
em2 <- corr.test(comp1.ct2, crc2b, use = "pairwise.complete.obs")

cor.ci <- bind_rows(fa1$ci, em1$ci, em2$ci, .id = "Model") %>% 
  bind_cols(tibble(component = rep(varlist, 3)))

# Plot data: raw correlation
library(ggplot2)

cor.ci %>%
  ggplot(aes(x=fct_inorder(component), y=r, shape = signature)) + 
  theme_bw() +
  geom_hline(yintercept = 0, linetype = "dashed") +
  geom_point(size = 2, position = position_dodge(width = 0.5)) +
  ylab("Pearson correlation") + xlab("Fatty acid or endogenous signature") +
  geom_errorbar(width=0.2, aes(ymin = lower, ymax = upper), 
                position= position_dodge(width = 0.5), colour="black") +
  theme(axis.title.x = element_blank())


# Old: facetted
ggplot(cor.ci, aes(x = Model, y = r)) + #geom_point() + 
  geom_hline(yintercept = 0, linetype = "dashed") +
  geom_errorbar(width=.2, aes(ymin=lower, ymax=upper), colour="black") +
  geom_point() + theme_minimal() +
  ylab("Pearson correlation") + xlab("Fatty acid or endogenous signature") +
  facet_grid(. ~ fct_inorder(component), scales = "free_x") +
  theme(panel.spacing=unit(1.5,"lines"))


# Supp table CVs----
# Formerly Table_all_cmpds.R

# Supplemental table of compounds and CVs
# Biocrates from IARC and Helmholtz, FAs from IARC
library(readxl)
library(janitor)
library(tidyverse)

bioc <- read_csv("Biocrates_cmpd_metadata.csv") 
facids <- read_csv("FA_compound_data.csv")

# Helmholtz Biocrates CVs. Read files, put together and calculate intrabatch mean
cvs <- list.files("helmoltz_CVs") %>% 
  map_df( ~ read_delim(paste("helmoltz_CVs", ., sep = "/"), delim = ";")) %>% 
  select(C3:H1) %>% summarise_all(funs(mean))

# Get list of compound names from full data
full <- read_delim("helmholtz_CV_data/2018-01-16_Conc_P160262_P-17044-BBMRI-P1_Kit-LOD.csv", 
                   skip = 1, delim = ";") %>% select(-contains("Status")) %>% select(C3:H1)

# Check names are the same and replace
df <- data.frame(n1 = colnames(cvs), n2 = colnames(full))
colnames(cvs) <- colnames(full)
helm <- gather(cvs, key = "shortname") %>% mutate(CV_Helm = round(value*100, 1))

# IARC Biocrates CVs. 2nd sheet for serum
cvbioc <- read_xlsx("Epic CRC_QC Results.xlsx", skip = 1, sheet = 2) %>% 
  filter(str_detect(Compounds, "CV% Interbatch")) %>% slice(2) %>% select(-(1:4))
iarc <- gather(cvbioc, key = "shortname") %>% slice(-(1:4)) %>% mutate(CV_Iarc = round(value, 1))

# Get table of coefficients for signatures
load("coefficient_tables.Rdata")

# Note: the following are excluded from signature because of absence from CRC
# 159 compounds reduced to 155
#"Biogenic_Nitro_Tyr" "Biogenic_Ac_Orn"    "Biogenic_Met_So"    "Biogenic_Total_Dma"

# Join two sets of CVs to Biocrates metadata
bioc1 <- bioc %>% left_join(iarc, by = "shortname")
bioc2 <- bioc1 %>% left_join(helm, by = "shortname")

pltdata1 <- pltdata %>% select(compound:Coefficient)

# Join coefficients table to biocrates metadata
bioc3 <- pltdata1 %>% left_join(bioc2, by = c("compound" = "displayname")) %>%
  select(class, compound, Coefficient, CV_Iarc, CV_Helm) %>% arrange(class)


# Fatty acids
cvfa <- read_xlsx("FA-QC-CV.xlsx", col_types = c("text", "numeric")) %>% 
  rename(displayname = "Fatty Acid") %>% mutate(CV = round(`CV%`, 1))

facidCV1 <- faplot %>% left_join(cvfa, by = c("compound" = "displayname")) %>%
  select(class, compound, Coefficient, CV_Iarc = CV) %>% arrange(class)

# Bind two tables together
allCVdat <- bind_rows(Endogenous.metabolites = bioc3, Fatty.acids = facidCV1, .id = "Platform")


# Get CVs for first two ACs (not supplied with data)
fullcvs <- list.files("helmholtz_CV_data") %>% 
  map_df( ~ read_delim(paste("helmholtz_CV_data", ., sep = "/"), delim = ";", skip = 1)) %>% 
  filter(`Sample Identification` == "Ref_Plasma-Hum_PK3") %>%
  select(KitBarcodeNr, C0, C2, C3) %>% mutate_at(vars(C0, C2, C3), funs(as.numeric))

# Overall CV (not used in table)
overallcvs <- fullcvs %>% summarise_all(~ sd(.)/mean(.))

#Batch average CVs
batchcvs <- fullcvs %>% group_by(KitBarcodeNr) %>% summarise_all(~ sd(.)/mean(.))
batchcvs %>% summarise_all(mean)
#C0: 7.0%, C2: 7.1%

# Calc score models----
# Formerly CRC_models_calc_score.R
# Models for questionnaire calculated score
# Models for colorectal cancer status and WCRF score. First source CRC_data_prep.R which preps the two CCs.
source("CRC_data_prep.R")
source("Metabolic_signature_WCRF.R")

library(haven)
library(tidyverse)
meta <- read_dta("clrt_caco.dta")

# Convert categorical variables to factors
var.list  <- c("Country", "Center", "Sex", "Match_Caseset", "Smoke_Stat", "L_School")
meta      <- meta %>% mutate_at(vars(var.list), as.factor)

# Get names of individual contributor variables to WCRF score
wcrf_vars <- c("Wcrf_Bmi", 
               "Wcrf_Pa", 
               "Wcrf_Fwg_Cal", 
               "Wcrf_Pf_Cal", 
               "Wcrf_Fv_Cal", 
               "Wcrf_Fibt_Cal",
               "Wcrf_Meat_Cal", 
               "Wcrf_C_Cal")

# Get variable descriptions (for forest plot)
scorecomp <- c("1. Body fatness", 
               "2. Physical activity", 
               "3. Energy density/sugary drinks", 
               "4. FV intake", 
               "5. Foods of plant origin", 
               "6. Fibre intake", 
               "7. Meat intake", 
               "       Overall WCRF score (cal.)")
scorecomp2 <- c(scorecomp, "       Signature metabolites")
xtitle     <- "Odds ratio (per unit increase in score)"
xtitle2    <- "Hazard ratio (per unit increase in score)"

library(survival)
library(broom)
# Function to run model case-control status from main WCRF score components and get tidy output
# Warning: takes a while to run
clr_crc <- function(dat, sig = "none") {
  
  # model with covariates
  mat <- dat %>% select(wcrf_vars)
  
  # bind predicted scores from Metabolic_signature_WCRF.R, print dimensions as check
  mat <- if(sig == "large") cbind(mat, large$score.2.comps) else if(sig == "small") cbind(mat, small$score.2.comps) else mat
  print(paste("dimensions:", dim(mat)))
  
  # Run model and apply across matrix
  fit <- function(x) clogit(Cncr_Caco_Clrt ~ x + Qe_Energy + L_School + Smoke_Stat + strata(Match_Caseset), data = dat)
  multifit <- apply(mat, 2, fit)
  
  # subset risk estimate for score only
  output <- map_df(multifit, tidy) %>% filter(term == "x")
}
t1 <- clr_crc(meta, sig = "none") %>% mutate(Study = scorecomp)
t2 <- clr_crc(crc1, sig = "small") %>% mutate(Study = scorecomp2)
t3 <- clr_crc(crc2, sig = "large") %>% mutate(Study = scorecomp2)

# Plot output of each study as forest
library(metafor)
par(mar=c(5,4,1,2))
# Full study ~ 14,000 subjects
forest(t1$estimate, ci.lb = t1$conf.low, ci.ub = t1$conf.high, refline = 1, xlab = xtitle, 
       transf = exp, pch = 18, psize = 1.5, slab = t1$Study) 
# Large study ~ 2370 subjects
forest(t2$estimate, ci.lb = t2$conf.low, ci.ub = t2$conf.high, refline = 1, xlab = xtitle, 
       transf = exp, pch = 18, psize = 1.5, slab = t2$Study) 
# SMall study  ~ 900 subjects
forest(t3$estimate, ci.lb = t3$conf.low, ci.ub = t3$conf.high, refline = 1, xlab = xtitle, 
       transf = exp, pch = 18, psize = 1.5, slab = t3$Study) 

hh <- par("usr")
text(hh[1], nrow(t2) + 2, "Component of score", pos = 4)
text(hh[2], nrow(t2) + 2, "OR [95% CI]", pos = 2)



# Plot only overall score and signature from each study
t4 <- bind_rows(t1[nrow(t1), ], t2[(nrow(t2)-1):nrow(t2), ], t3[(nrow(t3)-1):nrow(t3), ])

par(mar=c(5,4,1,2))
forest(t4$estimate, ci.lb = t4$conf.low, ci.ub = t4$conf.high, refline = 1, 
       xlab = "Odds ratio (per unit increase in score)", pch = 18, transf = exp, psize = 1.5,
       slab = c("WCRF score, full study", "WCRF score, small", "Metabolic signature, small", 
                "WCRF score, large", "Metabolic signature, large"))

hh <- par("usr")
text(hh[1], nrow(t4) + 2, "Parameter", pos = 4)
text(hh[2], nrow(t4) + 2, "OR [95% CI]", pos = 2)

# Fixed-effects meta-analysis of signatures
ma1 <- rma(estimate, sei = std.error, data=t4, method="FE", subset = 2:3)
# Random-effects meta-analysis of signatures
ma2 <- rma(estimate, sei = std.error, data=t4, method="REML", subset = 2:3)

# Of scores
ma3 <- rma(estimate, sei = std.error, data=t4, method="FE", subset = c(2,4))
ma4 <- rma(estimate, sei = std.error, data=t4, method="REML", subset = c(2,4))

par(mfrow = c(1,2))
slab <- c("WCRF small \nmetabolomics", "WCRF large \nmetabolomics")
forest(ma1, transf = exp, refline = 1, slab = slab, xlab = "OR")
forest(ma2, transf = exp, refline = 1, slab = slab, xlab = "OR")
forest(ma3, transf = exp, refline = 1, slab = slab, xlab = "OR")
forest(ma4, transf = exp, refline = 1, slab = slab, xlab = "OR")




#---- Prospective associations WCRF score and CRC

# Cox model for whole dataset
# Create the survival object. In stata, age at exit, age at birth, age at recruitment, and CRC yes/no are used
# survobj <- Surv(time = df$Age_Recr, time2 = ..., event = ...)
# join the WCRF scores
library(haven)
library(tidyverse)
wcrf <- read_dta("wcrf_score.dta")
fullepic <- read_dta("D:/full epic.dta")
fullepic <- fullepic %>% left_join(wcrf, by = "Idepic")

# Cox proportional hazards models
# model with covariates
mat <- fullepic %>% select(wcrf_vars)

library(survival)
# Create the survival object. In stata, age at exit, age at birth, age at recruitment, and CRC yes/no are used
survobj <- Surv(time = fullepic$Age_Recr, time2 = fullepic$Agexit, event = fullepic$Cncr_Evt_Clrt)

# Run model and apply across WCRF matrix
fit <- function(x) coxph(survobj ~ x + Smoke_Intensity + L_School + QE_ENERGY + strata(Sex, Center.x), data = fullepic)
multifit <- apply(mat, 2, fit)

# subset risk estimate for score only
library(broom)
t1 <- map_df(multifit, tidy) %>% filter(term == "x")

library(metafor)
par(mar=c(5,4,1,2))
forest(t1$estimate, ci.lb = t1$conf.low, ci.ub = t1$conf.high, refline = 1, xlab = xtitle2, 
       transf = exp, pch = 18, psize = 1.5, slab = scorecomp) 
hh <- par("usr")
text(hh[1], nrow(t1) + 2, "Component of score", pos = 4)
text(hh[2], nrow(t1) + 2, "HR [95% CI]", pos = 2)


# Variability----
# Formerly CRC_variability.R

# Find main sources of variability in EPIC controls with PCPR2
# Biocrates: First dataset, 3771 obs; updated November 2018 7191 obs
source("CRC_prep_data_rev.R")

# 1741 fasted controls without Greece
# Subset compounds (X data)
cmpds <- ctrl %>% 
  select(matches("Acylcarn_|Aminoacid_|Biogenic_|Glyceroph_|Sphingo_|Sugars_"), -starts_with("Outdq"))

zerocols <- apply(cmpds, 2, function(x) sum(x, na.rm = T)) != 0
concs <- cmpds[, zerocols]

# Subset metadata (Z data)
meta <- ctrl
meta$Center <- as.factor(meta$Center)
meta$Sex <- as.factor(meta$Sex)
meta$Study <- droplevels(meta$Study)
meta$BMI <- meta$Bmi_C
meta$Batch <- meta$batch_no
meta <- meta %>% select(Center, Batch, Sex, BMI, Study)

# Prepare controls matrix. Replace zero, impute with half mins, scale
concs[concs == 0] <- NA
library(zoo)
concs1 <- na.aggregate(concs, FUN = function(x) min(x)/2)
logconcs <- log2(concs1) %>% scale

# Run pcpr2
library(pcpr2)
obj <- runPCPR2(logconcs, meta)

# Adjust using residuals method and rerun
library(lme4)
adj   <- function(x) residuals(lmer(x ~ Center + Batch + Sex + (1|Study), data = meta))
adjmat <- apply(logconcs, 2, adj)
obj2 <- runPCPR2(adjmat, meta)

# Fatty acids
library(tidyverse)
# Exclude compounds with many missings. New version from Carine received 18/11/2018 with technical covariates
fa.ctrl <- readRDS("FA_WCRF_scores1.rds") %>% filter(STUDY != "Colorectum") %>% filter(Country != 6)
#fa.ctrl$N_Serie <- as.numeric(fa.ctrl$N_Serie)

epic.vars <- read.csv("full_epic_main_vars.csv")
fa.ctrl <- left_join(fa.ctrl, epic.vars, by = "Idepic", suffix = c("", "_1"))

concs <- fa.ctrl %>% select(P14_0 : PCLA_9t_11c, -P24_0, -P20_0) %>% as.matrix
concs[concs == 0] <- NA
library(zoo)
concs1 <- na.aggregate(concs, FUN = function(x) min(x)/2)

# categorical variables to factors
meta <- fa.ctrl #%>% select(Center, LABO, Sex, Bmi_C, STUDY)
meta$Laboratory <- meta$LABO
meta$BMI <- meta$Bmi_C
meta$Study <- meta$STUDY
meta$Laboratory <- meta$LABO
meta <- meta %>% select(Center, Laboratory, Sex, BMI, Study)

var.list <- c("Center", "Study", "Laboratory")
meta <- meta %>% mutate_at(vars(var.list), as.factor)

# Run pcpr2
output1 <- runPCPR2(concs1, meta)

# Adjust matrix for study, centre, batch, sex for Biocrates, subset calibrated scores 
library(lme4)
adj  <- function(x) residuals(lmer(x ~ Laboratory + Study + (1|Center), data = meta))
adjmat <- apply(concs1, 2, adj)
output2 <- runPCPR2(adjmat, meta)


# Plot figure for supp material (text positions determined by trial and error)
par(mfrow = c(2,2), oma = c(0, 0, 2, 0))
plot(obj, main = "Raw metabolite matrix")
plot(obj2, main = "After transformation")
mtext("A)", outer = TRUE, adj = 0.02)
mtext("B)", outer = FALSE, adj = -2, line = -20)

plot(output1, main = "Raw metabolite matrix")
plot(output2, main = "After transformation")

# Forest plots----
# Formerly CRC_more_forests.R

# Signature only for Biocrates small and large (all subjects in study) and fatty acids small
library(broom)
library(tidyverse)
library(metafor)
rows <- c("Score all", "    Score M", "    Score F", "Signature all", "    Sig M", "    Sig F", 
          "    Sex sp. sig M", "    Sex sp. sig F")
lab = "OR per unit increase"

# Biocrates A
# Score, score M, score F, Sig all, sig M, sex specific sig M, sig F, sex specific sig F
ll2 <- list(fit2, fit2m, fit2f, fit1, fit1m, fit1f, fit1ms, fit1fs)
t2 <- map_df(ll2, tidy) %>% filter(str_detect(term, "score.|Wcrf_"))

#par(mar=c(5,4,1,2))
forest(t2$estimate, ci.lb = t2$conf.low, ci.ub = t2$conf.high, refline = 1, psize = 1.5,
       xlab = lab, pch = 18, slab = rows, transf = exp, main = "CRC A Biocrates")

# Biocrates B
ll2 <- list(fit4, fit4m, fit4f, fit3, fit3m, fit3f, fit3ms, fit3fs)
t2 <- map_df(ll2, tidy) %>% filter(str_detect(term, "score.|Wcrf_"))

forest(t2$estimate, ci.lb = t2$conf.low, ci.ub = t2$conf.high, refline = 1, psize = 1.5,
       xlab = lab, pch = 18, slab = rows, transf = exp, main = "CRC B Biocrates")

# Fatty acids A
ll2 <- list(fit5, fit5m, fit5f, fit6, fit6m, fit6f, fit5fs)
t2 <- map_df(ll2, tidy) %>% filter(str_detect(term, "score.|Wcrf_"))

forest(t2$estimate, ci.lb = t2$conf.low, ci.ub = t2$conf.high, refline = 1, psize = 1.5,
       xlab = lab, pch = 18, slab = rows[-7], transf = exp, main = "CRC A fatty acids")

# Baseline char CC----
# Formerly CRC_baseline.R. Baseline characteristics for the 2 CRC metabolomics studies
# 8 unpaired samples have been removed for CRC1

source("CRC_prep_data.R")
library(tidyverse)

#Get predicted WCRF score
load("predscore_df_subsite1.Rdata")
crc.ph0 <- crc.ph %>% select(Idepic, comp1)
crc.ph1 <- crc3.ph %>% select(Idepic, comp2)

# Generate time to diagnosis and tumour site variables for both studies
crc1a <- crc1 %>% group_by(Match_Caseset) %>% filter(n() == 2) %>% ungroup() %>%
  select(Idepic, Match_Caseset, Sex, location, Age_Blood, Tfollowup, Height_C, Bmi_C, Waist_C, Qe_Energy, Bdg1clrt,
         Country, Pa_Mets, Smoke_Stat, Qe_Alc, Wcrf_C_Cal, Cncr_Caco_Clrt) %>% 
  mutate(Study = "CRC1", Bdg1clrt_hist = ifelse(Bdg1clrt %in% c(53,54,55,56,60,70), 1, 0)) %>%
  left_join(crc.ph0, by = "Idepic") %>%
  left_join(crc.ph1, by = "Idepic")

crc2a <- crc2 %>%
  select(Idepic, Match_Caseset, Sex, location, Age_Blood, Tfollowup, Height_C, Bmi_C, Waist_C, Qe_Energy, Bdg1clrt,
         Country, Pa_Mets, Smoke_Stat, Qe_Alc, Wcrf_C_Cal, Cncr_Caco_Clrt) %>% 
  mutate(Study = "CRC2", Bdg1clrt_hist = ifelse(Bdg1clrt %in% c(53,54,55,56,60,70), 1, 0)) %>%
  left_join(crc.ph0, by = "Idepic")

# Merge for revised submission and add predicted WCRF score
crc.sig <- bind_rows(crc1a, crc2a)
crc.sig$Tfollowup <- as.numeric(crc.sig$Tfollowup)

library(qwraps2)
options(qwraps2_markup = "markdown")

# Manually specified summary
crc.sum <-
  list("Sex" = 
         list("Male"       =  ~ n_perc0(Sex == 1, digits = 1),
              "Female"     =  ~ n_perc0(Sex == 2, digits = 1)),
       "Age at blood collection (years)" = 
         list("Mean"     =  ~ mean_sd(Age_Blood, digits = 1)),
       
       "Follow-up time to diagnosis (years)" =
         list("Mean"     =  ~ mean_sd(Tfollowup, digits = 1)),
       
       "Country" = 
         list("France"         = ~ n_perc0(Country == 1, digits = 1),
              "Italy"          = ~ n_perc0(Country == 2, digits = 1),
              "Spain"          = ~ n_perc0(Country == 3, digits = 1),
              "United Kingdom" = ~ n_perc0(Country == 4, digits = 1),
              "Netherlands"   = ~ n_perc0(Country == 5, digits = 1),
              #"Greece"        = ~ n_perc0(Country == 6, digits = 1),
              "Germany"        = ~ n_perc0(Country == 7, digits = 1),
              "Denmark"        = ~ n_perc0(Country == 9, digits = 1)),
       
       "Tumor site" =
         list("Proximal colon"  = ~ n_perc0(location == 1, na_rm = T, digits = 1),
              "Distal colon"    = ~ n_perc0(location == 2, na_rm = T, digits = 1),
              "Rectum"          = ~ n_perc0(location == 3, na_rm = T, digits = 1),
              "Other"           = ~ n_perc0(location == 4, na_rm = T, digits = 1),
              "Unknown"         = ~ n_perc0(is.na(location), digits = 1)),
       "Diagnosis with histological verification" =
         list("Yes"    = ~ n_perc0(Bdg1clrt_hist == 1, na_rm = T, digits = 1),
              "No"     = ~ n_perc0(Bdg1clrt_hist == 0, na_rm = T, digits = 1)),
       
       "Smoking status" =
         list("Non smoker"     = ~ n_perc0(Smoke_Stat == 1, digits = 1),
              "Never smoker"   = ~ n_perc0(Smoke_Stat == 2, digits = 1),
              "Smoker"         = ~ n_perc0(Smoke_Stat == 3, digits = 1)),
       #"Unknown"      =  ~ n_perc0(Smoke_Stat == 4, digits = 1)),
       
       "Height (cm)" =
         list("Mean (SD)" = ~ mean_sd(Height_C, digits = 1)),
       "BMI (kg/m2)" =
         list("Mean (SD)" = ~ mean_sd(Bmi_C, digits = 1)),
       "Waist circumference (cm)" =
         list("Mean (SD)" = ~ mean_sd(Waist_C, digits = 1)),
       "Total energy intake (kCal)" =
         list("Mean (SD)" = ~ mean_sd(Qe_Energy, na_rm = T, show_n = "never", digits = 0)),
       
       "Physical activity" =
         list("Mean (SD)" = ~ mean_sd(Pa_Mets, na_rm = T, show_n = "never", digits = 1)),
       "Alcohol intake (g/day)" =
         list("Mean (SD)" = ~ mean_sd(Qe_Alc, na_rm = T, show_n = "never", digits = 1)),
       
       "WCRF score" = 
         list("Mean (SD)" = ~ mean_sd(Wcrf_C_Cal, na_rm = T, show_n = "never")),
       "Fatty acid metabolic signature" = 
         list("Mean (SD)" = ~ mean_sd(comp1, na_rm = T, show_n = "never")),
       "Endogenous metabolic signature" = 
         list("Mean (SD)" = ~ mean_sd(comp2, na_rm = T, show_n = "never"))
  )

st1 <- summary_table(group_by(crc.sig, Cncr_Caco_Clrt), crc.sum)
print(st1, cnames = c("Controls", "Cases"))
# Copy and paste output into an Rmarkdown file and render to word/pdf etc
# Note: for some reason Tfollowup didn't work. Calculated mean and SD manually

# Tests for p-values for table
# McNemar and Wilcoxon signed rank test?

# Arrange df to pair cases and controls
df1 <- crc.sig %>% arrange(Cncr_Caco_Clrt, Match_Caseset.x)

# Case/control models: CRC A
t.test(df1$Height_C ~ df1$Cncr_Caco_Clrt, paired = T)$p.value
t.test(df1$Age_Blood ~ df1$Cncr_Caco_Clrt, paired = T)$p.value
t.test(df1$Bmi_C ~ df1$Cncr_Caco_Clrt, paired = T)$p.value
t.test(df1$Waist_C ~ df1$Cncr_Caco_Clrt, paired = T)$p.value
chisq.test(df1$Cncr_Caco_Clrt, df1$Smoke_Stat)$p.value


# Contain NAs
df1a <- df1 %>% group_by(Match_Caseset) %>% filter(!anyNA(Qe_Energy))
t.test(df1a$Qe_Energy ~ df1a$Cncr_Caco_Clrt, paired = T)$p.value

df1b <- df1 %>% group_by(Match_Caseset) %>% filter(!anyNA(Qe_Alc))
wilcox.test(df1b$Qe_Alc ~ df1b$Cncr_Caco_Clrt, paired = T)$p.value

df1c <- df1 %>% group_by(Match_Caseset) %>% filter(!anyNA(Pa_Mets))
t.test(df1c$Pa_Mets ~ df1c$Cncr_Caco_Clrt, paired = T)$p.value

df1d <- df1 %>% group_by(Match_Caseset) %>% filter(!anyNA(Wcrf_C_Cal))
t.test(df1d$Wcrf_C_Cal ~ df1d$Cncr_Caco_Clrt, paired = T)$p.value


df2 <- crc2a %>% arrange(Cncr_Caco_Clrt, Match_Caseset)

# Case/control models: CRC B
t.test(df2$Height_C ~ df2$Cncr_Caco_Clrt, paired = T)$p.value
t.test(df2$Age_Blood ~ df2$Cncr_Caco_Clrt, paired = T)$p.value
t.test(df2$Bmi_C ~ df2$Cncr_Caco_Clrt, paired = T)$p.value
t.test(df2$Waist_C ~ df2$Cncr_Caco_Clrt, paired = T)$p.value
chisq.test(df2$Cncr_Caco_Clrt, df2$Smoke_Stat)$p.value

# Contain NAs
df2a <- df2 %>% group_by(Match_Caseset) %>% filter(!anyNA(Qe_Energy))
t.test(df2a$Qe_Energy ~ df2a$Cncr_Caco_Clrt, paired = T)$p.value

df2b <- df2 %>% group_by(Match_Caseset) %>% filter(!anyNA(Qe_Alc))
wilcox.test(df2b$Qe_Alc ~ df2b$Cncr_Caco_Clrt, paired = T)$p.value

df2c <- df2 %>% group_by(Match_Caseset) %>% filter(!anyNA(Pa_Mets))
t.test(df2c$Pa_Mets ~ df2c$Cncr_Caco_Clrt, paired = T)$p.value

df2d <- df2 %>% group_by(Match_Caseset) %>% filter(!anyNA(Wcrf_C_Cal))
t.test(df2d$Wcrf_C_Cal ~ df2d$Cncr_Caco_Clrt, paired = T)$p.value

# Baseline char pool----
# Baseline characteristics for the pooled controls

source("CRC_prep_data.R")
library(haven)
fullepic <- read.csv("full_epic_main_vars.csv")
wcrf <- read_dta("wcrf_score.dta")

pooled1 <- ctrl %>%
  select(Study, Sex, Country, Age_Recr, Height_C, Bmi_C, Waist_C, Qe_Energy, 
         Pa_Mets, Smoke_Stat, Smoke_Intensity, Qe_Alc, Wcrf_C_Cal) %>% 
  mutate(Set = "Pooled controls endogenous")

pooled1$Study <- as.factor(pooled1$Study)

pooled2 <- fa.ctrl %>%
  left_join(fullepic, by = "Idepic", suffix = c("", "_1")) %>%
  select(Study = STUDY, Sex, Country, Age_Recr, Height_C, Bmi_C, Waist_C, Qe_Energy = QE_ENERGY, 
         Pa_Mets, Smoke_Stat, Smoke_Intensity, Qe_Alc = QE_ALC, Wcrf_C_Cal) %>% 
  mutate(Set = "Pooled controls fatty acids")

pooled2$Study <- as.factor(pooled2$Study)

library(qwraps2)
options(qwraps2_markup = "markdown")

# Manually specified
crc_sum <-
  list("Total subjects"   =
         list("N"   =   ~ n()),
       "Study" = 
         list("Breast"         = ~ n_perc0(Study == "Breast", digits = 1),
              "Kidney"          = ~ n_perc0(Study == "Kidney", digits = 1),
              "Ovary"          = ~ n_perc0(Study == "Ovary", digits = 1),
              "Pancreas" = ~ n_perc0(Study == "Pancreas", digits = 1),
              "Prostate"   = ~ n_perc0(Study == "Prostate", digits = 1),
              "Liver"        = ~ n_perc0(Study == "Liver", digits = 1)),
       "Sex" = 
         list("Male"   =  ~ n_perc0(Sex == 1, digits = 1),
              "Female"     =  ~ n_perc0(Sex == 2, digits = 1)),
       "Age at recruitment (years)" = 
         list("Mean"     =  ~ mean_sd(Age_Recr)),
       "Height (cm)" =
         list("Mean (SD)" = ~ mean_sd(Height_C)),
       "BMI (kg/m2)" =
         list("Mean (SD)" = ~ mean_sd(Bmi_C)),
       "Waist circumference (cm)" =
         list("Mean (SD)" = ~ mean_sd(Waist_C, na_rm = T)),
       "Total energy intake (kCal)" =
         list("Mean (SD)" = ~ mean_sd(Qe_Energy, na_rm = T, show_n = "never")),
       "Country" = 
         list("France"         = ~ n_perc0(Country == 1, digits = 1),
              "Italy"          = ~ n_perc0(Country == 2, digits = 1),
              "Spain"          = ~ n_perc0(Country == 3, digits = 1),
              "United Kingdom" = ~ n_perc0(Country == 4, digits = 1),
              "Netherlands"   = ~ n_perc0(Country == 5, digits = 1),
              #"Greece"        = ~ n_perc0(Country == 6, digits = 1),
              "Germany"        = ~ n_perc0(Country == 7, digits = 1),
              "Sweden"         = ~ n_perc0(Country == 8, digits = 1),
              "Denmark"        = ~ n_perc0(Country == 9, digits = 1),
              "Norway"         = ~ n_perc0(Country == "B", digits = 1)),
       "Physical activity" =
         list("Mean (SD)" = ~ mean_sd(Pa_Mets, na_rm = T, show_n = "never")),
       "Alcohol intake (g/day)" =
         list("Mean (SD)" = ~ mean_sd(Qe_Alc, na_rm = T, show_n = "never")),
       "Smoking status" =
         list("Non smoker"     = ~ n_perc0(Smoke_Stat == 1, digits = 1),
              "Never smoker"   = ~ n_perc0(Smoke_Stat == 2, digits = 1),
              "Smoker"         = ~ n_perc0(Smoke_Stat == 3, digits = 1)),
       #"Unknown"      =  ~ n_perc0(Smoke_Stat == 4, digits = 1)),
       "WCRF score" = 
         list("Mean (SD)" = ~ mean_sd(Wcrf_C_Cal, na_rm = T, show_n = "never"))
  )
st1 <- summary_table(pooled1, crc_sum)
st2 <- summary_table(pooled2, crc_sum)
both <- cbind(st1, st2)
print(both, cnames = c("Pooled controls endogenous", "Pooled controls fatty acids"))
