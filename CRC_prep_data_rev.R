# For revised submission to CGH. CRC1/A and CRC2/B are merged into one
library(tidyverse)
library(haven)
library(lubridate)

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