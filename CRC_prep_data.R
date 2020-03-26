# Preparation of CRC case-control datasets, and control only datasets 
# stored at \\inti\NME\EPIC_Projects\Epic_Colonrectum\Nested_CaCo_Study\2016
# Missings are already imputed
library(tidyverse)
library(haven)

# Small case-control subset----------
# Update 3/3/2020: remove Greece

# have to subset subjects with Biocrates data
crc <- read_sas("clrt_caco_metabo.sas7bdat")

# Metadata and WCRF scores (keep on local drive because big file)
wcrf <- read_dta("Wcrf_Score.dta") %>% select(ends_with("_Cal"), Idepic)

# Remove duplicated Idepics (with dplyr or base)
# Also get follow up time and colorectal site
meta <- read_dta("clrt_caco.dta") %>% 
  mutate(Tfollowup.days = D_Dgclrt - D_Bld_Coll, Tfollowup = Tfollowup.days/365.25, location = case_when(
    Case_Mal_Colon_Prox == 1 ~ 1,
    Case_Mal_Colon_Dist == 1 ~ 2,
    Case_Mal_Colon_Nos  == 1 ~ 4,
    Case_Mal_Rectum     == 1 ~ 3)) %>%
  select(-Match_Caseset, -Country, -Center, -Cncr_Caco_Clrt) %>%
  distinct(Idepic, .keep_all = T)

# Subset biocrates variables
concs <- crc %>% 
  select(matches("Acylcarn_|Aminoacid_|Biogenic_|Glyceroph_|Sphingo_|Sugars_"), -starts_with("Outdq"))

# Subset biocrates data by taking !is.na > 0
biocrates <- apply(concs, 1, function(x) sum(!is.na(x)) > 0)
crc1 <- crc[biocrates, ]

var.list <- c("Country", "Center", "Sex", "Match_Caseset", "Smoke_Stat", "L_School")
crc1 <- crc1 %>% inner_join(meta, by = "Idepic") %>% mutate_at(vars(var.list), as.factor) %>%
  filter(Country != 6)
  #group_by(Match_Caseset) %>% filter(n() == 2)

# Get colon cancer only (ungroup to stop Match_Caseset from being readded later)
colon1 <- crc1 %>% group_by(Match_Caseset) %>% 
  filter(max(location, na.rm = T) == 1 | max(location, na.rm = T) == 2) %>% ungroup(Match_Caseset)

# Subset male or female
crc1.ma <- crc1 %>% filter(Sex == 1)
crc1.fe <- crc1 %>% filter(Sex == 2)

# Large case-control subset (from Jelena)------------

library(lubridate)
crc2 <- read_csv("biocrates_p150.csv")
crc2$D_Bld_Coll <- dmy(crc2$D_Bld_Coll)
crc2$D_Dgclrt <- dmy(crc2$D_Dgclrt)

crc2 <- crc2 %>% mutate_at(vars(var.list), as.factor) %>%
  mutate_at(vars(starts_with("D_")), dmy) %>%
  mutate(Tfollowup.days = D_Dgclrt - D_Bld_Coll, Tfollowup = Tfollowup.days/365.25, location = case_when(
    Case_Mal_Colon_Prox == 1 ~ 1,
    Case_Mal_Colon_Dist == 1 ~ 2,
    Case_Mal_Colon_Nos  == 1 ~ 4,
    Case_Mal_Rectum     == 1 ~ 3)) %>%
  inner_join(wcrf, by = "Idepic") %>%
  filter(Country != 6)

# Get colon cancer or rectal cancer only
colon2 <- crc2 %>% group_by(Match_Caseset) %>% 
  filter(max(location, na.rm = T) == 1 | max(location, na.rm = T) == 2) %>% ungroup(Match_Caseset)

rectal2 <- crc2 %>% group_by(Match_Caseset) %>% 
  filter(max(location, na.rm = T) == 3) %>% ungroup(Match_Caseset)

# Subset male or female
crc2.ma <- crc1 %>% filter(Sex == 1)
crc2.fe <- crc1 %>% filter(Sex == 2)


# EPIC controls. First dataset, 3771 obs; updated November 2018 7191 obs
# Rename factor levels, split Batch_MetBio into 2 cols, extract numeric variable,
# Remove 1694 CRC controls

ctrl <- read_dta("obes_metabo.dta") %>% mutate(Study = 
  fct_recode(Study, Colonrectum_S = "Colonrectum", Colonrectum_L = "Colonrectum (bbmri)")) %>%
  separate(Batch_MetBio, into = c("batch", "rest")) %>%
  mutate(batch_no = as.numeric(flatten(str_extract_all(batch, "[0-9]+")))) %>% 
  filter(!(Study %in% c("Colonrectum_L", "Colonrectum_S")), Fasting_C == 2) %>%
  filter(Country != 6)

# Subgroup for sex-specific signature
ctrl.m <- ctrl %>% filter(Sex == 1)
ctrl.f <- ctrl %>% filter(Sex == 2)

print(paste(nrow(ctrl), "fasted controls read"))
# 1799 fasted subjects left
# 1741 after removal of Greece

# Get common compounds between CC and controls and get subset of compounds for each signature
select.ctrl.cmpds <- function(datalist, cor.data = F){

  # Sub-function to subset compounds only and remove zero columns
  get.cmpds <- function(dat) { 
    cmpds <- dat %>% 
      select(matches("Acylcarn_|Aminoacid_|Biogenic_|Glyceroph_|Sphingo_|Sugars_"), -starts_with("Outdq"))
    zerocols <- apply(cmpds, 2, function(x) sum(x, na.rm = T)) != 0
    cmpds <- cmpds[, zerocols]
  }
  cmpdlist <- lapply(datalist, get.cmpds)

  if(length(datalist) == 1) { 
    return(cmpdlist[[1]])
  } else {
    # Get common compounds between controls and CRC CC and select compounds
    common.cols <- intersect(colnames(cmpdlist[[1]]), colnames(cmpdlist[[2]]))
  }
  
  # Get data subset by common cols
  output <- cmpdlist[[1]] %>% select(one_of(common.cols))
  if(cor.data == F) return(output)
  
  # Make dataset for correlation with FAs
  common.cols1 <- sort(common.cols)
  dat.sorted.cols <- datalist[[2]] %>% select(Idepic, one_of(common.cols1))
}

ctrlA <- select.ctrl.cmpds(list(ctrl, crc1))
ctrlB <- select.ctrl.cmpds(list(ctrl, crc2))
ctrls <- select.ctrl.cmpds(list(ctrl))

# Data for sex-specific signatures
ctrlAm <- select.ctrl.cmpds(list(ctrl.m, crc1.ma))
ctrlAf <- select.ctrl.cmpds(list(ctrl.f, crc1.fe))
ctrlBm <- select.ctrl.cmpds(list(ctrl.m, crc2.ma))
ctrlBf <- select.ctrl.cmpds(list(ctrl.f, crc2.fe))

# Get compounds common to all 3 sets for comparison and subset df
common.all <- intersect(colnames(ctrlA), colnames(ctrlB))
ctrls0 <- select(ctrls, one_of(common.all))

# Fatty acids-------------

# Gets common compounds between CC and EPIC controls and puts them in the same order.

# Get CRC dataset from Elom and join WCRF scores. Convert categorical co-variates to factors
var.list <- c("L_School", "Smoke_Stat")
crc1fa <- read_dta("Database_Fatty acids.dta") %>% 
  left_join(wcrf, by = "Idepic") %>% 
  mutate_at(vars(var.list), as.factor) %>%
  filter(Country != 6)

crc1faf <- crc1fa %>% filter(Sex == 2)

# Get dataset for PLS modelling (all EPIC controls). Exclude compounds with many missings
# Note: new version from Carine received 18/11/2018 with technical covariates
fa.ctrl <- readRDS("FA_WCRF_scores1.rds") %>% filter(STUDY != "Colorectum") %>%
  filter(Country != 6)
fa.ctrl$N_Serie <- as.numeric(fa.ctrl$N_Serie)

# Get sex of participants
epic.vars <- read.csv("full_epic_sex.csv")
fa.ctrl <- left_join(fa.ctrl, epic.vars, by = "Idepic")
fa.ctrl.f <- fa.ctrl %>% filter(Sex == 2)
# Not used because 118 and 4121 M and F respectively

# categorical variables to factors
var.list <- c("Country", "Center", "STUDY", "LABO")
fa.ctrl <- fa.ctrl %>% mutate_at(vars(var.list), as.factor)

# Subset concentrations for CRC and controls
crcfa <- crc1fa  %>% select(P14_0 : PCLA_9t_11c) 
concs <- fa.ctrl %>% select(P14_0 : PCLA_9t_11c, -P24_0, -P20_0)
common.cols <- intersect(colnames(concs), colnames(crcfa))

concs.f <- fa.ctrl.f %>% select(P14_0 : PCLA_9t_11c, -P24_0, -P20_0)

# Subset only controls from FAs and CRC A
CRCfa.ctrl <- crc1fa %>% filter(Cncr_Caco_Clrt == 0) %>% select(Idepic, P14_0 : PCLA_9t_11c) 

# Remove unneeded data from workspace
rm(wcrf)
rm(meta)
rm(crc)

# Number of control profiles for biocrates and fatty acids
nrow(ctrl)
nrow(fa.ctrl)

# Number of control subjects , biocrates and fatty acids combined
length(intersect(ctrl$Idepic, fa.ctrl$Idepic))
intersect(ctrl$Idepic, fa.ctrl$Idepic)
