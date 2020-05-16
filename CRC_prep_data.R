# Preparation of CRC case-control datasets, and control only datasets 
# stored at \\inti\NME\EPIC_Projects\Epic_Colonrectum\Nested_CaCo_Study\2016
# Missings are already imputed
library(tidyverse)
library(haven)

# Metadata setup ----------

# Read metadata for whole CRC case-control 
# WARNING: do not use case-control status from this dataset. Use metabolomics datasets only.
# Remove duplicated Idepics (with dplyr or base). Also get follow up time and colorectal site
var.list <- c("Country", "Center", "Sex", "Match_Caseset", "L_School", #"Smoke_Int", 
              "Smoke_Stat", "Smoke_Intensity", "Fasting_C", "Menopause", "Phase_Mnscycle")

meta <- read_dta("clrt_caco.dta") %>% 
  # set D_Dgclrt of controls to that of corresponding cases, calculate followup time and 
  # and get cancer site variables
  mutate(Tfollowup.days = D_Dgclrt - D_Bld_Coll, Tfollowup = Tfollowup.days/365.25, 
         location = case_when(
    Case_Mal_Colon_Prox == 1 ~ 1, Case_Mal_Colon_Dist == 1 ~ 2,
    Case_Mal_Colon_Nos  == 1 ~ 4, Case_Mal_Rectum     == 1 ~ 3)) %>%
  group_by(Match_Caseset) %>% fill(c(D_Dgclrt, location), .direction = "downup") %>% ungroup() %>%
  select(-Match_Caseset, -Cncr_Caco_Clrt) %>%
  distinct(Idepic, .keep_all = T)

# Combine smoking intensity factor levels
#meta$Smoke_Int <- fct_collapse(meta$Smoke_Intensity, Other = c("8", "9", "10"))


# Small case-control subset----------
# Use glutamate to correctly subset biocrates data, join metadata. Update 3/3/2020: remove Greece

crc1 <- read_sas("clrt_caco_metabo.sas7bdat") %>% filter(!is.na(Aminoacid_Glu)) %>%
        left_join(meta, by = "Idepic", suffix = c("_1", "")) %>%
        mutate_at(vars(var.list), as.factor) %>% 
        mutate(Smoke_Int = fct_collapse(Smoke_Intensity, Other = c("8", "9", "10"))) %>% 
        filter(Country != 6)

# Get colon cancer only (ungroup to stop Match_Caseset from being readded later)
colon1 <- crc1 %>% group_by(Match_Caseset) %>% 
  filter(max(location, na.rm = T) == 1 | max(location, na.rm = T) == 2) %>% ungroup(Match_Caseset)

# Subset male or female
#crc1.ma <- crc1 %>% filter(Sex == 1)
#crc1.fe <- crc1 %>% filter(Sex == 2)


# Large case-control subset (from Jelena)------------

library(lubridate)
crc2 <- read_csv("biocrates_p150.csv") %>% 
  select(Match_Caseset, Cncr_Caco_Clrt, ends_with("Idepic"), matches("(carn|oacid|genic|roph|ingo|Sugars)[_]"), -contains("tdq")) %>%
  inner_join(meta, by = "Idepic") %>% mutate_at(vars(var.list), as.factor) %>% 
  mutate(Smoke_Int = fct_collapse(Smoke_Intensity, Other = c("8", "9", "10"))) %>%
  filter(Country != 6)

# Get colon cancer or rectal cancer subsets
colon2 <- crc2 %>% group_by(Match_Caseset) %>% 
  filter(max(location, na.rm = T) == 1 | max(location, na.rm = T) == 2) %>% ungroup(Match_Caseset)

rectal2 <- crc2 %>% group_by(Match_Caseset) %>% 
  filter(max(location, na.rm = T) == 3) %>% ungroup(Match_Caseset)

# Subset male or female
#crc2.ma <- crc1 %>% filter(Sex == 1)
#crc2.fe <- crc1 %>% filter(Sex == 2)


# EPIC pooled controls--------
# First dataset, 3771 obs; updated November 2018 7191 obs
# Rename factor levels, split Batch_MetBio into 2 cols, extract numeric variable,
# Remove 1694 CRC controls

ctrl <- read_dta("obes_metabo.dta") %>% mutate(Study = 
  fct_recode(Study, Colonrectum_S = "Colonrectum", Colonrectum_L = "Colonrectum (bbmri)")) %>%
  separate(Batch_MetBio, into = c("batch", "rest")) %>%
  mutate(batch_no = as.numeric(flatten(str_extract_all(batch, "[0-9]+")))) %>% 
  filter(!(Study %in% c("Colonrectum_L", "Colonrectum_S")), Fasting_C == 2) %>%
  filter(Country != 6)

#print(paste(nrow(ctrl), "fasted controls read"))
# 1799 fasted subjects left, 1741 after removal of Greece

# Get common compounds between CC and controls and get subset of compounds for each signature
select.ctrl.cmpds <- function(datalist, cor.data = F){
  
  # Sub-function to subset compounds only and remove zero columns
  get.cmpds <- function(dat) { 
    cmpds <- dat %>% select(matches("(carn|oacid|genic|roph|ingo|Sugars)[_]"), -contains("tdq")) %>%
      select_if(~ sum(., na.rm = T) != 0)
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

# Get compounds common to all 3 sets for comparison and subset df
common.all <- intersect(colnames(ctrlA), colnames(ctrlB))
ctrls0 <- select(ctrls, one_of(common.all))


# Fatty acids CRC dataset (from Elom)---------------

# Gets common compounds between CC and EPIC controls and puts them in the same order.

# Get CRC dataset from Elom and join WCRF scores. Convert categorical co-variates to factors
crc1fa <- read_dta("Database_Fatty acids.dta") %>% 
  left_join(meta, by = "Idepic", suffix = c("", "_1")) %>%
  mutate_at(vars(var.list), as.factor) %>%
  mutate(Smoke_Int = fct_collapse(Smoke_Intensity, Other = c("8", "9", "10"))) %>%
  filter(Country != 6)

#crc1faf <- crc1fa %>% filter(Sex == 2)

# Get dataset for PLS modelling (all EPIC controls). Exclude compounds with many missings
# Note: new version from Carine received 18/11/2018 with technical covariates
fa.ctrl <- readRDS("FA_WCRF_scores1.rds") %>% filter(STUDY != "Colorectum") %>% filter(Country != 6)
fa.ctrl$N_Serie <- as.numeric(fa.ctrl$N_Serie)

# Get sex of participants
epic.vars <- read.csv("full_epic_sex.csv")
fa.ctrl <- left_join(fa.ctrl, epic.vars, by = "Idepic")
#fa.ctrl.f <- fa.ctrl %>% filter(Sex == 2)
# Not used because 118 and 4121 M and F respectively

# Convert variables to factors
var.list <- c("Country", "Center", "STUDY", "LABO")
fa.ctrl <- fa.ctrl %>% mutate_at(vars(var.list), as.factor)

# Subset concentrations for CRC and controls
crcfa <- crc1fa  %>% select(P14_0 : PCLA_9t_11c) 
concs <- fa.ctrl %>% select(P14_0 : PCLA_9t_11c, -P24_0, -P20_0)
common.cols <- intersect(colnames(concs), colnames(crcfa))

#concs.f <- fa.ctrl.f %>% select(P14_0 : PCLA_9t_11c, -P24_0, -P20_0)

# Subset only controls from FAs and CRC A
CRCfa.ctrl <- crc1fa %>% filter(Cncr_Caco_Clrt == 0) %>% select(Idepic, P14_0 : PCLA_9t_11c) 

# Number of control profiles for biocrates and fatty acids
nrow(ctrl)
nrow(fa.ctrl)

# Number of control subjects , biocrates and fatty acids combined
length(intersect(ctrl$Idepic, fa.ctrl$Idepic))
#intersect(ctrl$Idepic, fa.ctrl$Idepic)

# Old subset for sex-specific signature
#ctrl.m <- ctrl %>% filter(Sex == 1)
#ctrl.f <- ctrl %>% filter(Sex == 2)


# Data for sex-specific signatures
#ctrlAm <- select.ctrl.cmpds(list(ctrl.m, crc1.ma))
#ctrlAf <- select.ctrl.cmpds(list(ctrl.f, crc1.fe))
#ctrlBm <- select.ctrl.cmpds(list(ctrl.m, crc2.ma))
#ctrlBf <- select.ctrl.cmpds(list(ctrl.f, crc2.fe))
