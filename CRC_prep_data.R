# First part of CRC_prep_data_rev.R that prepares case-control studies
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
crc1 <- read_sas("clrt_caco_metabo.sas7bdat") %>% filter(!is.na(Aminoacid_Glu)) %>%
  left_join(meta, by = "Idepic", suffix = c("_1", "")) %>%
  mutate_at(vars(var.list), as.factor) %>% 
  mutate(Smoke_Int = fct_collapse(Smoke_Intensity, Other = c("8", "9", "10"))) %>% 
  filter(Country != 6) # Greece removed

# Get colon cancer only (ungroup to stop Match_Caseset from being readded later)
colon1 <- crc1 %>% filter(location == 1 | location  == 2)

# Subsites
rectal1 <- crc1 %>% group_by(Match_Caseset) %>% filter(max(location, na.rm = T) == 3) %>% ungroup(Match_Caseset)
prox1 <- crc1 %>% group_by(Match_Caseset) %>% filter(max(location, na.rm = T) == 1) %>% ungroup(Match_Caseset)
dist1 <- crc1 %>% group_by(Match_Caseset) %>% filter(max(location, na.rm = T) == 2) %>% ungroup(Match_Caseset)

# Subset male, female, cases diagnosed after 2 years only
crc1m <- crc1 %>% filter(Sex == 1)
crc1f <- crc1 %>% filter(Sex == 2)
crc1t <- crc1 %>% group_by(Match_Caseset) %>% filter(max(Tfollowup, na.rm = T) > 2) %>% ungroup()

# Large case-control subset (from Jelena)
crc2 <- read_csv("biocrates_p150.csv") %>% 
  select(Match_Caseset, Cncr_Caco_Clrt, 
         ends_with("Idepic"), matches("(carn|oacid|genic|roph|ingo|Sugars)[_]"), -contains("tdq")) %>%
  inner_join(meta, by = "Idepic") %>% mutate_at(vars(var.list), as.factor) %>% 
  mutate(Smoke_Int = fct_collapse(Smoke_Intensity, Other = c("8", "9", "10"))) %>%
  filter(Country != 6)

# Get colon cancer or rectal cancer subsets
colon2 <- crc2 %>% group_by(Match_Caseset) %>% 
  filter(max(location, na.rm = T) == 1 | max(location, na.rm = T) == 2) %>% ungroup(Match_Caseset)

rectal2 <- crc2 %>% group_by(Match_Caseset) %>% filter(max(location, na.rm = T) == 3) %>% ungroup(Match_Caseset)
prox2 <- crc2 %>% group_by(Match_Caseset) %>% filter(max(location, na.rm = T) == 1) %>% ungroup(Match_Caseset)
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

