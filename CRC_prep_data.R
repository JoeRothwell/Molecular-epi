# First part of CRC_prep_data_rev.R that prepares case-control studies
library(tidyverse)
library(haven)
library(lubridate)

# Read metadata for whole CRC case-control 
# WARNING: do not use case-control status from this dataset. Use metabolomics datasets only.
# Remove duplicated Idepics (with dplyr or base). Also get follow up time and colorectal site
var.list <- c("Country", "Center", "Sex", "Match_Caseset", "L_School", #"Smoke_Int", 
              "Smoke_Stat", "Smoke_Intensity", "Fasting_C", "Menopause", "Phase_Mnscycle",
              "Alc_Drinker", "Pa_Total")

# set D_Dgclrt of controls to that of corresponding cases, calculate followup time and 
# and get colorectal subsite variables
meta <- read_dta("clrt_caco.dta") %>% 
  mutate(Tfollowup.days = D_Dgclrt - D_Bld_Coll, Tfollowup = Tfollowup.days/365.25, 
         location = case_when(
           Case_Mal_Colon_Prox == 1 ~ 1, Case_Mal_Colon_Dist == 1 ~ 2,
           Case_Mal_Colon_Nos  == 1 ~ 4, Case_Mal_Rectum     == 1 ~ 3)) %>%
  group_by(Match_Caseset) %>% fill(c(D_Dgclrt, location), .direction = "downup") %>% 
  ungroup() %>% select(-Match_Caseset, -Cncr_Caco_Clrt) %>%
  distinct(Idepic, .keep_all = T)

# Small case-control subset (p180)
# Use Batch_MetBio to correctly subset biocrates data, join metadata.
# 496 cases, 492 controls in this dataset. Delete if only 1 in the caseset.
crc1 <- read_sas("clrt_caco_metabo.sas7bdat") %>% 
  filter(!is.na(Batch_MetBio)) %>%
  left_join(meta, by = "Idepic", suffix = c("_1", "")) %>%
  mutate_at(vars(var.list), as.factor) %>% 
  mutate(Smoke_Int = fct_collapse(Smoke_Intensity, Other = c("8", "9", "10"))) %>% 
  group_by(Match_Caseset) %>% filter(n() == 2) %>% ungroup()
  #filter(Country != 6) # Greece removed

# 980 subjects left, checked 490 cases 490 controls

# Add categorical BMI
crc1$Bmi_Cat<- as.factor(cut(crc1$Bmi_C, c(0,25,30,99), labels=FALSE, right=FALSE))

# Amino acids study
# Get colon cancer only, fasting status (2) only
#colon1 <- crc1 %>% filter(!location %in% 3 & Fasting_C == 2)

# Get colon cancer by Jelena's list of IDs and join to data to leave 740 subjects
p180ids <- read.csv("p180_ids.csv")
colon1 <- crc1 %>% inner_join(p180ids, by = c("Idepic" = "ids_p180")) %>% filter(Country != 6)

# Subsites
rectal1 <- crc1 %>% group_by(Match_Caseset) %>% filter(max(location, na.rm = T) == 3) %>% ungroup(Match_Caseset)
prox1 <- crc1 %>% group_by(Match_Caseset) %>% filter(max(location, na.rm = T) == 1) %>% ungroup(Match_Caseset)
dist1 <- crc1 %>% group_by(Match_Caseset) %>% filter(max(location, na.rm = T) == 2) %>% ungroup(Match_Caseset)

# Subset male, female, cases diagnosed after 2 years only
crc1m <- crc1 %>% filter(Sex == 1)
crc1f <- crc1 %>% filter(Sex == 2)
crc1t <- crc1 %>% group_by(Match_Caseset) %>% filter(max(Tfollowup, na.rm = T) > 2) %>% ungroup()

# Large case-control subset (p150 from Jelena)
# 1185 cases, 1185 controls
crc2 <- read_csv("biocrates_p150.csv") %>% 
  select(Match_Caseset, Cncr_Caco_Clrt, 
         ends_with("Idepic"), matches("(carn|oacid|genic|roph|ingo|Sugars)[_]"), -contains("tdq")) %>%
  inner_join(meta, by = "Idepic") %>% 
  mutate_at(vars(var.list), as.factor) %>% 
  mutate(Smoke_Int = fct_collapse(Smoke_Intensity, Other = c("8", "9", "10"))) %>%
  group_by(Match_Caseset) %>% filter(n() == 2) %>% ungroup()
  #filter(Country != 6)

crc2$Bmi_Cat<- as.factor(cut(crc2$Bmi_C, c(0,25,30,99), labels=FALSE, right=FALSE))

# Get colon cancer for amino acids study
#colon2 <- crc2 %>% filter(!location %in% 3 & Fasting_C == 2) #%>% 
  #group_by(Match_Caseset) %>% filter(n() == 2) %>% ungroup()

# Get colon cancer by Jelena's list of IDs
p150ids <- read.csv("p150_ids.csv")
colon2 <- crc2 %>% inner_join(p150ids, by = c("Idepic" = "ids_p150")) %>% filter(Country != 6)

#Checks
#table(colon2$Cncr_Caco_Clrt)
#intersect(colon2$Idepic, as.character(p150ids$ids_p150))
#unique(droplevels(colon2$Match_Caseset))

rectal2 <- crc2 %>% group_by(Match_Caseset) %>% filter(max(location, na.rm = T) == 3) %>% ungroup(Match_Caseset)
prox2 <- crc2 %>% group_by(Match_Caseset) %>% filter(max(location, na.rm = T) == 1) %>% ungroup(Match_Caseset)
dist2 <- crc2 %>% group_by(Match_Caseset) %>% filter( max(location, na.rm = T) == 2) %>% ungroup(Match_Caseset)

# Get cases diagnosed after 2 years only
crc2t <- crc2 %>% group_by(Match_Caseset) %>% filter(mean(Tfollowup, na.rm = T) > 2) %>% ungroup()

# Subset male or female
crc2m <- crc2 %>% filter(Sex == 1)
crc2f <- crc2 %>% filter(Sex == 2)

### 11 september 2020, reviewers revisions for CGH ###
# Merge crc1 and crc2 to make complete dataset (now moved to CRC_prep_data_rev)
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

