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
meta <- read_dta("clrt_caco.dta") %>% 
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


# Large case-control subset (from Jelena)------------

crc2 <- read_csv("biocrates_p150.csv") %>% mutate_at(vars(var.list), as.factor) %>%
  inner_join(wcrf, by = "Idepic") %>%
  filter(Country != 6)

# EPIC controls. First dataset, 3771 obs; updated November 2018 7191 obs
# Rename factor levels
ctrl <- read_dta("obes_metabo.dta") %>% mutate(Study = 
  fct_recode(Study, Colonrectum_S = "Colonrectum", Colonrectum_L = "Colonrectum (bbmri)")) %>%

# Split Batch_MetBio into two columns, then extract numeric variable
# Remove 1694 CRC controls
  separate(Batch_MetBio, into = c("batch", "rest")) %>%
  mutate(batch_no = as.numeric(flatten(str_extract_all(batch, "[0-9]+")))) %>% 
  filter(!(Study %in% c("Colonrectum_L", "Colonrectum_S")), Fasting_C == 2) %>%
  filter(Country != 6)

print(paste(nrow(ctrl), "fasted controls read"))
# 1799 fasted subjects left
# 1741 after removal of Greece

# Get common compounds between CC and controls and order
select.ctrl.cmpds <- function(crc, no.subset = F, cor.data = F){
  
  library(tidyverse)
  # Gets common compounds between CC studies and EPIC controls. Puts them in the same order.
  # Subset biocrates compounds
  controls <- ctrl %>% 
    select(matches("Acylcarn_|Aminoacid_|Biogenic_|Glyceroph_|Sphingo_|Sugars_"), -starts_with("Outdq"))
  
  zerocols <- apply(controls, 2, function(x) sum(x, na.rm = T)) != 0
  controls <- controls[, zerocols]
  colnames(controls) %>% length # 159 variables
  
  if(no.subset == T) return(controls)
  
  # ----
  
  print(paste("Subjects in case-control:", nrow(crc)))
  
  cmpds <- crc %>% 
    select(matches("Acylcarn_|Aminoacid_|Biogenic_|Glyceroph_|Sphingo_|Sugars_"), -starts_with("Outdq"))
  
  print(paste("No. of compounds:", ncol(cmpds)))
  
  # Convert variables to factors
  var.list <- c("Country", "Center", "Sex")
  crc <- crc %>% mutate_at(vars(var.list), as.factor)

  zerocols <- apply(cmpds, 2, function(x) sum(x, na.rm = T)) != 0
  cmpds <- cmpds[, zerocols]
  colnames(cmpds) %>% length # 147 variables
  
  # Get common compounds between controls and CRC CC and select compounds
  common.cols <- intersect(colnames(controls), colnames(cmpds)) # 146 cpds
  controls <- controls %>% select(one_of(common.cols))
  
  # Make dataset for correlation with FAs
  common.cols1 <- sort(common.cols)
  Biocrates <- crc %>% select(Idepic, one_of(common.cols1))
  
  if(cor.data == T) return(Biocrates)
  return(controls)
}
ctrlA <- select.ctrl.cmpds(crc1)
ctrlB <- select.ctrl.cmpds(crc2)
ctrls <- select.ctrl.cmpds(no.subset = T)

# Get compounds common to all 3 sets for comparison and subset df
common.all <- intersect(colnames(ctrlA), colnames(ctrlB))
ctrls0 <- select(ctrls, one_of(common.all))


# Fatty acids-------------

# Gets common compounds between CC and EPIC controls and puts them in the same order.

# Get CRC dataset from Elom and join WCRF scores. Convert categorical co-variates to factors
var.list <- c("L_School", "Smoke_Stat")
CRCfa1 <- read_dta("Database_Fatty acids.dta") %>% 
  left_join(wcrf, by = "Idepic") %>% 
  mutate_at(vars(var.list), as.factor) %>%
  filter(Country != 6)

# Get dataset for PLS modelling (all EPIC controls). Exclude compounds with many missings
# Note: new version from Carine received 18/11/2018 with technical covariates
fa.ctrl <- readRDS("FA_WCRF_scores1.rds") %>% filter(STUDY != "Colorectum") %>%
  filter(Country != 6)
fa.ctrl$N_Serie <- as.numeric(fa.ctrl$N_Serie)

# categorical variables to factors
var.list <- c("Country", "Center", "STUDY", "LABO")
fa.ctrl <- fa.ctrl %>% mutate_at(vars(var.list), as.factor)

# Subset concentrations for CRC and controls
CRCfa <- CRCfa1  %>% select(P14_0 : PCLA_9t_11c) 
concs <- fa.ctrl %>% select(P14_0 : PCLA_9t_11c, -P24_0, -P20_0)
common.cols <- intersect(colnames(concs), colnames(CRCfa))
  
# return(list(fa.ctrl, common.cols))



