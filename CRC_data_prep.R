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
  meta <- read_dta("clrt_caco.dta") %>% 
    select(-Match_Caseset, -Country, -Center, -Cncr_Caco_Clrt) %>%
    distinct(Idepic, .keep_all = T)
  
  # meta.dedup <- meta[!duplicated(meta["Idepic"]), ]
  
  # Converted to dta
  # write_dta(data, "clrt_caco_metabo.dta")
  
  # Subset biocrates concentrations
  concs <- crc %>% 
    select(matches("Acylcarn_|Aminoacid_|Biogenic_|Glyceroph_|Sphingo_|Sugars_"), -starts_with("Outdq"))
  
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
  
  # Recodes batch and filters fasted samples
  
  if(data == "new") {
  
    library(haven)
    ctrl <- read_dta("obes_metabo.dta")
    
    # Convert Study variable to factor and rename levels
    ctrl$Study <- fct_recode(ctrl$Study, Colonrectum_S = "Colonrectum", Colonrectum_L = "Colonrectum (bbmri)")
    
    # Split Batch_MetBio into two columns, then extract numeric variable
    # Remove 1694 CRC controls
    ctrl <- ctrl %>% 
      separate(Batch_MetBio, into = c("batch", "rest")) %>%
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
    
    ctrl <- ctrl %>% 
      separate(Batch_MetBio, into = c("batch", "rest")) %>%
      mutate(batch_no = as.numeric(flatten(str_extract_all(batch, "[0-9]+")))) %>% 
      filter(Fasting_C == "Yes", Study != "Colonrectum")
  
  }
  
}
ctrl <- prepctrl(data = "new")

get.common.Bioc <- function(fasting = T){
  
  library(tidyverse)
  # ets common compounds between CC studies and EPIC controls. Puts them in the same order.G
  # First dataset, 3771 obs; updated November 2018 7191 obs
  print(paste(nrow(ctrl), "Controls read"))
  
  # Subset biocrates compounds
  controls <- ctrl %>% 
    select(matches("Acylcarn_|Aminoacid_|Biogenic_|Glyceroph_|Sphingo_|Sugars_"), -starts_with("Outdq"))
  
  zerocols <- apply(controls, 2, function(x) sum(x, na.rm = T)) != 0
  controls <- controls[, zerocols]
  
  colnames(controls) %>% length # 147 variables
  
  # Large CRC metabolomics subset
  
  # (from Jelena, ~ 1200 case-control pairs), call it crc2
  # Join scores and filter out non-fasted samples
  crc2 <- if(fasting == T) crc2 %>% filter(Fasting_C == 2) else crc2
  print(paste("Subjects in larger case-control: ", nrow(crc2)))
  # Variables were converted to factors in CRC_data_prep.R
  
  setA <- crc2 %>% 
    select(matches("Acylcarn_|Aminoacid_|Biogenic_|Glyceroph_|Sphingo_|Sugars_"), -starts_with("Outdq"))
  zerocols1 <- apply(setA, 2, function(x) sum(x, na.rm = T)) != 0
  setA <- setA[, zerocols1]
  colnames(setA) %>% length
  # 163 variables
  
  # Get common compounds between controls and CRC CC
  common_cols <- intersect(colnames(controls), colnames(setA))   # 135 compounds in common
  
  # Small CRC metabolomics subset
  # (from Bertrand, ~ 490 case-control pairs), call it crc1
  print(paste("Subjects in smaller case-control:", nrow(crc1)))
  
  setB <- crc1 %>% 
    select(matches("Acylcarn_|Aminoacid_|Biogenic_|Glyceroph_|Sphingo_|Sugars_"), -starts_with("Outdq"))
  colnames(setB) %>% length
  
  # Convert variables to factors
  var.list <- c("Country", "Center", "Sex")
  crc1 <- crc1 %>% mutate_at(vars(var.list), as.factor)
  #common_cols <- intersect(colnames(controls), colnames(setB))   # 145 compounds in common
  
  # Get common cols between all three datasets (controls, small CC, big CC)
  common_cols2 <- intersect(common_cols, colnames(setB))   # Glyceroph_Lysopc_A_C24_0 is removed
  
  # subset and reorder both datasets to get the same 126 compounds in the same order
  controls <- controls %>% select(one_of(common_cols2))
  setA <- setA %>% select(one_of(common_cols2))
  #setB <- setB %>% select(one_of(common_cols2))
  
  # check colnames are the same for both sets
  identical(colnames(controls), colnames(setA))
  return(controls)
}
controls <- get.common.Bioc()

get.common.FAs  <- function(){
  
  # Gets common compounds between CC and EPIC controls and puts them in the same order.
  # Outputs the controls dataset and ordered list of compounds
  library(haven)
  library(tidyverse)
  library(Amelia)
  
  # Get WCRF scores
  wcrf <- read_dta("Wcrf_Score.dta") %>% select(Idepic, Wcrf_C_Cal)
  
  # Get CRC dataset from Elom and join WCRF scores. Convert categorical co-variates to factors
  var.list <- c("L_School", "Smoke_Stat")
  CRCfa1 <- read_dta("Database_Fatty acids.dta") %>% 
    left_join(wcrf, by = "Idepic") %>% 
    mutate_at(vars(var.list), as.factor)
  
  # Get dataset for PLS modelling (all EPIC controls). Exclude compounds with many missings
  # Note: new version from Carine received 18/11/2018 with technical covariates
  fa.ctrl <- readRDS("FA_WCRF_scores1.rds") %>% filter(STUDY != "Colorectum")
  fa.ctrl$N_Serie <- as.numeric(fa.ctrl$N_Serie)
  
  # categorical variables to factors
  var.list <- c("Country", "Center", "STUDY", "LABO")
  fa.ctrl <- fa.ctrl %>% mutate_at(vars(var.list), as.factor)
  
  # Subset concentrations for CRC and controls
  CRCfa <- CRCfa1  %>% select(P14_0 : PCLA_9t_11c) 
  concs <- fa.ctrl %>% select(P14_0 : PCLA_9t_11c, -P24_0, -P20_0)

  #missmap(concs)
  #missmap(CRCfa)
  
  length(colnames(concs))
  length(colnames(CRCfa))
  common_cols <- intersect(colnames(concs), colnames(CRCfa))
  
  return(list(fa.ctrl, common_cols))
}
output <- get.common.FAs()





