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




