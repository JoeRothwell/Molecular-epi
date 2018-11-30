# Preparation of CRC case-control datasets

# Read in small CRC case control dataset (1) and subset the 146 biocrates compounds
# stored at \\inti\NME\EPIC_Projects\Epic_Colonrectum\Nested_CaCo_Study\2016
# Missings are already imputed

prepcrc1 <- function(){

  library(tidyverse)
  library(haven)
  crc <- read_sas("clrt_caco_metabo.sas7bdat")
  
  # Metadata and WCRF scores (keep on local drive because big file)
  wcrf <- read_dta("Wcrf_Score.dta") %>% select(Idepic, Wcrf_C_Cal)
  
  # Remove duplicated Idepics (with dplyr or base)
  meta <- read_dta("D:/clrt_caco.dta") %>% select(-Match_Caseset, -Country, -Center, -Cncr_Caco_Clrt) %>%
    distinct(Idepic, .keep_all = T)
  
  # meta.dedup <- meta[!duplicated(meta["Idepic"]), ]
  
  # Converted to dta
  # write_dta(data, "clrt_caco_metabo.dta")
  
  # Subset biocrates data by taking !is.na > 0
  concs <- crc %>% select(Acylcarn_C0 : Sugars_H1, -Batch_MetBio)
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

library(haven)
prepcontrol <- function(){
  
  ctrl <- read_dta("obes_metabo.dta")
  
  # Convert Study variable to factor and rename levels
  ctrl$Study <- fct_recode(ctrl$Study, Colonrectum_S = "Colonrectum", Colonrectum_L = "Colonrectum (bbmri)")
  
  # Split Batch_MetBio into two columns, then extract numeric variable
  # Remove 1694 CRC controls
  ctrl <- ctrl %>% separate(Batch_MetBio, into = c("batch", "rest")) %>%
    mutate(batch_no = as.numeric(flatten(str_extract_all(batch, "[0-9]+")))) %>% 
    filter(!(Study %in% c("Colonrectum_L", "Colonrectum_S")), Fasting_C == 2)
  # 1799 fasted subjects left
  
}
ctrl <- prepcontrol()



