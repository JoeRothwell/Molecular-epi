library(readxl)
library(janitor)
dat <- read_xlsx("CRC Metabolomics MasterFile only Biocrates.xlsx", skip = 1) %>% 
  clean_names() %>% remove_empty("rows")

table(dat$pathology_summary)
