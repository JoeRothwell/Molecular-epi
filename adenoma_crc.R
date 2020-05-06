library(readxl)
library(janitor)
library(tidyverse)
dat <- read_xlsx("CRC Metabolomics MasterFile only Biocrates.xlsx", skip = 1, na = c(".", "<LOD")) %>% 
  clean_names() %>% remove_empty("rows")

table(dat$pathology_summary)

# Make pathology groups
dat <- dat %>% mutate(path.group = case_when(
               pathology_summary == 1 ~ "crc",
               pathology_summary %in% 2:3 ~ "adenoma",
               pathology_summary == 4 ~ "polyp",
               pathology_summary %in% 5:6 ~ "normal"
               ))

table(dat$path.group)

# Binary variable for adenoma
dat <- dat %>% mutate(ct = case_when(
               pathology_summary %in% 2:3 ~ 1,
               pathology_summary %in% 5:6 ~ 0
               ))
                      
table(dat$ct)

# Metabolite selection (From Magda sheet) (377 subjects)
# Remove 11 missings (no peak detected)
# Remove those with more than 21% missings
# Recode 998 and 999 with missing, remove missing > 40%

# PCPR2
mat <- dat %>% select(path.group, lyso_pc_a_c16_0:c9) #%>% remove_constant(na.rm = T)
class(mat)
library(Amelia)
missmap(mat, rank.order = F, x.cex = 1)
