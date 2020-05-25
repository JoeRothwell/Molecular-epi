library(readxl)
library(janitor)
library(tidyverse)
library(Amelia)
library(zoo)

# Prep data: replace NAs, clean up names, remove empty rows, rename sarcosine_178 to sarcosine
dat <- read_xlsx("CRC Metabolomics MasterFile only Biocrates.xlsx", skip = 1, na = c(".", "<LOD")) %>% 
  clean_names() %>% remove_empty("rows") %>% mutate(pathsum = pathology_summary) %>%
  rename(sarcosine = sarcosine_178)

table(dat$pathology_summary)

# Make pathology groups
dat <- dat %>% rename(sex = sex_f_female_m_male) %>%
  mutate(path.group = case_when(
         pathsum == 1 ~ "crc", pathsum %in% 2:3 ~ "adenoma",
         pathsum == 4 ~ "polyp", pathsum %in% 5:6 ~ "normal"
               ))

table(dat$path.group)

# Binary variable for adenoma
dat <- dat %>% mutate(ct = case_when(pathsum %in% 2:3 ~ 1, pathsum %in% 5:6 ~ 0))
table(dat$ct)

# Metabolite selection (From Magda sheet) (377 subjects)
# Remove 11 missings (no peak detected)
# Remove those with more than 21% missings
# Recode 998 and 999 with missing, remove missing > 40%

# Get Biocrates data only (use glutamate)
#mat <- dat %>% select(path.group, lyso_pc_a_c16_0:c9) %>% filter(!is.na(glu))
mat <- dat %>% select(path.group, country:age, lyso_pc_a_c16_0:c9) %>% filter(!is.na(glu))
missmap(mat, rank.order = F, x.cex = 1)

# Convert character columns to numeric
mat <- mat %>% mutate_at(.vars = vars(lyso_pc_a_c16_0:c9), .funs = as.numeric)
mat1 <- mat %>% select_if(~ sum(is.na(.)) < 80)
missmap(mat1, rank.order = F, x.cex = 1)

# Impute half min value
mat2 <- na.aggregate(as.matrix(mat1[, -c(1:4)]), function(x) min(x)/2)
library(pca3d)
pca <- prcomp(mat2, scale. = T)
dev.off()
pca2d(pca, group = mat1$path.group, legend = "bottomright")
box(which = "plot", lty = "solid")

# mat1 is the unimputed matrix with metadata, mat2 is the imputed matrix without metadata
# PC-PR2
library(pcpr2)
props <- runPCPR2(mat2, mat1[, 1:4])
plot(props, col = "red")

# Get standardised names for matrix (see match_cmpd_names.R)
mat2a <- mat2
colnames(mat2a) <- df1$Compound

# Get matrices for normal and adenoma, normal and CRC
adenoma <- mat2[mat1$path.group %in% c("adenoma", "normal"), ] %>% log2 %>% scale
crc     <- mat2[mat1$path.group %in% c("crc", "normal"), ] %>% log2 %>% scale

# Get compound overlap between adenoma, crc and EPIC CRC datasets. Then place in same order
# Run CRC_prep_data
cmpd_meta <- read_csv("Biocrates_cmpd_metadata.csv")

# Model for transformation to residuals
adj   <- function(x) residuals(lm(x ~ country, data = mat1))



