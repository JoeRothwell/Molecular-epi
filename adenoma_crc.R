library(readxl)
library(janitor)
library(tidyverse)
library(Amelia)
library(zoo)

# Prep data: replace NAs, clean up names, remove empty rows, rename sarcosine_178 to sarcosine
dat <- read_xlsx("CRC Metabolomics MasterFile only Biocrates.xlsx", skip = 1, na = c(".", "<LOD")) %>% 
  clean_names() %>% remove_empty("rows") %>% mutate(pathsum = pathology_summary) %>%
  rename(sarcosine = sarcosine_178)

# Pathology codings
# 1 = Cancer; 2 = High Grade Dysplasia; 3 = Adenoma (TA, TVA, VA);  4 = Polyp (HP or small TA); 
# 5 = normal (after colonoscopy); 6 = blood donor control
table(dat$pathology_summary)

# Make pathology groups
dat <- dat %>% 
  rename(sex = sex_f_female_m_male, 
         norm.class = normal_classification_for_other_minor_diagnoses_for_colonoscopy_normals) %>%
  mutate(path.group = case_when(
         pathsum == 1 ~ "crc", pathsum %in% 2:3 ~ "adenoma",
         pathsum == 4 ~ "polyp", pathsum %in% 5:6 ~ "normal"
               )) %>%
  mutate(path = if_else(pathsum %in% 1:4, 1, 0))

table(dat$path.group)
#adenoma crc  normal   polyp 
#60      153     103      73

# Remove the normals with inflammatory conditions
dat <- dat %>% filter(pathsum %in% c(1,2,3,4,6) | norm.class == 1)

# Binary variable for adenoma
dat <- dat %>% mutate(ct = case_when(pathsum %in% 2:3 ~ 1, pathsum %in% 5:6 ~ 0))
table(dat$ct)

# Metabolite selection (From Magda sheet) (377 subjects)
# Recode 998 and 999 with missing, remove missing > 40%

# Get Biocrates data only (use glutamate)
# Remove 11 missings (no peak detected)
mat <- dat %>% select(path.group, country:age, lyso_pc_a_c16_0:c9) %>% filter(!is.na(glu))
missmap(mat, rank.order = F, x.cex = 1)

# Convert character columns to numeric, remove metabolites with more than 21% missings (80),
# code normal = 0, pathology = 1
mat <- mat %>% mutate_at(.vars = vars(lyso_pc_a_c16_0:c9), .funs = as.numeric)
mat1 <- mat %>% select_if(~ sum(is.na(.)) < 80) %>% mutate(ct = ifelse(path.group == "normal", 0, 1))
missmap(mat1, rank.order = F, x.cex = 1)
# Add numeric case-control variable

# Impute half min value
mat2 <- na.aggregate(as.matrix(mat1[, -c(1:4, ncol(mat1))]), function(x) min(x)/2)
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



# Get standardised compound names for PLS
cmpd.meta <- read_csv("Biocrates_cmpd_metadata.csv")

# Make data frame of Biocrates compound names
#allnames <- dat %>% select(lyso_pc_a_c16_0:c9)
matnames <- data.frame(cmpd.low = colnames(mat2), ord = 1:length(colnames(mat2)))

# Make new names that match those in adenoma dataset
cmpd.meta2 <- cmpd.meta %>% 
  separate(Compound, into = c("Compound.cl", "Compound.str"),
            remove = F, extra = "merge", fill = "right") %>%
  mutate(cmpd.low = str_to_lower(Compound.str)) %>% 
  mutate(cmpd.low = str_replace(cmpd.low, "lyso", "lyso_")) %>%
  mutate(cmpd.low = str_replace(cmpd.low, "c4_oh_", "c4_oh"))

library(fuzzyjoin)
df1 <- stringdist_right_join(cmpd.meta2, matnames, by = "cmpd.low", max_dist = 0.5)


# Replace matrix names with standardised ones (see match_cmpd_names.R)
mat2a <- mat2
colnames(mat2a) <- df1$Compound

# Get matrices for normal and adenoma, normal and CRC
adenoma <- mat2a[mat1$path.group %in% c("adenoma", "normal"), ] %>% log2 %>% scale
crc     <- mat2a[mat1$path.group %in% c("crc", "normal"), ] %>% log2 %>% scale
polyp   <- mat2a[mat1$path.group %in% c("polyp", "normal"), ] %>% log2 %>% scale

adenoma.meta <- mat1[mat1$path.group %in% c("adenoma", "normal"), ]
# Table sex and country
#   1  2
#F 58 30
#M 55 17


crc.meta <- mat1[mat1$path.group %in% c("crc", "normal"), ]
polyp.meta <- mat1[mat1$path.group %in% c("polyp", "normal"), ]


