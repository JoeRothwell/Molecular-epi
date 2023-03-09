library(readxl)
library(janitor)
library(tidyverse)
library(Amelia)
library(zoo)

# Prep data. Get standardised compound names for PLS that match those in adenoma data
cmpd.meta <- read_csv("Biocrates_cmpd_metadata.csv") %>%
  separate(Compound, into = c("Compound.cl", "Compound.str"),
           remove = F, extra = "merge", fill = "right") %>%
  mutate(cmpd.low = str_to_lower(Compound.str)) %>% 
  mutate(cmpd.low = str_replace(cmpd.low, "lyso", "lyso_")) %>%
  mutate(cmpd.low = str_replace(cmpd.low, "c4_oh_", "c4_oh")) 


# Main biocrates data: replace NAs, clean up names, remove empty rows, rename sarcosine_178 to sarcosine
bioc <- read_xlsx("CRC Metabolomics MasterFile only Biocrates.xlsx", skip = 1, na = c(".", "<LOD")) %>% 
  clean_names() %>% remove_empty("rows") %>% mutate(pathsum = pathology_summary) %>%
  rename(sarcosine = sarcosine_178)

# SCFA data to possibly add to analysis: clean and prepare for join to biocrates data
# Note: tidied spreadsheet (manually changed N/A to NA, etc)
scfa <- read_xlsx("Irish Czech_Data_file_2020_JR_DH_mod.xlsx", sheet = 3, skip = 1) %>%
  remove_empty("rows") %>% rename(id = `SCFA Plasma ID`) %>% clean_names()

missmap(scfa, rank.order = F, x.cex = 1)

dat <- bioc %>% left_join(scfa, by = "id") #%>% select(-id, -pathology)

# Pathology codings
# 1 = Cancer; 2 = High Grade Dysplasia; 3 = Adenoma (TA, TVA, VA);  4 = Polyp (HP or small TA); 
# 5 = normal (after colonoscopy); 6 = blood donor control
table(dat$pathology_summary)
table(dat$country, dat$pathology_summary)
# 217 subjects in Ireland, 172 in CR. CR only has colorectal cancer as pathology (125)

# Make pathology groups
dat <- dat %>% 
  rename(sex = sex_f_female_m_male, 
         norm.class = normal_classification_for_other_minor_diagnoses_for_colonoscopy_normals,
         batch = batch_plasma_metabolomics,
         smoke = smoking) %>%
  mutate(path.group = case_when(
         pathsum == 1 ~ "crc", pathsum %in% 2:3 ~ "adenoma",
         pathsum == 4 ~ "polyp", pathsum %in% 5:6 ~ "normal"
               )) ##%>%
  #mutate(path = if_else(pathsum %in% 1:4, 1, 0))

table(dat$path.group)
#adenoma crc  normal   polyp 
#60      153     103      73
table(dat$path.group, dat$country)

# Remove a) normals with inflammatory conditions, b) polyps that are not hyperplastic polyps
#dat <- dat %>% filter(pathsum %in% c(1,2,3,4,6) | norm.class == 1)
#dat <- dat %>% filter(pathsum %in% c(1,2,3,5,6) | histology_of_adenoma %in% 7)

# Add binary variable for adenoma
dat <- dat %>% mutate(ct = case_when(pathsum %in% 2:3 ~ 1, pathsum %in% 5:6 ~ 0))
table(dat$ct)

# Metabolite selection (From Magda sheet) (377 subjects)
# Recode 998 and 999 with missing, remove missing > 40%

# Get main variables only incl. Biocrates and SCFA data.
dat <- dat %>% select(pathsum, path.group, bmi, diabetes, country:age, batch, scfa_analysis_batch,
                      alcohol, alcohol_drinks_week, smoke, lyso_pc_a_c16_0 : c9) %>% 
#Remove samples with no glutamate measurement (11 missings, no peak detected)
               filter(!is.na(glu))

# Or: get main variables only incl. SCFA data.
datsc <- dat %>% select(pathsum, path.group, bmi, diabetes, country:age, batch, scfa_analysis_batch,
                      alcohol, alcohol_drinks_week, smoke, acetic_acid_m_m:valeric_acid_m_m)

missmap(dat, rank.order = F, x.cex = 1)
missmap(datsc, rank.order = F, x.cex = 1)

### Analysis of Biocrates compounds only (for SCFA, go to end of script)
 
# Convert character columns to numeric, remove metabolites with more than 21% missings (80)
#dat1 <- dat %>% mutate_at(.vars = vars(lyso_pc_a_c16_0 : c9), .funs = as.numeric)

# Remove metabolites with more than 21% missings (missings == T OR metadata == T)
# Logical vectors for non-missings (warning, participant data hard coded in cols 1-12)
missvec <- sapply(dat, function(x) sum(is.na(x)) < 80)
metvec <- 1:ncol(dat) %in% 1:12

# code normal = 0, pathology = 1 for logistic regression
dat1 <- dat %>% select_if(missvec == T | metvec == T) %>% mutate(ct = ifelse(pathsum %in% 5:6, 0, 1))

missmap(dat1, rank.order = F, x.cex = 1)

# Impute half min value and subset full matrix
mat <- dat1 %>% select(lyso_pc_a_c16_0 : c9)
mat2 <- na.aggregate(as.matrix(dat1[ , -c(1:10, ncol(dat1))]), function(x) min(x)/2)

# Replace matrix names with standardised ones (see match_cmpd_names.R)
matnames <- data.frame(cmpd.low = colnames(mat2), ord = 1:length(colnames(mat2)))

library(fuzzyjoin)
df1 <- stringdist_right_join(cmpd.meta, matnames, by = "cmpd.low", max_dist = 0.5)
mat2a <- mat2 %>% log2 %>% scale
colnames(mat2a) <- df1$Compound

# Subsets for all, normal vs adenoma, CRC or polyp
adenoma <- mat1$pathsum %in% c(2:3, 5:6)
crc     <- mat1$pathsum %in% c(1, 5:6)
polyp   <- mat1$pathsum %in% c(4, 5:6)

# Table sex and country: F, 58 and 30 in 1 and 2, M 55 and 17 in 1 and 2

# Plot PCA and PCPR2
library(pca3d)
pca <- prcomp(mat2, scale. = T)
#dev.off()
pca2d(pca, group = mat1$path.group, legend = "bottomright")
box(which = "plot", lty = "solid")

# mat1 is the unimputed matrix with metadata, mat2 is the imputed matrix without metadata
# PC-PR2
library(pcpr2)
props <- runPCPR2(mat2, mat1[, c("path.group", "country", "sex", "age", "batch")]) #need to get right columns
plot(props, col = "red")

# Cluster dendrogram of individuals
scalemat <- scale(mat2)
colnames(scalemat) <- NULL
hh <- hclust(dist(scalemat))
plot(hh) # 2 main clusters

# Or with dendextend (from Stack Overflow post)
library(dendextend)
dend <- hclust(dist(scalemat)) %>% as.dendrogram
dend <- dend %>% #set("labels_colors", mat1$country, order_value = TRUE) %>%
  set("labels_colors", as.numeric(as.factor(mat1$path.group)), order_value = TRUE) %>%
  set("labels_cex", .5)
plot(dend, horiz = F)



### SCFA analysis

# Get metadata and SCFA only. Remove 11 missings using NA glu (no peak detected)
mat <- dat %>% select(pathsum, path.group, bmi, diabetes, country:age, batch, alcohol, alcohol_drinks_week, 
                      smoke, lyso_pc_a_c16_0:c9, 
                      acetic_acid_m_m:valeric_acid_m_m) %>%
  filter(if_all(acetic_acid_m_m:valeric_acid_m_m, ~ !is.na(.)))


# Convert character columns to numeric, remove metabolites with more than 21% missings (80),
mat <- mat %>% mutate_at(.vars = vars(acetic_acid_m_m:valeric_acid_m_m), .funs = as.numeric)
mat <- mat %>% mutate_at(.vars = vars(lyso_pc_a_c16_0:valeric_acid_m_m), .funs = as.numeric)

# Get correlations for all compounds or just SCFA
library(corrplot)
cormat <- cor(mat[, -(1:10)], use = "pairwise.complete.obs")
corrplot(cormat, method = "square", tl.col = "black", order = "hclust")
# 2-methylbutyric acid, butyric acid and valeric acid strongly correlated

# Remove metabolites with more than 21% missings (missings == T OR metadata == T)
missvec <- sapply(mat, function(x) sum(is.na(x)) < 80)
metvec <- 1:ncol(mat) %in% 1:10
mat1 <- mat %>% select_if(missvec == T | metvec == T) %>% mutate(ct = ifelse(pathsum %in% 5:6, 0, 1))

#missmap(mat1, rank.order = F, x.cex = 1)

# Impute half min value and subset full matrix
mat2 <- na.aggregate(as.matrix(mat1[, -c(1:10, ncol(mat1))]), function(x) min(x)/2)

