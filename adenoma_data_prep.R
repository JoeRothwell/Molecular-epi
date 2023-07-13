library(readxl)
library(janitor)
library(tidyverse)
library(Amelia)
library(zoo)

# Put datasets together: compound metadata, main dataset and SCFA
# Get standardised compound names for PLS that match those in adenoma data
cmpd.meta <- read_csv("Biocrates_cmpd_metadata.csv") %>%
  separate(Compound, into = c("Cmpd.cl", "Cmpd.str"), remove = F, extra = "merge", fill = "right") %>%
  mutate(cmpd.low = str_to_lower(Cmpd.str)) %>% 
  mutate(cmpd.low = str_replace(cmpd.low, "lyso", "lyso_")) %>%
  mutate(cmpd.low = str_replace(cmpd.low, "c4_oh_", "c4_oh")) 


# Main biocrates data: replace NAs, clean up names, remove empty rows, rename sarcosine_178 to sarcosine
alldata <- read_xlsx("CRC Metabolomics MasterFile only Biocrates.xlsx", skip = 1, na = c(".", "<LOD")) %>% 
  clean_names() %>% remove_empty("rows") %>% mutate(pathsum = pathology_summary) %>%
  rename(sarcosine = sarcosine_178)

# SCFA data to possibly add to analysis: clean and prepare for join to biocrates data
# Note: tidied spreadsheet (manually changed N/A to NA, etc)
scfa <- read_xlsx("Irish Czech_Data_file_2020_JR_DH_mod.xlsx", sheet = 3, skip = 1) %>%
  remove_empty("rows") %>% rename(id = `SCFA Plasma ID`) %>% clean_names()

missmap(scfa, rank.order = F, x.cex = 1)

# Add the SCFA data
dat0 <- alldata %>% left_join(scfa, by = "id") #%>% select(-id, -pathology)


### Pathology codings
table(dat0$pathology_summary)

# 1 = Cancer; to be coded as "CRC"
# 2 = High grade dysplasia (HGD); 3 = adenoma (TA, TVA, VA); to be coded as "adenoma"
# 4 = Polyp (hyperplastic polyp HP or small tubulovillous adenoma); to be coded as "polyp"
# 5 = Normal after colonoscopy, 6 = blood donor control; to be coded as "normal"

table(dat0$country, dat0$pathology_summary)
# 217 subjects in Ireland, 172 in CR. CR only has colorectal cancer as pathology (125)

# Rename and add categorical pathology group variable
dat <- dat0 %>% 
  rename(sex = sex_f_female_m_male, 
         norm.class = normal_classification_for_other_minor_diagnoses_for_colonoscopy_normals,
         batch = batch_plasma_metabolomics,
         smoke = smoking) %>%
  mutate(path.group = case_when(
         pathsum == 1 ~ "crc", pathsum %in% 2:3 ~ "adenoma",
         pathsum == 4 ~ "polyp", pathsum %in% 5:6 ~ "normal"
               )) ##%>%
  #mutate(path = if_else(pathsum %in% 1:4, 1, 0))

# Specify cols to convert to factor
varlist <- c("country", "sex", "batch", "diabetes", "smoke")

table(dat$path.group)
#adenoma crc  normal   polyp 
#60      153     103      73
table(dat$path.group, dat$country) # only CRC and normal for country 2 (CR)


# Sensitivity analyses:
# Remove a) normals with inflammatory conditions, b) polyps that are not hyperplastic polyps
dat.reset <- dat
dat <- dat %>% filter(pathsum %in% c(1:6) | norm.class == 1)
dat <- dat %>% filter(pathsum %in% c(1:6) | histology_of_adenoma %in% 7)


# Convert to factor and add binary variable for adenoma
dat <- dat %>% mutate(across(all_of(varlist), as.factor)) #%>% 
  #mutate(ct = case_when(pathsum %in% 2:3 ~ 1, pathsum %in% 5:6 ~ 0))

#table(dat$ct)

# Metabolite selection (From Magda sheet) (377 subjects)
# Recode 998 and 999 with missing, remove missing > 40%

# Get main variables only incl. Biocrates data.
# Remove samples with no biocrates measurements using glutamate (11 missings, no peak detected)
dat <- dat %>% select(pathsum, path.group, bmi, diabetes, country:age, batch, 
                      scfa_analysis_batch, alcohol, alcohol_drinks_week, smoke, lyso_pc_a_c16_0 : c9) %>% 
                      filter(!is.na(glu))

colnames(dat)

# Or: get main variables only incl. SCFA data.
datsc <- dat %>% select(pathsum, path.group, bmi, diabetes, country:age, batch, scfa_analysis_batch,
                      alcohol, alcohol_drinks_week, smoke, acetic_acid_m_m:valeric_acid_m_m)

missmap(dat, rank.order = F, x.cex = 1)
missmap(datsc, rank.order = F, x.cex = 1)

### Analysis of Biocrates compounds only (for SCFA, go to end of script)
 
# Convert character columns to numeric, remove metabolites with more than 21% missings (80)
dat1 <- dat %>% mutate(across(lyso_pc_a_c16_0 : c9, as.numeric))

# Remove metabolites with more than 21% missings (missings == T OR metadata == T)
# Logical vectors for non-missings (warning, participant data hardcoded in cols 1-12)
missvec <- sapply(dat1, function(x) sum(is.na(x)) < 80)

# Logical vector of metabolite/non-metabolite columns
metab <- 1:ncol(dat1) %in% 1:12

# Keep column if participant data or metabolite with less than 21% missing
# code normal = 0, pathology = 1 for logistic regression
dat1 <- dat1 %>% select_if(missvec == T | metab == T) %>% mutate(ct = ifelse(pathsum %in% 5:6, 0, 1))

missmap(dat1, rank.order = F, x.cex = 1)

# Impute half min value and subset full matrix
mat <- dat1 %>% select(lyso_pc_a_c16_0 : c5)
# mat <- dat1 %>% select(lyso_pc_a_c16_0 : c9)
mat2 <- na.aggregate(as.matrix(dat1[ , -c(1:12, ncol(dat1))]), function(x) min(x)/2)
mat2a <- mat2 %>% scale

# Replace matrix names with standardised ones (see match_cmpd_names.R)
matnames <- data.frame(cmpd.low = colnames(mat2), ord = 1:length(colnames(mat2)))
library(fuzzyjoin)
df1 <- stringdist_right_join(cmpd.meta, matnames, by = "cmpd.low", max_dist = 0.5)
colnames(mat2a) <- df1$Compound

# Subsets for all, normal vs adenoma, CRC or polyp
adenoma <- dat1$pathsum %in% c(2:3, 5:6)
crc     <- dat1$pathsum %in% c(1, 5:6)
crc.cr  <- dat1$pathsum %in% c(1, 5:6) & dat1$country == 2
polyp   <- dat1$pathsum %in% c(4, 5:6)

# Table sex and country: F, 58 and 30 in 1 and 2, M 55 and 17 in 1 and 2

# Perform and plot PCA
pca <- prcomp(mat2, scale. = T)
dev.off()

scores <- pca$x[, 1:2]
plot(scores[ , 1], scores[ , 2], col = as.factor(dat1$path.group))

# pc3d package was removed from CRAN
#library(pca3d)
#pca2d(pca, group = dat1$path.group, legend = "bottomright")
#box(which = "plot", lty = "solid")

# CRCs are at visibly higher values along PC1


# dat1 is the unimputed matrix with metadata, mat2 is the imputed matrix without metadata
# PC-PR2
library(pcpr2)
props <- runPCPR2(mat2, dat1[, c("path.group", "country", "sex", "age", "batch")]) #need to get right columns
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
                      smoke, lyso_pc_a_c16_0 : c9, 
                      acetic_acid_m_m : valeric_acid_m_m) %>%
  filter(if_all(acetic_acid_m_m : valeric_acid_m_m, ~ !is.na(.)))


# Convert character columns to numeric, remove metabolites with more than 21% missings (80),
mat <- mat %>% mutate_at(.vars = vars(acetic_acid_m_m : valeric_acid_m_m), .funs = as.numeric)
mat <- mat %>% mutate_at(.vars = vars(lyso_pc_a_c16_0 : valeric_acid_m_m), .funs = as.numeric)

# Get correlations for all compounds or just SCFA
library(corrplot)
cormat <- cor(mat[, -(1:10)], use = "pairwise.complete.obs")
corrplot(cormat, method = "square", tl.col = "black", order = "hclust")
# 2-methylbutyric acid, butyric acid and valeric acid strongly correlated

# Remove metabolites with more than 21% missings (missings == T OR metadata == T)
missvec <- sapply(mat, function(x) sum(is.na(x)) < 80)
metab <- 1:ncol(mat) %in% 1:10
mat1 <- mat %>% select_if(missvec == T | metab == T) %>% mutate(ct = ifelse(pathsum %in% 5:6, 0, 1))

#missmap(mat1, rank.order = F, x.cex = 1)

# Impute half min value and subset full matrix
mat2 <- na.aggregate(as.matrix(mat1[, -c(1:10, ncol(mat1))]), function(x) min(x)/2)



# N for adenoma/CRc models

dim(dat1)
adenoma.df <- dat1[adenoma, ]
crc.df <- dat1[crc, ]
crc.cr.df <- dat1[crc.cr, ]
polyp.df <- dat1[polyp, ]

table(path = adenoma.df$path.group)
table(country = adenoma.df$country, path = adenoma.df$path.group)

table(path = crc.df$path.group)
table(country = crc.df$country, path = crc.df$path.group)

table(path = crc.cr.df$path.group)
table(country = crc.cr.df$country, path = crc.cr.df$path.group)

table(path = polyp.df$path.group)
table(country = polyp.df$country, path = polyp.df$path.group)



