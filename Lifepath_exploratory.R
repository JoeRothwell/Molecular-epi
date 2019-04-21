# LIFEPATH data exploratory

library(tidyverse)
library(readxl)

# Final data files seem to be the following:

# 1623 observations of 44 intensity variables. Looks scaled version of dat 4 and final prepared data
ints <- read_tsv("C:/J_ROTHWELL/1507_XMetabolite_std_cpmg_E3N.txt")
ints <- read_tsv("1507_XMetabolite_std_cpmg_E3N.txt")

# metadata (from XL or csv)
meta <- read_xlsx("C:/J_ROTHWELL/1603_MatriceY_CohorteE3N_Appar.xlsx", sheet = 6)
meta <- read.csv("Lifepath_meta.csv")

# subset for baseline characteristics table
meta0 <- meta %>% 
  select(CT, AGE, BMI, HANCHE, MENOPAUSE, SMK, DIABETE, Life_Alcohol_Pattern_1, BP, Trait_Horm, 
         CO, CENTTIMECat1, FASTING, STOCKTIME, BEHAVIOUR, SUBTYPE, HR, Estro_THM, Pg_seul, 
         SBR, GRADE, STADE, DIAGSAMPLING) %>% 
  mutate_at(vars(-AGE, -BMI, -HANCHE, -DIAGSAMPLING, -STOCKTIME), as.factor)

# ---- Other files

# Rawest feature data
raw <- read_delim("C:/J_ROTHWELL/X_AlignedCohorteE3NData_cpmg_ssCitPEG_0612.txt", delim = ";")

# Excel files with intermediate steps
# List of 54 compounds and IDs
xl1 <- read_xlsx("C:/J_ROTHWELL/1505_E3N_Identification.xlsx")

# Previous version of intensities above without scaling
dat2 <- read_tsv("C:/J_ROTHWELL/1507_XMetaboliteE3N_cpmg.txt")

# List of metabolites and corresponding clusters
xl2 <- read_xlsx("C:/J_ROTHWELL/1507_ClusterAnnotation_E3N.xlsx")

# Metadata and intensities
xl4 <- read_xlsx("C:/J_ROTHWELL/1603_MatriceY_CohorteE3N_Appar.xlsx", sheet = 6)

ints <- read_tsv("C:/J_ROTHWELL/1507_XMetabolite_std_cpmg_E3N.txt")
ints1 <- read_tsv("C:/J_ROTHWELL/1507_XMetaboliteE3N_cpmg.txt", skip = 1)

pca <- prcomp(ints, scale.=F)
pca1 <- prcomp(ints, scale.=F)
plot(pca)
plot(pca1)

library(pca3d)
pca2d(pca)
pca2d(pca1)

# Food intake data ---- 
path <- "Y:/RepPerso/Fabienne WILM/02_Demandes_ponctuelles/10_LIFEPATH/TABLES"
list.files(path)

# Food intake and other data
meta1 <- read_csv("D01_20171031_LIFEPATH.csv")
meta2 <- read_csv("D01_20161018_LIFEPATH.csv")
meta3 <- read_csv("D01_20150917_LIFEPATH.csv")






