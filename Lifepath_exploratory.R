# LIFEPATH data exploratory

library(tidyverse)
library(readxl)

# Path to food intake data
path <- "Y:/RepPerso/Fabienne WILM/02_Demandes_ponctuelles/10_LIFEPATH/TABLES"
list.files(path)

# Food intake and other data
meta1 <- read_csv("D01_20171031_LIFEPATH.csv")
meta2 <- read_csv("D01_20161018_LIFEPATH.csv")
meta3 <- read_csv("D01_20150917_LIFEPATH.csv")

# Rawest feature data
raw <- read_delim("C:/J_ROTHWELL/X_AlignedCohorteE3NData_cpmg_ssCitPEG_0612.txt", delim = ";")

# Open other data files
# List of 54 compounds and IDs
dat1 <- read_xlsx("C:/J_ROTHWELL/1505_E3N_Identification.xlsx")

# List of metabolites and corresponding clusters
dat2 <- read_xlsx("C:/J_ROTHWELL/1507_ClusterAnnotation_E3N.xlsx")

# 1623 observations of 44 intensity variables. Looks scaled version of dat 4 and final prepared data
dat3 <- read_tsv("C:/J_ROTHWELL/1507_XMetabolite_std_cpmg_E3N.txt")

# Previous version of file above without scaling
dat4 <- read_tsv("C:/J_ROTHWELL/1507_XMetaboliteE3N_cpmg.txt")

# metadata
dat5 <- read_xlsx("C:/J_ROTHWELL/1603_MatriceY_CohorteE3N_Appar.xlsx", sheet = 6)


# Metadata and intensities
meta <- read_xlsx("C:/J_ROTHWELL/1603_MatriceY_CohorteE3N_Appar.xlsx", sheet = 6)
ints <- read_tsv("C:/J_ROTHWELL/1507_XMetabolite_std_cpmg_E3N.txt")
ints1 <- read_tsv("C:/J_ROTHWELL/1507_XMetaboliteE3N_cpmg.txt", skip = 1)

pca <- prcomp(ints, scale.=F)
pca1 <- prcomp(ints, scale.=F)
plot(pca)
plot(pca1)

library(pca3d)
pca2d(pca)
pca2d(pca1)






