# Supplemental table of compounds and CVs
# Biocrates from IARC and Helmholtz, FAs from IARC
library(readxl)
library(janitor)
library(tidyverse)

bioc <- read_csv("Biocrates_cmpd_metadata.csv") 
facids <- read_csv("FA_compound_data.csv")

# Helmholtz Biocrates CVs. Read files, put together and calculate intrabatch mean
cvs <- list.files("helmoltz_CVs") %>% 
  map_df( ~ read_delim(paste("helmoltz_CVs", ., sep = "/"), delim = ";")) %>% 
  select(C3:H1) %>% summarise_all(funs(mean))

# Get list of compound names from full data
full <- read_delim("helmholtz_CV_data/2018-01-16_Conc_P160262_P-17044-BBMRI-P1_Kit-LOD.csv", 
                 skip = 1, delim = ";") %>% select(-contains("Status")) %>% select(C3:H1)

# Check names are the same and replace
df <- data.frame(n1 = colnames(cvs), n2 = colnames(full))
colnames(cvs) <- colnames(full)
helm <- gather(cvs, key = "shortname") %>% mutate(CV_Helm = round(value*100, 1))

# IARC Biocrates CVs. 2nd sheet for serum
cvbioc <- read_xlsx("Epic CRC_QC Results.xlsx", skip = 1, sheet = 2) %>% 
  filter(str_detect(Compounds, "CV% Interbatch")) %>% slice(2) %>% select(-(1:4))
iarc <- gather(cvbioc, key = "shortname") %>% slice(-(1:4)) %>% mutate(CV_Iarc = round(value, 1))

# Get table of coefficients for signatures
load("coefficient_tables.Rdata")

# Note: the following are excluded from signature because of absence from CRC
# 159 compounds reduced to 155
#"Biogenic_Nitro_Tyr" "Biogenic_Ac_Orn"    "Biogenic_Met_So"    "Biogenic_Total_Dma"

# Join two sets of CVs to Biocrates metadata
bioc1 <- bioc %>% left_join(iarc, by = "shortname")
bioc2 <- bioc1 %>% left_join(helm, by = "shortname")

pltdata1 <- pltdata %>% select(compound:Coefficient)

# Join coefficients table to biocrates metadata
bioc3 <- pltdata1 %>% left_join(bioc2, by = c("compound" = "displayname")) %>%
  select(class, compound, Coefficient, CV_Iarc, CV_Helm) %>% arrange(class)


# Fatty acids
cvfa <- read_xlsx("FA-QC-CV.xlsx", col_types = c("text", "numeric")) %>% 
  rename(displayname = "Fatty Acid") %>% mutate(CV = round(`CV%`, 1))

facidCV1 <- faplot %>% left_join(cvfa, by = c("compound" = "displayname")) %>%
  select(class, compound, Coefficient, CV_Iarc = CV) %>% arrange(class)

# Bind two tables together
allCVdat <- bind_rows(Endogenous.metabolites = bioc3, Fatty.acids = facidCV1, .id = "Platform")


# Get CVs for first two ACs (not supplied with data)
fullcvs <- list.files("helmholtz_CV_data") %>% 
  map_df( ~ read_delim(paste("helmholtz_CV_data", ., sep = "/"), delim = ";", skip = 1)) %>% 
  filter(`Sample Identification` == "Ref_Plasma-Hum_PK3") %>%
  select(KitBarcodeNr, C0, C2, C3) %>% mutate_at(vars(C0, C2, C3), funs(as.numeric))

# Overall CV (not used in table)
overallcvs <- fullcvs %>% summarise_all(~ sd(.)/mean(.))

#Batch average CVs
batchcvs <- fullcvs %>% group_by(KitBarcodeNr) %>% summarise_all(~ sd(.)/mean(.))
batchcvs %>% summarise_all(mean)
#C0: 7.0%, C2: 7.1%
