# PCPR2 on the BC metabolomics study

library(tidyverse)

# 1623 observations of 44 intensity variables. Looks scaled version of dat 4 and final prepared data
Xdata <- read_tsv("1507_XMetabolite_std_cpmg_E3N.txt") %>% filter(CODBMB != 11094738)
meta <- read_csv("Lifepath_meta.csv") %>% filter(CODBMB != 11094738)
alldata <- inner_join(meta, Xdata, by = "CODBMB")
# 1109473 is not in the metadata anyway

X_DataMatrixScaled <- dplyr::select(alldata, `3Hydroxybutyrate`:Succinate) %>% as.matrix
Z_Meta <- alldata %>%
  select(CT, WEEKS, PLACE, AGE, BMI, MENOPAUSE, FASTING, SMK, DIABETE, CENTTIMECat1, SAMPYEAR, STOCKTIME) %>%
  mutate_at(vars(-AGE, -BMI, -WEEKS, -STOCKTIME), as.factor)

library(pcpr2)
props <- runPCPR2(X_DataMatrixScaled, Z_Meta)
plotProp(props)
title("Sources of variation in BC cases and controls", font.main = 1)