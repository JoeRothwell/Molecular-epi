# Exploratory analysis of NMR compound datasets

library(tidyverse)
library(readxl)

# Original data, 1623 observations of 44 compounds - appears to be Pareto scaled
ints <- read_tsv("1507_XMetabolite_std_cpmg_E3N.txt")

# Update 11/10/19: Updated with unscaled data (will scale to unit variance). Added "_unscaled" to filename.
# Note: NAC1 and NAC2 are merged and only 1582 observations
ints0 <- read_tsv("1510_XMetaboliteE3N_cpmg_unscaled.txt") #%>% as.matrix
ints1 <- scale(ints0)

# Lifestyle data. Subset variables needed
meta <- read_csv("Lifepath_meta.csv", na = "9999") %>% select(CODBMB, CT, MATCH, PLACE, AGE, BMI, 
        BP, RTH, ALCOHOL, MENOPAUSE, FASTING, SMK, DIABETE, CENTTIMECat1, CENTTIME, SAMPYEAR, STOCKTIME, DURTHSDIAG) 

# Exploratory of different datasets

# CLR models to get odds ratios for metabolites
dat <- left_join(meta, ints, by = "CODBMB")
dat0 <- cbind(meta, ints0)

library(MetabolAnalyze)

# 1. Original data (1507)
# 2. New data (no scaling)
# 3. New data (unit scaled)
# 4. New data (Pareto scaled)

first <- dat %>% select(-(CODBMB:DURTHSDIAG)) %>% as.matrix
new <- ints0 %>% as.matrix
unit    <- scale(ints0)
pareto1 <- scaling(ints0, type = "pareto")

pca.first <- prcomp(first, scale. = F, center = F)
pca.new <- prcomp(new, scale. = F, center = F)
pca.unit <- prcomp(unit, scale. = F, center = F)
pca.pareto1 <- prcomp(pareto1, scale. = F, center = F)

library(pca3d)
par(mfrow = c(2,2))
pca2d(pca.first, col = "hotpink")
box(which = "plot")
title(main = "1507_XMetabolite_std_cpmg_E3N.txt")

pca2d(pca.new, col = "limegreen")
box(which = "plot")
title(main = "1510_XMetaboliteE3N_cpmg.txt")

pca2d(pca.unit, col = "orange")
box(which = "plot")
title(main = "1510_XMetaboliteE3N_cpmg.txt, unit")

pca2d(pca.pareto1, col = "dodgerblue")
box(which = "plot")
title(main = "1510_XMetaboliteE3N_cpmg.txt, Pareto")

# Conclusion: the old dataset seems to be a unit variance scaled version of the new dataset
# (although they are not exactly the same)
