# Test to find differences between received datasets 

library(tidyverse)
library(readxl)

# Original data, 1623 observations of 44 compounds - appears to be Pareto scaled
ints <- read_tsv("1507_XMetabolite_std_cpmg_E3N.txt")

# Update 11/10/19: Updated with unscaled data (will scale to unit variance). Added "_unscaled" to filename.
# Note: NAC1 and NAC2 are merged and only 1582 observations
ints0 <- read_tsv("1510_XMetaboliteE3N_cpmg_unscaled.txt") #%>% as.matrix
ints1 <- scale(ints0)

# Lifestyle data. Subset variables needed
meta <- read_csv("Lifepath_meta.csv", na = "9999")
samples <- ints$CODBMB %in% meta$CODBMB
ints <- ints[samples, -1]

# Exploratory of different datasets
library(MetabolAnalyze)

# 1. Original data (1507); 2. New data (no scaling); 3. New data (unit scaled); New data (Pareto scaled)
new     <- ints0 %>% as.matrix
unit    <- scale(ints0)
pareto1 <- scaling(ints0, type = "pareto")

ll <- lapply(list(ints, new, unit, pareto1), function(x) prcomp(x, scale. = F, center = F))
par(mfrow = c(2,2))
lapply(ll, function(x,y) { pca2d(x)
  title(main = paste("Scores plot"))
  box(which = "plot")
})

#how to make hotpink, limegreen, orange, dodgerblue and give different titles?

# Conclusion: the old dataset seems to be a unit variance scaled version of the new dataset
# (although they are not exactly the same; the new datasets has merged NAC1 and NAC2)
# Decision: use unscaled dataset and apply unit scaling, discard old scaled dataset
