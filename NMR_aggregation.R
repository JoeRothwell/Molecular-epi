# Calculation of compound table by aggregation of 107 NMR features corresponding to 43 compounds
# Raw unscaled data received by email 11/10/2019

# 3 compounds have values < 0: formate, hypoxanthine, inosine
# 130 clusters were used to make the 43 compounds
library(readxl)
library(tidyverse)
library(reshape2)

# Get cluster data
dat <- read_xlsx("1510_ClusterAnnotation_E3N.xlsx", sheet = 3, skip = 1)
dat0 <- data.frame(t(dat)) %>% rownames_to_column("cluster")

# Get compound names (with spaces)
cmpds <- read_xlsx("1510_ClusterAnnotation_E3N.xlsx", sheet = 3, col_names = F, n_max = 1) %>% 
  gather() %>% select(value)

# Bind together, drop the sums to leave just the features
feats <- cbind(cmpds, dat0) %>% 
  fill(value) %>%
  filter(!str_detect(cluster, "SUM")) %>%
  select(-cluster)

# Get minimum values
ff <- which(apply(feats[, -1], 1, min) < 0)
feats$value[ff]

# Recalculate sums from groups
cmpds <- feats %>% group_by(value) %>% summarise_all(sum) %>% column_to_rownames("value")
# Get minimum values
which(apply(cmpds, 1, min) < 0)
