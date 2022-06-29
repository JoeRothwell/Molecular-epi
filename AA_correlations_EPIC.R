# Amino acid correlations in EPIC, manuscript figure 1

source("CRC_prep_data.R")
ctrl <- colon$Cncr_Caco_Clrt == 0

# For 13 P150 AAs only
mat0 <- colon %>% select(contains("Aminoacid_")) %>% select_if(~ sum(is.na(.)) < (nrow(colon)*0.31))

# For 21 P150+180 AAs (removing leucine + isoleucine)
mat0 <- colon %>% select(contains("Aminoacid_"), -Aminoacid_Xleu) 

mat <- colon %>% select(colnames(mat0)) 
mat1 <- mat[ctrl, ]

# Remove Aminoacid from labels. Any character before _
colnames(mat1) <- str_replace(colnames(mat1), ".*_", "")
# Aminoacid before _
colnames(mat1) <- str_replace(colnames(mat1), "Aminoacid_", "")

cormat <- cor(mat1, use = "pairwise.complete.obs")
library(corrplot)
corrplot(cormat, method = "square", tl.col = "black", order = "hclust")
