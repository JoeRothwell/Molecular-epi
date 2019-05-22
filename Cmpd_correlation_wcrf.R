# Correlations between Biocrates compounds and FAs in small CC
source("CRC_data_prep.R")
source("Metabolic_signatures.R")

FAs.ID <- get.FA.sig(cor.data = T)
Bioc.ID <- get.common.Bioc(cor.data = T)

dat <- inner_join(FAs.ID, Bioc.ID, by = "Idepic") %>% select(-Idepic) 
cormat <- cor(dat, use = "pairwise.complete.obs")
colnames(cormat) <- NULL
library(corrplot)
corrplot(cormat, tl.col = "black", tl.cex = 0.5, order = "hclust", method = "square")

library(reshape2)
all.correlations <- melt(cormat) %>% filter(value != 1) %>% arrange(desc(value))
all.correlations %>% filter(Var1 == "P17_0")
