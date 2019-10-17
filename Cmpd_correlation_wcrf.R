# Correlations between Biocrates compounds and FAs in small CC
source("CRC_data_prep.R")
source("Metabolic_signatures.R")

FAs.ID <- get.plsdata(CRCfa1, cor.data = T)
Bioc.ID <- select.ctrl.cmpds(crc1, cor.data = T)

dat <- inner_join(FAs.ID, Bioc.ID, by = "Idepic") %>% select(-Idepic) 

meta.bioc <- read.csv("Biocrates_cmpd_metadata.csv") %>% select(Compound, displayname)
meta.fa   <- read.csv("FA_compound_data.csv") %>% select(Compound, displayname)
meta <- bind_rows(meta.bioc, meta.fa)

mm <- data.frame(Compound = colnames(dat))

all <- left_join(mm, meta, by = "Compound")

cormat <- cor(dat, use = "pairwise.complete.obs")
rownames(cormat) <- all$displayname

library(corrplot)
corrplot(cormat, tl.col = "black", tl.cex = 0.5, order = "hclust", method = "square")

#clustmat <- hclust(dist(cormat))
#dend1 <- as.dendrogram(clustmat)

library(dendextend)
dend1 <- cormat %>% dist %>% hclust %>% as.dendrogram %>%
  set("branches_k_color", k=4) %>% 
  #set("labels_col", value = c(1,2,3,4), k=4) %>%
  set("labels_cex", c(0.6)) #%>%
#hang.dendrogram(hang_height = 0.3)

circlize_dendrogram(dend1, dend_track_height = 0.85)

library(circlize)
circlize

library(reshape2)
all.correlations <- melt(cormat) %>% filter(value != 1) %>% arrange(desc(value))
all.correlations %>% filter(Var1 == "P17_0")
