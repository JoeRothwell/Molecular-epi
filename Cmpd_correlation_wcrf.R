# Correlations between Biocrates compounds and FAs in small CC
source("CRC_prep_data.R")
source("CRC_get_signatures.R")

# Get compounds and idepics from CRC1 and join together by Idepic
FAs.ID <- get.plsdata(CRCfa1, cor.data = T)
Bioc.ID <- select.ctrl.cmpds(crc1, cor.data = T)

dat <- inner_join(FAs.ID, Bioc.ID, by = "Idepic") %>% select(-Idepic) 

# Get compound metadata for short names
meta.bioc <- read_csv("Biocrates_cmpd_metadata.csv") %>% select(Compound, displayname)
meta.fa   <- read_csv("FA_compound_data.csv") %>% select(Compound, displayname)
meta <- bind_rows(meta.bioc, meta.fa)

mm <- data.frame(Compound = colnames(dat))

all <- left_join(mm, meta, by = "Compound")

# Calculate correlations
cormat <- cor(dat, use = "pairwise.complete.obs")
rownames(cormat) <- all$displayname
colnames(cormat) <- all$displayname

# Correlation matrix (large)
library(corrplot)
corrplot(cormat, tl.col = "black", tl.cex = 0.5, order = "hclust", method = "square")

# Dendrogram. First get order of clustered compounds
dend <- cormat %>% dist %>% hclust %>% as.dendrogram

# Make ordered numeric vector for label colour
ord <- order.dendrogram(dend)
cmpds <- c(rep(2, 34), rep(1, 180-34))[ord]

library(dendextend)
dend1 <- cormat %>% dist %>% hclust %>% as.dendrogram %>%
  set("branches_k_color", k=4) %>% 
  set("labels_col", cmpds) %>%
  set("labels_cex", c(0.6)) #%>%
#hang.dendrogram(hang_height = 0.3)

circlize_dendrogram(dend1, dend_track_height = 0.85)

library(circlize)
circlize

# Get table of correlations in descending order
library(reshape2)
all.correlations <- melt(cormat) %>% filter(value != 1) %>% arrange(desc(value))
all.correlations %>% filter(Var1 == "P17_0")
t1 <- all.correlations %>% filter(Var1 %in% meta.fa$displayname & Var2 %in% meta.bioc$displayname)

# Get correlations > 0.5 in correct order
t2 <- melt(cormat) %>% 
  filter(value > 0.5, Var1 %in% meta.fa$displayname & Var2 %in% meta.bioc$displayname)
colnames(t2) <- c("displayname2", "displayname", "value")


# Heatmap Biocrates vs fatty acids only
mat <- acast(t1, Var2 ~ Var1, value.var = "value")

library(pheatmap)
pheatmap(mat, fontsize = 8, color = colorRampPalette(rev(brewer.pal(n = 7, name ="RdYlBu")))(100),)

# Make data for annotations
cmpd_meta <- read_csv("Biocrates_cmpd_metadata.csv")
df <- data_frame(displayname = rownames(mat))
annotation_df <- inner_join(df, cmpd_meta, by = "displayname")

t3 <- inner_join(t2, annotation_df, by = "displayname")

row_ha <- rowAnnotation(Class = annotation_df$class)

# With Complex Heatmap
# https://jokergoo.github.io/ComplexHeatmap-reference/book/a-single-heatmap.html
library(ComplexHeatmap)
library(circlize)
library(RColorBrewer)
Heatmap(mat, col = colorRampPalette(rev(brewer.pal(n = 7, name ="RdYlBu")))(100),
        #column_title = "Fatty acids", 
        cluster_rows = F,
        #row_title = "Endogenous metabolites", 
        #column_km = 3, 
        #row_km = 3, 
        row_split = annotation_df$class,
        row_title = NULL,
        row_names_side = "right", show_row_names = T,
        row_names_gp = gpar(fontsize = 10),
        column_names_gp = gpar(fontsize = 10),
        left_annotation = row_ha,
        #row_gap = unit(0, "mm"),
        border = F)

ha <- rowAnnotation(foo = anno_mark(at = c(10, 20, 50, 100), labels = c("metabo1", "metabo2", "metabo3", "metabo4")))
Heatmap(mat, right_annotation = ha)














