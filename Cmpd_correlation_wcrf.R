# Correlations between Biocrates compounds and FAs in small CC
source("CRC_prep_data.R")
source("CRC_get_signatures.R")

# Whole case-control 1
# Get compounds and idepics from CRC1 and join together by Idepic (controls only)
FAs.ID <- crc1fa[, c("Idepic", intersect(colnames(concs), colnames(crcfa)))]
Bioc.ID <- crc1[, c("Idepic", intersect(colnames(crc1p), colnames(ctrls)))] %>% distinct()
dat <- inner_join(FAs.ID, Bioc.ID, by = "Idepic") %>% select(-Idepic) 

# Or

# Controls from CC1 and controls dataset
# Get compounds and idepics from CRC1 and join together by Idepic (controls only)
#FAs.ID <- CRCfa.ctrl %>% select(Idepic, one_of(common.cols))
#FAs.IDa <- CRCfa.ctrl[, c("Idepic", intersect(colnames(concs), colnames(crcfa)))]
FAs.IDb <- crc1fa[crc1fa$Cncr_Caco_Clrt == 0, c("Idepic", intersect(colnames(concs), colnames(crcfa)))]

#Bioc.ID <- select.ctrl.cmpds(list(ctrl, crc1), cor.data = T) %>% distinct() # remove dupes
Bioc.ID <- crc1[crc1$Cncr_Caco_Clrt == 0, c("Idepic", intersect(colnames(crc1p), colnames(ctrls)))]
dat <- inner_join(FAs.ID, Bioc.ID, by = "Idepic") %>% select(-Idepic) 

# Get compounds and idepics from controls datasets and join by Idepic
FAs.ID2 <- fa.ctrl %>% select(Idepic, one_of(common.cols))
Bioc.ID2 <-  ctrl %>% select(Idepic, one_of(sort(colnames(ctrlA))))
dat2 <- inner_join(FAs.ID2, Bioc.ID2, by = "Idepic") %>% select(-Idepic)

allcordat <- rbind(dat, dat2)
dim(allcor)
# Total of 438 + 241 observations for correlation
# (Did not give very results, perhaps from combining different labs)


#----

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

# Get table of correlations in descending order
library(reshape2)
cordf <- melt(cormat) %>% filter(value != 1) %>% arrange(desc(value))

# Keep Biocrates - fatty acids correlations only
t1 <- cordf %>% filter(Var1 %in% meta.fa$displayname & Var2 %in% meta.bioc$displayname)
# or
# Keep Biocrates compounds with max correlation < 0.25 only (remove some low correlations from heatmap)
t1 <- t1 %>% group_by(Var2) %>% filter(abs(max(value)) > 0.3)

#----

# Get correlations > 0.5 in correct order (for supplemental table)
t2 <- t1 %>% filter(value > 0.55)
colnames(t2) <- c("displayname2", "displayname", "value")

# Heatmap Biocrates vs fatty acids only
mat <- acast(t1, Var2 ~ Var1, value.var = "value")

# Annotated heatmap of correlations
# https://jokergoo.github.io/ComplexHeatmap-reference/book/a-single-heatmap.html

# Make data for annotations
cmpd_meta <- read_csv("Biocrates_cmpd_metadata.csv")
df <- tibble(displayname = rownames(mat))
anno_df <- inner_join(df, cmpd_meta, by = "displayname")

t3 <- inner_join(t2, anno_df, by = "displayname")
corr.cmpds <- unique(t3$displayname)
cmpd.pos <- which(anno_df$displayname %in% corr.cmpds)

library(ComplexHeatmap)
library(circlize)
library(RColorBrewer)

# Make object and legend for class annotations (remove AAs, BCs and monosaccharides if filtering)
row_ha <- rowAnnotation(Class = anno_df$class, annotation_legend_param = list(
                        labels = c("Acylcarnitines", "Amino acids", "Biogenic amines", 
                                   "LysoPC", "Monosacchaaride", 
                                   "PC (acyl-acyl)", "PC (acyl-alkyl)", "Sphingolipids")))

# Make object for annotated compounds on right
ha <- rowAnnotation(foo = anno_mark(at = cmpd.pos, labels = corr.cmpds, 
                                    labels_gp = gpar(fontsize = 9)))

# Colour ramp palette is taken from pheatmap
rownames(mat) <- NULL
Heatmap(mat, col = colorRampPalette(rev(brewer.pal(n = 7, name ="RdYlBu")))(100),
        #column_title = "Fatty acids", row_title = "Endogenous metabolites", 
        #cluster_rows = T,
        #cluster_columns = T,
        #row_split = anno_df$class,
        row_title = NULL,
        row_names_side = "right", show_row_names = T,
        row_names_gp = gpar(fontsize = 10),
        column_names_gp = gpar(fontsize = 10),
        left_annotation = row_ha,
        right_annotation = ha,
        heatmap_legend_param = list(title = "Correlation"),
        row_gap = unit(0, "mm"),
        border = F)


# Correlation matrix (large) (supplemental data)
library(corrplot)
corrplot(cormat, tl.col = "black", tl.cex = 0.3, order = "hclust", method = "color")



# Not in manuscript

#library(pheatmap)
#pheatmap(mat, fontsize = 8, color = colorRampPalette(rev(brewer.pal(n = 7, name ="RdYlBu")))(100),)

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

library(circlize)
circlize_dendrogram(dend1, dend_track_height = 0.85)














