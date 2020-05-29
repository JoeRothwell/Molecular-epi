# Correlations between Biocrates compounds and FAs in small CC
load("studyA_data_cor.Rdata")
library(tidyverse)

# Get compound metadata for short names
meta.bioc <- read_csv("Biocrates_cmpd_metadata.csv")
cmpd.meta <- meta.bioc %>% select(Compound, displayname)
meta.fa   <- read_csv("FA_compound_data.csv") %>% select(Compound, displayname)
meta <- bind_rows(meta.bioc, meta.fa)

# Get compounds and idepics from CRC1 and join together by Idepic (controls only)
# Fatty acids, just case-control 1 or just controls
FAs.ID <- crc1fa[, c("Idepic", intersect(colnames(concs), colnames(crcfa)))]
FAs.ID <- crc1fa[crc1fa$Cncr_Caco_Clrt == 0, c("Idepic", intersect(colnames(concs), colnames(crcfa)))]

# Biocrates, just case-control 1 or just controls
Bioc.ID <- crc1[, c("Idepic", intersect(colnames(crc1p), colnames(ctrls)))] %>% distinct()
Bioc.ID <- crc1[crc1$Cncr_Caco_Clrt == 0, c("Idepic", intersect(colnames(crc1p), colnames(ctrls)))]

# Get compounds and idepics from CRC1 and join together by Idepic
join.concs <- inner_join(FAs.ID, Bioc.ID, by = "Idepic") %>% select(-Idepic) 

#----

# Add 241 obs from pooled controls dataset
# Get compounds and idepics from controls datasets and join by Idepic
FAs.ID2 <- fa.ctrl[, c("Idepic", intersect(colnames(concs), colnames(crcfa))) ]
Bioc.ID2 <-  ctrl[ , c("Idepic", intersect(colnames(crc1p), colnames(ctrls)))]
join.concs1 <- inner_join(FAs.ID2, Bioc.ID2, by = "Idepic") %>% select(-Idepic)
join.concs2 <- rbind(join.concs, join.concs2)
dim(allcordat)
# Total of 439 + 241 = 680 observations for correlation
# Changed results strangely: lots of correlations for 15:1


#----

names.df <- data.frame(Compound = colnames(join.concs))
all <- left_join(names.df, meta, by = "Compound")

# Calculate correlations and set names to display names
cormat <- cor(join.concs, use = "pairwise.complete.obs")
rownames(cormat) <- all$displayname
colnames(cormat) <- all$displayname

# Melt cormat to get table of correlations in descending order
#library(reshape2)
#cordf <- melt(cormat) %>% filter(value != 1) %>% arrange(desc(value))

library(corrr)
cordf <- cormat %>% as_cordf() %>% stretch(na.rm = T)

# Keep Biocrates - fatty acids correlations only
t1 <- cordf %>% filter(x %in% meta.fa$displayname & y %in% meta.bioc$displayname)
# or
# Keep Biocrates compounds with max correlation < 0.25 only (remove some low correlations from heatmap)
t1 <- t1 %>% group_by(y) %>% filter(abs(max(r)) > 0.3)

#----

# Supplemental table of correlations: get main correlations in descending order
t2 <- t1 %>% filter(r > 0.55)
colnames(t2) <- c("displayname2", "displayname", "value")

# Heatmap Biocrates vs fatty acids only
mat <- acast(t1, y ~ x, value.var = "r")

# Annotated heatmap of correlations
# https://jokergoo.github.io/ComplexHeatmap-reference/book/a-single-heatmap.html

# Make data for annotations
df <- tibble(displayname = rownames(mat))
anno.df <- inner_join(df, meta.bioc, by = "displayname")

t3 <- inner_join(t2, anno.df, by = "displayname")
corr.cmpds <- unique(t3$displayname)
cmpd.pos <- which(anno.df$displayname %in% corr.cmpds)

library(ComplexHeatmap)
library(circlize)
library(RColorBrewer)

# Make object and legend for class annotations (remove AAs, BCs or monosaccharides if filtering)
row_ha <- rowAnnotation(Class = anno.df$class, annotation_legend_param = list(
                        labels = c("Acylcarnitines", "Amino acids", #"Biogenic amines", 
                                   "LysoPC", #"Monosacchaaride", 
                                   "PC (diacyl)", "PC (acyl-alkyl)", "Sphingolipids")))

# Make object for annotated compounds on right
ha <- rowAnnotation(foo = anno_mark(at = cmpd.pos, labels = corr.cmpds, 
                                    labels_gp = gpar(fontsize = 9)))
#ha2 <- rowAnnotation(foo = anno_mark(at = cmpd.pos, labels = corr.cmpds, 
 #                                   labels_gp = gpar(fontsize = 9)))

# Colour ramp palette is taken from pheatmap
rownames(mat) <- NULL
Heatmap(mat, col = colorRampPalette(rev(brewer.pal(n = 7, name ="RdYlBu")))(100),
        column_title = "Fatty acids", row_title = "Endogenous metabolites", 
        row_title_gp = gpar(fontsize = 11), column_title_gp = gpar(fontsize = 11),
        row_title_side = "right", column_title_side = "bottom",
        #cluster_rows = T,
        #row_title = NULL,
        #row_names_side = "right", 
        #show_row_names = F, # affects the spacing of annotations
        row_names_gp = gpar(fontsize = 10),
        column_names_gp = gpar(fontsize = 10),
        left_annotation = row_ha,
        right_annotation = ha,
        #bottom_annotation = ha,
        heatmap_legend_param = list(title = "Correlation"),
        row_gap = unit(0, "mm"),
        border = F)


# With ggplot
# First get fatty acid clusters
ord <- hclust(dist(t(mat)))$order
ord2 <- hclust(dist(mat))$order


library(RColorBrewer)

# borrowing pheatmap's colour palette:
library(RColorBrewer)
phtpalette <- colorRampPalette(rev(brewer.pal(n = 7, name ="RdYlBu")))(100)

ggplot(t1, aes(x = y, y = x, fill = r)) + geom_tile() +
  theme(axis.text.x = element_text(angle = 90, hjust = 1),
        axis.title = element_blank()) +
  scale_fill_gradientn(colours = phtpalette)# +
  #scale_y_discrete(limits = levels(t1$y)[ord]) +
  #scale_x_discrete(limits = levels(t1$Var2)[ord2])


# Correlation matrix (large) (supplemental data)
library(corrplot)
corrplot(cormat, tl.col = "black", tl.cex = 0.3, order = "hclust", method = "color")



# Not in manuscript

#library(pheatmap)
#pheatmap(mat, color = colorRampPalette(rev(brewer.pal(n = 7, name ="RdYlBu")))(100),)

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














