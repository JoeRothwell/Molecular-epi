# Amino acid correlations in EPIC, manuscript figure 1

source("CRC_prep_data.R")
ctrl <- colon$Cncr_Caco_Clrt == 0

# For 13 P150 AAs only
mat0 <- colon %>% select(contains("Aminoacid_")) %>% select_if(~ sum(is.na(.)) < (nrow(colon)*0.31))

# For 21 P150+180 AAs (removing sum of leucine + isoleucine)
mat0 <- colon %>% select(contains("Aminoacid_"), -Aminoacid_Xleu) 

mat <- colon %>% select(colnames(mat0)) 
mat1 <- mat[ctrl, ]

# Remove Aminoacid from labels. Any character before _
colnames(mat1) <- str_replace(colnames(mat1), ".*_", "")
# Aminoacid before _
colnames(mat1) <- str_replace(colnames(mat1), "Aminoacid_", "")


#Full names
colnames(mat1) <- c("Alanine", "Arginine", "Asparagine", "Aspartate", "Citrulline", "Glutamate", "Glutamine", 
                    "Glycine", "Histidine", "Isoleucine", "Leucine", "Lysine", "Methionine", "Ornithine", 
                    "Phenylalanine", "Proline", "Serine", "Threonine", "Tryptophan", "Tyrosine", "Valine")

cormat <- cor(mat1, use = "pairwise.complete.obs")
library(corrplot)
cp <- corrplot(cormat, method = "square", tl.col = "black", order = "hclust", addrect = 10)

# Get names in cluster order
#dimnames(cp$corr)[[1]]


# Get mean concentrations for each amino acid
# First add proper names
colnames(mat0) <- c("Alanine", "Arginine", "Asparagine", "Aspartate", "Citrulline", "Glutamate", "Glutamine", 
                    "Glycine", "Histidine", "Isoleucine", "Leucine", "Lysine", "Methionine", "Ornithine", 
                    "Phenylalanine", "Proline", "Serine", "Threonine", "Tryptophan", "Tyrosine", "Valine")

AAmeans <- apply(mat0, 2, mean, na.rm = T)
AAsd <- apply(mat0, 2, sd, na.rm = T)

data <- data.frame(AAmeans)
data$sd <- AAsd
data1 <- rownames_to_column(data, "compound")

# Put in cluster order
AAord <- AAmeans[rev(c("Aspartate", "Glutamate", "Glycine", "Serine", "Tyrosine", "Arginine", 
                       "Methionine", "Tryptophan", "Alanine", "Proline", "Lysine", "Valine", 
                       "Isoleucine", "Leucine", "Ornithine", "Histidine", "Phenylalanine", "Glutamine", 
                       "Citrulline", "Asparagine", "Threonine"))]
AAsd <- apply(mat0, 2, sd, na.rm = T)
AAord1 <- AAsd[rev(c("Aspartate", "Glutamate", "Glycine", "Serine", "Tyrosine", "Arginine", 
                     "Methionine", "Tryptophan", "Alanine", "Proline", "Lysine", "Valine", 
                     "Isoleucine", "Leucine", "Ornithine", "Histidine", "Phenylalanine", "Glutamine", 
                     "Citrulline", "Asparagine", "Threonine")) ]

# Old names
#AAord <- AAmeans[rev(c("Aminoacid_Asp", "Aminoacid_Glu", "Aminoacid_Gly", "Aminoacid_Ser", "Aminoacid_Tyr", "Aminoacid_Arg", 
#                   "Aminoacid_Met", "Aminoacid_Trp", "Aminoacid_Ala", "Aminoacid_Pro", "Aminoacid_Lys", "Aminoacid_Val", 
#                   "Aminoacid_Ile", "Aminoacid_Leu", "Aminoacid_Orn", "Aminoacid_His", "Aminoacid_Phe", "Aminoacid_Gln", 
#                   "Aminoacid_Cit", "Aminoacid_Asn", "Aminoacid_Thr"))]
#AAsd <- apply(mat0, 2, sd, na.rm = T)
#AAord1 <- AAsd[rev(c("Aminoacid_Asp", "Aminoacid_Glu", "Aminoacid_Gly", "Aminoacid_Ser", "Aminoacid_Tyr", "Aminoacid_Arg", 
#                 "Aminoacid_Met", "Aminoacid_Trp", "Aminoacid_Ala", "Aminoacid_Pro", "Aminoacid_Lys", "Aminoacid_Val", 
#                 "Aminoacid_Ile", "Aminoacid_Leu", "Aminoacid_Orn", "Aminoacid_His", "Aminoacid_Phe", "Aminoacid_Gln", 
#"Aminoacid_Cit", "Aminoacid_Asn", "Aminoacid_Thr")) ]

AAmeans - AAsd
AAmeans + AAsd

# Histogram
library(gplots)
barplot2(AAord, horiz = T, xlab = "Mean concentration (SD), uM)", axisnames = F, ci.l = AAord-AAord1, ci.u = AAord+AAord1, 
         plot.ci = T)

data <- data.frame(name=names(AAord),value=AAord, sd=AAord1)
data <- data.frame(name=names(AAord),value=AAord, sd=AAord1)[order(AAord), ]

AAs <- c("Threonine", "Asparagine", "Citrulline", "Glutamine", "Phenylalanine", "Histidine",
         "Ornithine", "Leucine", "Isoleucine", "Valine", "Lysine", "Proline", "Alanine", "Tryptophan",
         "Methionine", "Arginine", "Tyrosine", "Serine", "Glycine", "Glutamate", "Aspartate")

data1 <- data.frame(name=names(AAord),value=AAord, sd=AAord1)[order(AAord), ]

data1$cmpd <- AAs

library(ggplot2)
ggplot(data1) +
  geom_bar(aes(y=reorder(compound, AAmeans), x=AAmeans), stat="identity", colour = "black", fill = "lightblue") + 
  theme_bw() +
  geom_errorbar(aes(y=compound, xmin=AAmeans-sd, xmax=AAmeans+sd), width = 0.5) +
  xlab(expression(paste("Mean concentration (SD), ", mu, "M"))) +
  theme_classic() +
  theme(panel.grid.minor = element_blank(), axis.title.y = element_blank())




# Using complex heatmap (needs more work)
library(ComplexHeatmap)
Heatmap(t(mat1), 
        #column_title = "Fatty acids", 
        row_title = "Amino acid", 
        row_title_gp = gpar(fontsize = 11), column_title_gp = gpar(fontsize = 11),
        row_title_side = "left", column_title_side = "top",
        cluster_rows = T,
        #row_title = NULL,
        row_names_side = "right", 
        show_row_names = T, # affects the spacing of annotations
        row_names_gp = gpar(fontsize = 10),
        column_names_gp = gpar(fontsize = 10),
        #left_annotation = row_ha,
        #right_annotation = ha,
        #bottom_annotation = ha,
        heatmap_legend_param = list(title = "Concentration"),
        row_gap = unit(0, "mm"),
        border = F)
