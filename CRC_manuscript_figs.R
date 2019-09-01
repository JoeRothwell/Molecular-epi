source("CRC_prep_data.R")
source("CRC_get_signatures.R")

# Venn diagram for compounds
library(VennDiagram)
venn.diagram(list(Controls = colnames(ctrlA), CRC1 = colnames(ctrlB), CRC2 = colnames(ctrls)), 
             imagetype = "png", filename = "compound_venn.png", 
             height=150, width=150, units="mm", cat.fontfamily = "",
             fontfamily = "")


# Coefficient plots

df1$Class <- factor(df1$class, levels = rev(c("Amino acids", "Biogenic amines", "Monosaccharides", "Acylcarnitines", 
                                              "Lysophosphatidylcholine", "Phosphatidylcholine (acyl-acyl)", "Phosphatidylcholine (acyl-alkyl)", "Sphingolipids")))

library(ggplot2)
ggplot(df1, aes(y = Class, x = Coefficient, colour = Class)) + 
  geom_jitter(height = 0.1) + 
  theme_bw() + geom_vline(xintercept = 0, linetype = "dashed") +
  theme(legend.position = "none", axis.title.y = element_blank()) +
  xlab("Coefficient from first PLS latent variable") +
  ggtitle("A") +
  scale_y_discrete(labels = c("Phosphatidylcholine (acyl-alkyl)" = "Phosphatidylcholine\n(acyl-alkyl)",
                              "Phosphatidylcholine (acyl-acyl)" = "Phosphatidylcholine\n(acyl-acyl)",
                              "Lysophosphatidylcholine" = "Lysophosphatidyl-\ncholine"))

df2$Class <- factor(df2$class, levels = rev(c("Saturated", "Monounsaturated", "Polyunsaturated", "Natural trans", 
                                              "Industrial trans")))

ggplot(df2, aes(y = Class, x = Coefficient, colour = Class)) + 
  geom_jitter(height = 0) + 
  xlab("Coefficient from first PLS latent variable") +
  theme_bw() + geom_vline(xintercept = 0, linetype = "dashed") +
  theme(legend.position = "none", axis.title.y = element_blank()) +
  ggtitle("B")