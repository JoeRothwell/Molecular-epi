source("CRC_prep_data.R")
source("CRC_get_signatures.R")

# Venn diagram for compounds
library(VennDiagram)
venn.diagram(list(Controls = colnames(ctrlA), CRC1 = colnames(ctrlB), CRC2 = colnames(ctrls)), 
             imagetype = "png", filename = "compound_venn.png", 
             height=150, width=150, units="mm", cat.fontfamily = "",
             fontfamily = "")


# Coefficient plots

pltdata$Class <- factor(pltdata$class, 
              levels = rev(c("Amino acids", "Biogenic amines", "Monosaccharides", "Acylcarnitines", 
              "Lysophosphatidylcholine", "Phosphatidylcholine (acyl-acyl)", 
              "Phosphatidylcholine (acyl-alkyl)", "Sphingolipids")))

library(ggplot2)
p1 <- ggplot(pltdata, aes(y = Class, x = Coefficient, colour = Class)) + 
  geom_jitter(height = 0.1) + 
  theme_bw() + geom_vline(xintercept = 0, linetype = "dashed") +
  theme(legend.position = "none", axis.title.y = element_blank()) +
  xlab("Coefficient from first PLS latent variable") +
  ggtitle("A") +
  scale_y_discrete(labels = c("Phosphatidylcholine (acyl-alkyl)" = "Phosphatidylcholine\n(acyl-alkyl)",
                              "Phosphatidylcholine (acyl-acyl)" = "Phosphatidylcholine\n(acyl-acyl)",
                              "Lysophosphatidylcholine" = "Lysophosphatidyl-\ncholine"))

faplot$Class <- factor(faplot$class, levels = rev(c("Saturated", "Monounsaturated", 
                  "Polyunsaturated", "Natural trans", "Industrial trans")))

p2 <- ggplot(faplot, aes(y = Class, x = Coefficient, colour = Class)) + 
  geom_jitter(height = 0) + 
  xlab("Coefficient from first PLS latent variable") +
  theme_half_open() + geom_vline(xintercept = 0, linetype = "dashed") +
  theme(legend.position = "none", axis.title.y = element_blank()) +
  ggtitle("B")

# Generate aligned plots
library(cowplot)
both <- align_plots(p1, p2, align = "hv", axis = "tblr")
p1x <- ggdraw(both[[1]])
p2x <- ggdraw(both[[2]])



# Biocrates and FAs together in facetted plot (couldn't put classes in right order)
all <- bind_rows(Endogenous = pltdata, `Fatty acids` = faplot, .id = "Metabolites")

ggplot(all, aes(y = Class, x = Coefficient, colour = Class)) + 
  geom_jitter(height = 0.05) + 
  xlab("Coefficient from first PLS latent variable") +
  theme_bw() + 
  geom_vline(xintercept = 0, linetype = "dashed") +
  theme(legend.position = "none", axis.title.y = element_blank(),
        axis.line.y = element_line()) +
  facet_grid(Metabolites ~ ., scales = "free", space = "free", switch = "y")



# Old coefficient scatter (top and bottom percentiles)

# Vector of black and grey for plot points
vec <- c( rep("black", n.low), rep("grey", nrow(dat) - nrow(dat1)), rep("black", n.high) )

# Now plot data, adding text
plot(sort(coeff$value), pch = 17, col=vec, xlab = "", ylab = "Coefficient",
     main = paste(nrow(mod$scores), "fasted subjects, optimal dimensions =", lv))
# High and low labels
text(nrow(dat) : (nrow(dat) - n_top), df1$Coefficient, df1$compound, pos=2, cex = 0.6)
text(1:nrow(df2), df2$Coefficient, df2$compound, pos=4, cex=0.6)
abline(a=0, b=0, lty = "dotted")


