source("CRC_prep_data.R")
source("CRC_get_signatures.R")

# Coefficient plots
library(ggplot2)
library(ggrepel)

# To change legend order, reorder factor levels
corr <- cor(Bioc0[, -1], Bioc0[, 1])
pltdata1 <- pltdata %>% mutate(compound1 = ifelse(abs(Coefficient) > 0.022, compound, NA), corr)

p1 <- 
  ggplot(pltdata1, aes(x = Coefficient, y = corr, shape = str_wrap(class, 20))) + 
  geom_point() + theme_bw() +
  theme(panel.grid.major = element_blank(), legend.title = element_blank(),
        legend.position = c(0.85, 0.25)) +
  scale_shape_manual(values = c(15, 1, 19, 25, 22, 14, 3, 4)) +
  xlab("Coefficient 1st PLSR latent variable") + 
  ylab("Correlation WCRF score - metabolite") +
  geom_hline(yintercept = 0, linetype = "dashed") +
  geom_vline(xintercept = 0, linetype = "dashed") +
  geom_text_repel(aes(label = compound1), size = 3) #+
  #ggtitle("A")

corr1 <- cor(FAdata[, -1], FAdata[, 1])
faplot1 <- faplot %>% mutate(compound1 = ifelse(abs(Coefficient) > 0.045, compound, NA), corr1)

p2 <- ggplot(faplot1, aes(x = Coefficient, y = corr1, shape = class)) + geom_point() + 
  theme_bw() + xlim(-0.19, 0.19) + ylim(-0.22, 0.22) +
  theme(panel.grid.major = element_blank(), legend.title = element_blank(),
        legend.position = c(0.85, 0.2)) +
  scale_shape_manual(values = c(6, 16, 1, 4, 3)) +
  xlab("Coefficient 1st PLSR latent variable") + 
  ylab("Correlation WCRF score - metabolite") +
  geom_hline(yintercept = 0, linetype = "dashed") +
  geom_vline(xintercept = 0, linetype = "dashed") +
  geom_text_repel(aes(label = compound1), size = 3) #+
  #ggtitle("B")

# Generate aligned plots
library(cowplot)
plot_grid(p1, p2, nrow = 2, labels = c("A", "B"), label_size = 12) %>%
  save_plot("s_plots.png")



both <- align_plots(p1, p2, align = "hv", axis = "tblr")
ggdraw(both[[1]])
ggsave("endogenous.png", height = 100, width = 180, units = "mm")
ggdraw(both[[2]])
ggsave("FAs.png", height = 100, width = 180, units = "mm")

# Venn diagram for compounds
library(VennDiagram)
venn.diagram(list(Controls = colnames(ctrlA), CRC1 = colnames(ctrlB), CRC2 = colnames(ctrls)), 
             imagetype = "png", filename = "compound_venn.png", 
             height=150, width=150, units="mm", cat.fontfamily = "", fontfamily = "")

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
