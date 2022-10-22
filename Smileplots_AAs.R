# Smile plots for AA associations

# EPIC
source("CRC_aminoacid_pooled.R")
library(tidyverse)
library(ggplot2)
library(ggrepel)

# Get data for EPIC and UK Biobank and split columns
results1 <- results %>% separate(compound, into = c("aminoacid", "compound"), sep = "_")

AAresults <- read.csv("AAresults.csv", header = F)
colnames(AAresults) = c("compound", "estimate", "se", "coeff", "p.value", "conf.low", "conf.high")

# Make ggplots
plotA <- ggplot(results1, aes(x = estimate, y = log10(p.value))) + geom_point() + 
  scale_y_reverse(limits = c(0, -3), breaks = c(0:-2), labels = function(x) 10^x) +
  theme_bw() + geom_hline(yintercept = log10(0.05), linetype = "dashed") +
  geom_vline(xintercept = 1, size = 0.2, colour = "grey60") + 
  xlab("OR per SD concentration") + ylab("Raw P-value") +
  theme(panel.grid.minor = element_blank(), panel.grid.major = element_blank()) +
  geom_text(aes(label = compound), vjust = 1.5) +
  labs(title = "EPIC")

plotB <- ggplot(AAresults, aes(x = estimate, y = log10(p.value))) + geom_point() + 
  scale_y_reverse(limits = c(0, -2), breaks = c(0:-2), labels = function(x) 10^x) +
  theme_bw() + geom_hline(yintercept = log10(0.05), linetype = "dashed") +
  geom_vline(xintercept = 1, size = 0.2, colour = "grey60") + 
  xlab("HR per SD concentration") + ylab("Raw P-value") +
  theme(panel.grid.minor = element_blank(), panel.grid.major = element_blank()) +
  geom_text(aes(label = compound), vjust = 1.5) +
  labs(title = "UK Biobank")

library(cowplot)
plot_grid(plotA, plotB, ncol = 2)


# UK Biobank p-value adjustment (entered results manually from spreadsheet)

# p-values UK Biobank 2021
pval <- c(ala = 0.597, gln = 0.087, gly = 0.404, hist = 0.028, iso = 0.844,
  leu = 0.342, val = 0.653, phe = 0.586, tyr = 0.754)

p.adjust(pval, method = "fdr")
#   ala   gln   gly  hist   iso   leu   val   phe  tyr 
# 0.840 0.392 0.840 0.252 0.844 0.840 0.840 0.840 0.84 

# Adjusted pvalues UK biobank 2022 update
pval <- c(ala = 0.45, gln = 0.086, gly = 0.389, hist = 0.019, iso = 0.819,
  leu = 0.153, val = 0.421, phe = 0.732, tyr = 0.646)

p.adjust(pval, method = "fdr")
#   ala   gln   gly  hist   iso   leu   val   phe   tyr 
# 0.675 0.387 0.675 0.171 0.819 0.459 0.675 0.819 0.819 