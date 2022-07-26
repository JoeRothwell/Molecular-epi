# Smile plots for AA associations

# EPIC
source("CRC_aminoacid_pooled.R")
ggplot(results, aes(x = estimate, y = -log10(p.value))) + geom_point() + 
  theme_bw() + geom_hline(yintercept = -log10(0.05), linetype = "dashed")

# UK Biobank (need to get results from spreadsheet)

pval <- c(ala = 0.597,
gln = 0.087,
gly = 0.404,
hist = 0.028,
iso = 0.844,
leu = 0.342,
val = 0.653,
phe = 0.586,
tyr = 0.754)

p.adjust(pval, method = "fdr")
