# Barplots alcohol intake and ethanol
library(tidyverse)

# Read 1623 observations of 44 intensity variables (appears to be final scaled data) and metadata
ints <- read_tsv("1507_XMetabolite_std_cpmg_E3N.txt")

# Lifestyle data. Subset variables needed
meta <- read_csv("Lifepath_meta.csv", na = "9999")

data <- left_join(meta, ints, by = "CODBMB")

mat <- data %>% group_by(CT, MENOPAUSE) %>% summarise(alcohol = mean(ALCOHOL)) %>% 
  spread(MENOPAUSE, alcohol) %>% as.matrix

barplot2(mat[ , -1], beside = T, col = c("grey12", "grey82"), #plot.ci = T,
         ylim = c(0, 20),
         names.arg = c("Pre-menopausal", "Post-menopausal"), legend = c("Controls", "Cases"))

pval.pre <- wilcox.test(ALCOHOL ~ CT, data = data, subset = MENOPAUSE == 0)$p.value
pval.post <- wilcox.test(ALCOHOL ~ CT, data = data, subset = MENOPAUSE == 1)$p.value







boxplot(data$Ethanol ~ data$CT + data$MENOPAUSE, varwidth = T, outline = F,
        names = c("Control, pre", "Case, pre", "Control, post", "Case, post"),
        col = "dodgerblue", ylab = "Plasma ethanol conc (scaled)")



mat <- data %>% group_by(CT, MENOPAUSE) %>% summarise(ethanol = mean(Ethanol)) %>% 
  spread(MENOPAUSE, ethanol) %>% as.matrix

barplot2(mat[ , -1], beside = T, col = c("grey12", "grey82"), #plot.ci = T,
         #ylim = c(0, 20),
         names.arg = c("Pre-menopausal", "Post-menopausal"), legend = c("Controls", "Cases"))
