#ob.m <- read.csv("D:/obes_metabo.csv") 
#saveRDS(ob.m, "Biocrates data controls.rds")
ob.m <- readRDS("Biocrates data controls.rds")
cof <- readRDS("Coffee intakes all EPIC.rds")

# 3773 observations
# Wcrf_C_Cal: Total WCRF/AICR score (categ.) (calibrated)
# Wcrf_Cal: Total WCRF/AICR score (calibrated)
# Model scores from Biocrates data ----

library(tidyverse)

# code 0 and 1 as categories 1 and 2 vs 4 and 5, remove cat 3
ob <- ob.m %>% filter(Wcrf_C != 3) %>% mutate(score_cat = ifelse(Wcrf_C %in% c(1,2), 0, 1))
table(ob$score_cat)

# Make matrix of metabolomics data. 146 compounds
metabo <- ob %>% select(Acylcarn_C10:Sugars_H1) %>% as.matrix
library(zoo)
metabo1 <- na.aggregate(metabo, FUN = function(x) min(x)/2)
metabo.log <- log2(metabo1)

dim(metabo)
# Check missings
missings <- apply(metabo, 2, function(x) sum(is.na(x)))
which(missings > 1500)
library(zoo)
metabo1 <- na.aggregate(metabo, FUN = function(x) min(x)/2)

# Subset only covariates needed from data
ob <- ob %>% select(score_cat, Sex, Center, Study, Smoke_Intensity) %>% as.tibble
ob$Sex <- droplevels(ob$Sex)
ob$Study <- droplevels(ob$Study)

# fn for logistic regression on categories 1 and 2 vs 4 and 5
glm.ob <- function(x) glm(score_cat ~ x + Sex + Center + Smoke_Intensity + Study,
                          data = ob, family = "binomial")

# apply logistic regression to each metabolite
multi.fit <- apply(metabo1, 2, glm.ob)
# extract p-values
p.all <- sapply(multi.fit, function(f) summary(f)$coefficients[ , 4])
p <- p.all[2 , ] 
p.adj <- p.adjust(p, method = "BH")

# calculation of fold changes for volcano plot
df <- data.frame(ob$score_cat, metabo1)
means <- aggregate(. ~ ob.score_cat, data = df, mean) %>% t
# compute fold change to have higher fold change for higher score
mean.fc <- means[, 2]/means[, 1]

volcanodf <- data_frame(Cmpd = colnames(metabo1), Fold_change = mean.fc[-1], p.value = p.adj) %>%
  mutate(points2label = ifelse(p.adj < 0.001 & abs(1 - Fold_change) > 0.05, Cmpd, NA))

ggplot(volcanodf, aes(x=Fold_change, y = -log10(p.value))) + geom_point() + theme_bw() +
  xlim(0.64, 1.34) +
  geom_text(aes(label = points2label), hjust = -0.1, vjust = 0, size = 2)

# Correlation with coffee intakes ----------------------------------------------------------------

ob.caf <- ob.m %>% left_join(cof, by="Idepic")
metabo <- ob.caf %>% select(Acylcarn_C10:Sugars_H1) %>% as.matrix
cormat <- cor(log2(ob.caf$QGE130301 + 1), metabo, use = "pairwise.complete.obs")
log.cormat <- cor(log2(ob.caf$QGE130301 + 1), metabo, use = "pairwise.complete.obs")
