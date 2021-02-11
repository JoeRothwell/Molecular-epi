# CRC amino acids study. Get data from CRC_prep_data and remove unneeded objects
source("CRC_prep_data.R")
rm(list = ls(pattern = "crc|dist|prox|rect"))

# Subset whole colon study. Get compounds and filter those with more than 31% NAs
mat0 <- colon %>% select(contains("Aminoacid_")) %>% select_if(~ sum(is.na(.)) < (nrow(colon)*0.31)) 

# Separate studies for Jelena's paper
# p180 (small) had 740 subjects, p150 (large) had 556
mat1 <- colon1 %>% select(colnames(mat0)) %>% na_if(0) #738, 698 w/o Greece
mat2 <- colon2 %>% select(colnames(mat0)) %>% na_if(0) #564, 556 w/o Greece

library(zoo)
scalemat1 <- mat1 %>% na.aggregate(FUN = function(x) min(x)/2) %>% scale
scalemat2 <- mat2 %>% na.aggregate(FUN = function(x) min(x)/2) %>% scale

# Plot distributions
plot.ts(mat1[, 1:6], type = "p", main = "Arg to Orn")
plot.ts(mat1[, 7:13], type = "p", main = "Phe to Val")

### Continuous models per SD increase
# Define function to apply across quartiles (already matched by lab)

library(survival)
multiclr <- function(x, dat) { 
  clogit(Cncr_Caco_Clrt ~ x + Bmi_C + Qe_Energy + L_School + Smoke_Stat + Smoke_Int + Height_C +
           Qe_Alc + Qge0701 + Qge0704 + strata(Match_Caseset), data = dat) 
}

# Original co-variates

multiclr <- function(x, dat) { 
  clogit(Cncr_Caco_Clrt ~ x + Bmi_Cat + Smoke_Stat + Alc_Drinker + Pa_Total + 
           strata(Match_Caseset), data = dat) 
}


library(broom)
mods1 <- apply(scalemat1, 2, multiclr, dat = colon1) %>% 
  map_df( ~ tidy(., exponentiate = T, conf.int = T)) %>% filter(grepl("x", term)) %>% 
  mutate(p.adj = p.adjust(p.value, method = "fdr"), compound = colnames(mat1))

mods2 <- apply(scalemat2, 2, multiclr, dat = colon2) %>% 
  map_df( ~ tidy(., exponentiate = T, conf.int = T)) %>% filter(grepl("x", term)) %>% 
  mutate(p.adj = p.adjust(p.value, method = "fdr"), compound = colnames(mat2))

#pFDR <- mods1 %>% filter(p.adj <= 0.05) %>% select(p.value) %>% max
#pFDR <- mods1 %>% filter(p.adj <= 0.05) %>% select(p.value) %>% max

# Vector of ORs to paste into Excel sheet
paste(round(mods1$estimate, 2), " (", round(mods1$conf.low, 2), "-", 
      round(mods1$conf.high, 2), ")", sep = "") %>% as.tibble()

paste(round(mods2$estimate, 2), " (", round(mods2$conf.low, 2), "-", 
      round(mods2$conf.high, 2), ")", sep = "") %>% as.tibble()

# Meta analysis by nesting
library(metafor)
cont <- bind_rows(mods1, mods2) %>% group_by(compound) %>% nest() %>% 
  mutate(mods = lapply(data, function(df) rma(estimate, sei = std.error, data=df, method="REML")))
ma.cont <- lapply(cont$mods, "[", c("b", "ci.lb", "ci.ub", "se", "I2", "QEp"))
results1 <- map_df(ma.cont, bind_rows) %>% mutate_all(~round(., 2))


### Categorical analysis. Need to alter function to get categories with cutpoints based on controls
ctrl <- colon1$Cncr_Caco_Clrt == 0
cutct <- function(x) cut(x, breaks = quantile(x[ctrl]), include.lowest = T, labels = 1:4)
mat3 <- apply(mat1, 2, cutct)

fits1 <- apply(mat3, 2, multiclr, dat = colon1) %>% 
  map_df( ~tidy(., exponentiate = T, conf.int = T)) %>%
  filter(grepl("x", term)) %>% 
  mutate(compound = rep(colnames(mat1), each = 3))
# mutate(p.adj = p.adjust(p.value, method = "fdr"), compound = colnames(mat1))

hrci <- paste(round(fits1$estimate, 2), " (", round(fits1$conf.low, 2), "-", 
        round(fits1$conf.high, 2), ")", sep = "") %>% as.tibble()

# Get p-trend by entering quartile cutoffs into the model as continuous variables
library(Hmisc)
cutpts <- function(x) as.numeric(as.character(cut2(x, cuts = quantile(x[ctrl]), levels.mean = T)))
mat3b <- map_df(mat1, cutpts)
trends <- apply(mat3b, 2, multiclr, dat = colon1) %>% map_df( ~ tidy(.)) %>% filter(grepl("x", term))

#mat4 <- mat2 %>% mutate_all(~cut_number(., n = 4, labels = 1:4))
cutct <- function(x) cut(x, breaks = quantile(x[colon2$Cncr_Caco_Clrt == 0]), include.lowest = T, labels = 1:4)
mat4 <- apply(mat2, 2, cutct)

fits2 <- apply(mat4, 2, multiclr, dat = colon2) %>% 
  map_df( ~tidy(., exponentiate = T, conf.int = T)) %>%
  filter(grepl("x", term)) %>% 
  mutate(compound = rep(colnames(mat1), each = 3))
  mutate(p.adj = p.adjust(p.value, method = "fdr"), compound = colnames(mat2))
  
hrci <- paste(round(fits1$estimate, 2), " (", round(fits1$conf.low, 2), "-", 
        round(fits1$conf.high, 2), ")", sep = "") %>% as.tibble() #%>%
        #mutate(quartile = rep(c("X2", "X3", "X4"), 13)) %>%
        #pivot_wider(names_from = quartile, values_from = value)

ctrl <- colon2$Cncr_Caco_Clrt == 0
cutpts <- function(x) as.numeric(as.character(cut2(x, cuts = quantile(x[ctrl]), levels.mean = T)))
mat4b <- map_df(mat2, cutpts)
trends <- apply(mat4b, 2, multiclr, dat = colon2) %>% map_df( ~ tidy(.)) %>% filter(grepl("x", term))

# Meta-analysis
cats <- bind_rows(fits1, fits2) %>% group_by(compound, term) %>% nest() %>% 
  mutate(mods = lapply(data, function(df) rma(estimate, sei = std.error, data=df, method="REML")))
ma.cat <- lapply(cats$mods, "[", c("b", "ci.lb", "ci.ub", "se", "I2", "QEp"))
results2 <- map_df(ma.cat, bind_rows) %>% mutate_all(~round(., 2)) %>% 
  mutate(quartile = cats$term) %>% arrange(quartile)

# Bind continuous and categorical results together
p180 <- bind_rows(Continuous = mods1, Categorical = fits1, .id = "Analysis")
p150 <- bind_rows(Continuous = mods2, Categorical = fits2, .id = "Analysis")

library(ggplot2)
ggplot(df)

ggplot(all, aes(x = estimate, y = fct_inorder(rev(group)), 
               colour = Append, xmin = ci.low, xmax = ci.high)) + 
  geom_pointrange() + theme_bw() +
  geom_errorbar(aes(xmin=ci.low, xmax=ci.high), width=0.5, cex=1) + 
  ylab('Group') + xlab("Geometric mean biomarker concentration") +
  facet_grid(biomarker ~ ., scales = "free_y") +
  ggtitle("Biomarker measurements in appendectomy\nand non-appendectomy groups") +
  theme(axis.title.y =element_blank())



