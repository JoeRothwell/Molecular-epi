# CRC lipids study
# 112 lipids in 3223 participants, no missing values

# Get data from CRC_prep_data_rev and remove unneeded objects
source("CRC_prep_data_rev.R")
rm(list = ls(pattern = "1|2|bmi"))

# Join GRS data
snps <- read_dta("clrt_gwas_gecco_snps_GRS.dta")
csnp <- inner_join(crc, snps, by = "Idepic")
crc1 <- left_join(crc, snps, by = "Idepic")
table(CT = csnp$Cncr_Caco_Clrt, lab = csnp$lab)
table(CT = csnp$Cncr_Caco_Clrt)

library(zoo)

# Continuous models (needs imputation and log transformation)
# Remove non-matched subjects
#crc1 <- crc %>% group_by(Match_Caseset) %>% filter(n() == 2) %>% ungroup()


# Select the lipids only (112) and convert zeros to NA
mat1 <- crc %>% select(Acylcarn_C10:Acylcarn_C8, Glyceroph_Lysopc_A_C16_0:Sphingo_Sm_C26_1) %>% na_if(0)
mat2 <- colon %>% select(Acylcarn_C10:Acylcarn_C8, Glyceroph_Lysopc_A_C16_0:Sphingo_Sm_C26_1) %>% na_if(0)
mat3 <- prox %>% select(Acylcarn_C10:Acylcarn_C8, Glyceroph_Lysopc_A_C16_0:Sphingo_Sm_C26_1) %>% na_if(0)
mat4 <- dist %>% select(Acylcarn_C10:Acylcarn_C8, Glyceroph_Lysopc_A_C16_0:Sphingo_Sm_C26_1) %>% na_if(0)
mat5 <- rectal %>% select(Acylcarn_C10:Acylcarn_C8, Glyceroph_Lysopc_A_C16_0:Sphingo_Sm_C26_1) %>% na_if(0)

# Check for missings
hist(colSums(is.na(mat)), breaks = 30)
scalemat1 <- mat1 %>% na.aggregate(FUN = function(x) min(x)/2) %>% log2 %>% scale
scalemat2 <- mat2 %>% na.aggregate(FUN = function(x) min(x)/2) %>% log2 %>% scale
scalemat3 <- mat3 %>% na.aggregate(FUN = function(x) min(x)/2) %>% log2 %>% scale
scalemat4 <- mat4 %>% na.aggregate(FUN = function(x) min(x)/2) %>% log2 %>% scale
scalemat5 <- mat5 %>% na.aggregate(FUN = function(x) min(x)/2) %>% log2 %>% scale

library(Amelia)
missmap(as_tibble(scalemat))
heatmap.2(scalemat1, trace = "none", col = colpalette1)

library(corrplot)
cormat <- cor(scalemat1)
corrplot(cormat, method = "color", order = "hclust", tl.col = "black", tl.cex = 0.7)

### Continuous models per SD increase
# Define function to apply across quartiles (already matched by lab)
library(survival)
multiclr <- function(x, dat) { 
  clogit(Cncr_Caco_Clrt ~ x + Bmi_C + Qe_Energy + L_School + Smoke_Stat + Smoke_Int + Height_C +
    Qe_Alc + Qge0701 + Qge0704 + strata(Match_Caseset), data = dat) 
  }

library(broom)
# Apply models by subsite (use 2nd fn parameter as an option)
mods1 <- apply(scalemat1, 2, multiclr, dat = crc) %>% map_df( ~ tidy(., exponentiate = T)) %>% 
  filter(grepl("x", term)) %>% mutate(p.adj = p.adjust(p.value, method = "fdr"))

mods2 <- apply(scalemat2, 2, multiclr, dat = colon) %>% map_df( ~ tidy(., exponentiate = T)) %>% 
  filter(grepl("x", term)) %>% mutate(p.adj = p.adjust(p.value, method = "fdr"))

mods3 <- apply(scalemat3, 2, multiclr, dat = prox) %>% map_df( ~ tidy(., exponentiate = T)) %>% 
  filter(grepl("x", term)) %>% mutate(p.adj = p.adjust(p.value, method = "fdr"))

mods4 <- apply(scalemat4, 2, multiclr, dat = dist) %>% map_df( ~ tidy(., exponentiate = T)) %>% 
  filter(grepl("x", term)) %>% mutate(p.adj = p.adjust(p.value, method = "fdr"))

mods5 <- apply(scalemat5, 2, multiclr, dat = rectal) %>% map_df( ~ tidy(., exponentiate = T)) %>% 
  filter(grepl("x", term)) %>% mutate(p.adj = p.adjust(p.value, method = "fdr"))


#crccon <- map_df(mods, ~ tidy(., exponentiate = T)) %>% filter(grepl("x", term)) %>% 
#  mutate_if(is.numeric, ~round(., 4)) %>% mutate(p.adj = p.adjust(p.value, method = "fdr")) %>%
#  unite(OR.CI, estimate, conf.low, conf.high, sep = "-") %>% cbind(cmpd = colnames(scalemat))


### Categorical associations (doesn't need imputation or scaling)
# First need to run CRC_data_prep to get crc and crc3 objects
library(ggplot)
library(broom)

# Biocrates compounds. Split into quartiles with cut_number
mat6 <- mat1 %>% mutate_all(~cut_number(., n = 4, labels = 1:4))
mat7 <- mat2 %>% mutate_all(~cut_number(., n = 4, labels = 1:4))
mat8 <- mat3 %>% mutate_all(~cut_number(., n = 4, labels = 1:4))
mat9 <- mat4 %>% mutate_all(~cut_number(., n = 4, labels = 1:4))
mat10 <- mat5 %>% mutate_all(~cut_number(., n = 4, labels = 1:4))


fits1a <- apply(mat6, 2, multiclr, dat = crc) %>% map_df( ~tidy(., exponentiate = T)) %>% 
  filter(grepl("x", term)) %>% mutate(p.adj = p.adjust(p.value, method = "fdr")) %>%
  mutate_if(is.numeric, ~round(., 3)) #%>%
  #unite(OR.CI, estimate, conf.low, conf.high, sep = "-") %>%
  #cbind(cmpd = colnames(df3)) %>% group_by(cmpd) %>% filter(min(p.value) < 0.05)
  
fits2 <- apply(mat7, 2, multiclr, dat = colon) %>% map_df( ~tidy(., exponentiate = T)) %>% 
  filter(grepl("x4", term)) %>% mutate(p.adj = p.adjust(p.value, method = "fdr")) %>%
  mutate_if(is.numeric, ~round(., 3))
  
fits3 <- apply(mat8, 2, multiclr, dat = prox) %>% map_df( ~tidy(., exponentiate = T)) %>% 
  filter(grepl("x4", term)) %>% mutate(p.adj = p.adjust(p.value, method = "fdr")) %>%
  mutate_if(is.numeric, ~round(., 3))

fits4 <- apply(mat9, 2, multiclr, dat = dist) %>% map_df( ~tidy(., exponentiate = T)) %>% 
  filter(grepl("x4", term)) %>% mutate(p.adj = p.adjust(p.value, method = "fdr")) %>%
  mutate_if(is.numeric, ~round(., 3)) 

fits5 <- apply(mat10, 2, multiclr, dat = rectal) %>% map_df( ~tidy(., exponentiate = T)) %>% 
  filter(grepl("x4", term)) %>% mutate(p.adj = p.adjust(p.value, method = "fdr")) %>%
  mutate_if(is.numeric, ~round(., 3))


  
# Polygenic risk scores
snps <- read_dta("clrt_gwas_gecco_snps_GRS.dta")
csnp <- inner_join(crc, snps, by = "Idepic")
table(CT = csnp$Cncr_Caco_Clrt, lab = csnp$lab)
table(CT = csnp$Cncr_Caco_Clrt)
plot(csnp$GRS)
csnp$GRScat <- cut_number(csnp$GRS, n = 4, labels = 1:4)

table(CT = csnp$Cncr_Caco_Clrt, Q = csnp$GRScat)
chisq.test(csnp$Cncr_Caco_Clrt, csnp$GRScat)

ggplot(csnp, aes(x = GRS, group = as.factor(Cncr_Caco_Clrt))) + 
  geom_density(aes(fill = as.factor(Cncr_Caco_Clrt)), alpha = 0.4)

# Metabolite associations with covariates
# Subset control matrix
ints <- mat1[crc$Cncr_Caco_Clrt == 0, ]
meta1 <- crc1[crc1$Cncr_Caco_Clrt == 0, ]



lm1 <- function(x) lm(as.numeric(Fasting_C) ~ x + Sex + lab, data = meta1)
lm2 <- function(x) lm(as.numeric(Sex) ~ x + Fasting_C + lab, data = meta1)
lm3 <- function(x) lm(as.numeric(lab) ~ x + Fasting_C + Sex, data = meta1)
lm4 <- function(x) lm(Qe_Alc ~ x + Fasting_C + Sex + lab, data = meta1)
lm5 <- function(x) lm(as.numeric(Smoke_Stat) ~ x + Fasting_C + Sex + lab, data = meta1)
lm6 <- function(x) lm(as.numeric(Smoke_Int) ~ x + Fasting_C + Sex + lab, data = meta1)
lm7 <- function(x) lm(Age_Blood ~ x + Bmi_C + Fasting_C + Sex + lab, data = meta1)
lm8 <- function(x) lm(Bmi_C ~ x + Fasting_C + Sex + lab, data = meta1)
lm9 <- function(x) lm(Waist_C ~ x + Fasting_C + Sex + lab, data = meta1)
lm10 <- function(x) lm(Qge0701 ~ x + Fasting_C + Sex + lab, data = meta1)
lm11 <- function(x) lm(Qge0704 ~ x + Fasting_C + Sex + lab, data = meta1)
lm12 <- function(x) lm(GRS ~ x + Fasting_C + Sex + lab, data = meta1)

library(broom)
fits1 <- apply(ints, 2, lm1) %>% map_df(tidy) #%>% filter(str_detect(term, "x"))
fits2 <- apply(ints, 2, lm2) %>% map_df(tidy) #%>% filter(str_detect(term, "x"))
fits3 <- apply(ints, 2, lm3) %>% map_df(tidy) #%>% filter(str_detect(term, "x"))
fits4 <- apply(ints, 2, lm4) %>% map_df(tidy) #%>% filter(str_detect(term, "x"))
fits5 <- apply(ints, 2, lm5) %>% map_df(tidy) #%>% filter(str_detect(term, "x"))
fits6 <- apply(ints, 2, lm6) %>% map_df(tidy) #%>% filter(str_detect(term, "x"))
fits7 <- apply(ints, 2, lm7) %>% map_df(tidy) #%>% filter(str_detect(term, "x"))
fits8 <- apply(ints, 2, lm8) %>% map_df(tidy) #%>% filter(str_detect(term, "x"))
fits9 <- apply(ints, 2, lm9) %>% map_df(tidy) #%>% filter(str_detect(term, "x"))
fits10 <- apply(ints, 2, lm10) %>% map_df(tidy) #%>% filter(str_detect(term, "x"))
fits11 <- apply(ints, 2, lm11) %>% map_df(tidy) #%>% filter(str_detect(term, "x"))
fits12 <- apply(ints, 2, lm12) %>% map_df(tidy) #%>% filter(str_detect(term, "x"))

all <- bind_rows(
                 #fits1, fits2, fits3, 
                 fits4, fits5, fits6, fits7, fits8, fits9, fits10,
                 fits11, fits12) %>% filter(term == "x") %>%
  mutate(p.adj = p.adjust(p.value, method = "fdr")) %>%
  bind_cols(cmpd = rep(colnames(ints), 9))

library(RColorBrewer)
# Different colour palettes
set2 <- rep(brewer.pal(8, "Set2"), each = 112, length.out = nrow(all))
accent <- rep(brewer.pal(8, "Accent") , each = 112, length.out = nrow(all))
rain <- rep(rainbow(12), each = 112, length.out = nrow(all))

# Plot (add col palette as necessary)
#plot(-log10(all$p.value), col = set2, pch = 19, cex = 0.6)

# x axis width
x = 1:nrow(all)

# draw empty plot
plot(NULL, xlim=c(0, nrow(all)), ylim=c(0, max(-log10(all$p.value))), xaxt='n',
     ylab='-log10(p-value)', xlab='')
points(x, -log10(all$p.value), pch=19, col=set2, cex = 0.6)

# axis labels
labs <- c("Alc", "Smoke", "Smoke", "Age", "BMI", "Waist", "Red meat", "Proc meat", 
          "Genetic risk")
abline(h = -log(0.05), col = "grey")
abline(h = -log(0.016), col = "red")
axis(1, at = c(21, 64, 107, 150, 193, 236, 279, 322, 365, 408, 451, 494), labels = labs, 
     las=1, cex.axis = 0.7)
text(350, 22.5, "Choline", pos = 2, cex = 0.7)
text(360, 21.5, "Glycerophosphocholine", pos = 2, cex = 0.7)

