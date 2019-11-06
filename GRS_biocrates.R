source("CRC_prep_data.R")
# Associations between metabolites and genetic risk score of CRC

# Source data is allele dosages in Stata file clrt_gwas_gecco_snps.dta. Non-integers are probably imputed SNPs. 
# PRS were first calculated from Todd's STATA file CRC_GRS_Rothwell.do and the output named clrt_gwas_gecco_snps_GRS.dta

library(haven)
snps <- read_dta("clrt_gwas_gecco_snps_GRS.dta")

prep.data <- function(dat, group = 0, metabs = T) {

  library(tidyverse)
  
  # Find which observations have SNP data
  intersect(dat$Idepic, crc1$Idepic) %>% length
  intersect(dat$Idepic, crc2$Idepic) %>% length
  
  # Join small CRC to SNPs
  crc1 <- inner_join(crc1, dat, by="Idepic")
  crc2 <- inner_join(crc2, dat, by="Idepic")
  # table(crc1$Cncr_Caco_Clrt)
  # 444 controls 442 cases in small, 296 controls 439 cases in large
  
  # For either case-control model or associations with metabolites
  if(metabs == F) {
    crc1 <- crc1 %>% group_by(Match_Caseset) %>% filter(n() == 2) %>% mutate(study = "small") %>% ungroup()
    crc2 <- crc2 %>% group_by(Match_Caseset) %>% filter(n() == 2) %>% mutate(study = "large") %>% ungroup()
  } else {
    crc1 <- crc1 %>% filter(Cncr_Caco_Clrt == group) %>% mutate(study = "small")
    crc2 <- crc2 %>% filter(Cncr_Caco_Clrt == group) %>% mutate(study = "large")
  }
  
  # subset Biocrates data from the two datasets and arrange by the common columns
  mat1 <- crc1 %>% select(matches("Acylcarn_|Aminoacid_|Biogenic_|Glyceroph_|Sphingo_|Sugars_"), 
                 -starts_with("Outdq"))
  mat2 <- crc2 %>% select(matches("Acylcarn_|Aminoacid_|Biogenic_|Glyceroph_|Sphingo_|Sugars_"), 
                 -starts_with("Outdq"))
  
  common.cols <- intersect(colnames(mat1), colnames(mat2))
  
  # Prepare matrix of metabolite concentrations to apply model across
  mat1a <- crc1 %>% select(one_of(common.cols))
  mat2a <- crc2 %>% select(one_of(common.cols))
  
  mat <- rbind(mat1a, mat2a)
  mat[mat == 0] <- NA
  library(zoo)
  mat.impute <- na.aggregate(mat, FUN = function(x) min(x)/2)
  logmat <- log2(mat.impute)

  # Prepare metadata in same order
  meta1 <- crc1 %>% select(Match_Caseset, Sex, Country.x, Cncr_Caco_Clrt, study, GRS)
  meta2 <- crc2 %>% select(Match_Caseset, Sex, Country.x, Cncr_Caco_Clrt, study, GRS)
  meta <- rbind(meta1, meta2)
  
  meta$Country <- as.factor(meta$Country.x)
  meta$study <- as.factor(meta$study)
  
  if(metabs == F) return(meta)
  
  # Get cutpoints for lowest and highest tertiles
  q <- quantile(meta$GRS, probs = c(0, 0.33, 0.67, 1)) # cutpoints are 3.39 and 3.81
  meta <- meta %>% mutate(GRS1 = ifelse(GRS > q[2] & GRS < q[3], NA, 1), GRSgroup= ifelse(GRS > 3.5, 1, 0))
  output <- data.frame(meta, logmat)
  
}  

# Data for associations between GRS and metabolites for 1. Controls and 2. Cases
controls <- prep.data(snps, group = 0)
cases    <- prep.data(snps, group = 1)
meta     <- prep.data(snps, metabs = F)

# Function to do logistic regression and get volcano plot data
volcano.data <- function(dat){
  
  df <- dat %>% filter(!is.na(GRS1))
  table(dat$GRSgroup)
  
  # Subset metadata (not necessary but makes model objects smaller) and concs matrix
  metaGRS <- select(df, Match_Caseset:GRSgroup)
  logmat1 <- select(df, -(Match_Caseset:GRSgroup)) %>% as.matrix
  
  # logistic regression scores low (0) vs high (0), apply across metabo matrix
  glm.fa <- function(x) glm(GRSgroup ~ x + Sex + study + Country, data = metaGRS, family = "binomial")
  multifit <- apply(logmat1, 2, glm.fa)
  
  # Extract FDR adjusted p-values and put in data frame with compound names
  library(broom)
  p <- map_df(multifit, tidy) %>% filter(term == "x") #%>% select(p.value) %>% pull
  p.adj <- p.adjust(p$p.value, method = "fdr")
  
  # Read in Biocrates compound metadata
  cmpd.meta <- read.csv("Biocrates_cmpd_metadata.csv")
  alldf <- data_frame(p.adj, Compound = colnames(logmat1)) %>% left_join(cmpd.meta, by = "Compound")
  
  # calculation of fold changes for volcano plot (need to use imputed matrix)
  df    <- data.frame(grp = metaGRS$GRSgroup, logmat1)
  means <- aggregate(. ~ grp, data = df, mean) %>% t
  
  # Get fold change of high score over low score
  alldf$meanfc   <- (means[, 2] - means[, 1])[-1]
  
  # If meanfc  > 0, conc in high GRS > conc in low GRS
  output <- alldf %>% mutate(direction = ifelse(meanfc > 0, "high", "low")) %>%
    separate(Compound, into = c("subclass", "rest"), sep = "_", extra = "merge", remove = F)
}

ctrls <- volcano.data(controls)
cases <- volcano.data(controls)

# Plot results for metabolites vertically
library(ggplot2)
ggplot(ctrls, aes(y = reorder(displayname, p.adj), x = log10(p.adj), shape = direction, colour = direction)) + 
  theme_minimal(base_size = 10) +
  geom_point(show.legend = F) + 
  geom_vline(xintercept = log10(0.05), linetype = "dashed") +
  xlab("-log10(FDR-adjusted p-value)") + ylab("") +
  facet_grid(subclass ~ ., scales = "free_y", space = "free_y", switch= "x") +
  theme(strip.text.y = element_blank())
#axis.text.x(angle = 0, size=7, hjust = 0.95, vjust = 0.5)) +
#ggtitle("Metabolite associations with WCRF score (cal)") +
#ggsave("WCRF score associations FAs.svg")

# or horizontally
library(ggplot2)
ggplot(alldf, aes(x = reorder(displayname, p.adj), y = log10(p.adj), shape = direction, colour = direction)) + 
  theme_minimal(base_size = 10) +
  geom_point(show.legend = F) + 
  geom_hline(yintercept = log10(0.05), linetype = "dashed") +
  ylab("-log10(FDR-adjusted p-value)") + xlab("") +
  facet_grid(. ~ subclass, scales = "free_x", space = "free_x", switch= "y") +
  theme(strip.text.y = element_blank(),
        axis.text.x = element_text(angle = 90, size=7, hjust = 0.95, vjust = 0.5))


# Same as above for fatty acids (formerly GRS_fatty acids.R)
volcano.FAs <- function(CaCo = 0, metabs = T) {
  
  library(tidyverse)
  fa.scores <- readRDS("FA_WCRF_scores.rds")
  # CRC case-control fatty acids dataset (from Elom) Convert categorical co-variates to factors
  var.list <- c("L_School", "Smoke_Stat")
  CRCfa1 <- read_dta("Database_Fatty acids.dta") %>% mutate_at(vars(var.list), as.factor)
  
  # Get CC obs which have SNP data
  intersect(snps$Idepic, CRCfa1$Idepic) %>% length
  crcsnps1 <- inner_join(CRCfa1, snps, by="Idepic")
  table(crcsnps1$Cncr_Caco_Clrt)
  # 422 controls, 416 cases in small dataset
  
  # For either case-control model or associations with metabolites
  crcsnps1 <- if(metabs == F) {
    crcsnps1 %>% group_by(Match_Caseset) %>% filter(n() == 2) %>% ungroup()
  } else {
    crcsnps1 %>% filter(Cncr_Caco_Clrt == CaCo)
  }
  
  # subset Biocrates data from the Fatty acids dataset and arrange by the common columns. Exclude FAs with many missings
  concs <- crcsnps1 %>% select(P14_0 : PCLA_9t_11c, -P24_0, -P20_0) %>% as.matrix
  concs[concs == 0] <- NA
  library(zoo)
  mat.impute <- na.aggregate(concs, FUN = function(x) min(x)/2)
  logmat <- log2(mat.impute)
  
  # Prepare metadata in same order
  meta <- crcsnps1 %>% select(Match_Caseset, Country.x, Cncr_Caco_Clrt, GRS)
  meta$Country <- as.factor(meta$Country.x)
  
  if(metabs == F) return(meta)
  
  # Get cutpoints for lowest and highest tertiles
  q <- quantile(meta$GRS, probs = c(0, 0.33, 0.67, 1)) # cutpoints are 3.39 and 3.81
  
  meta <- meta %>% mutate(GRS1 = ifelse(GRS > q[2] & GRS < q[3], NA, 1), GRSgroup= ifelse(GRS > 3.5, 1, 0))
  df <- data.frame(meta, logmat) %>% filter(!is.na(GRS1))
  table(df$GRSgroup)
  metaGRS <- select(df, Match_Caseset:GRSgroup)
  
  logmat1 <- select(df, -(Match_Caseset:GRSgroup)) %>% as.matrix
  
  #logistic regression scores low (0) vs high (0), apply across metabo matrix
  glm.fa <- function(x) glm(GRSgroup ~ x + Country, data = metaGRS, family = "binomial")
  multifit <- apply(logmat1, 2, glm.fa)
  
  # Extract FDR adjusted p-values and put in data frame with compound names
  library(broom)
  p <- map_df(multifit, tidy) %>% filter(term == "x") #%>% select(p.value) %>% pull
  p.adj <- p.adjust(p$p.value, method = "fdr")
  
  cmpd.meta <- read.csv("FA_compound_data.csv")
  alldf <- data_frame(p.adj, Compound = colnames(concs)) %>% inner_join(cmpd.meta, by = "Compound")
  
  # calculation of fold changes for volcano plot (need to use imputed matrix)
  df       <- data.frame(grp = metaGRS$GRSgroup, logmat1)
  means    <- aggregate(. ~ grp, data = df, mean) %>% t
  
  # Get fold change of high score over low score
  alldf$meanfc   <- (means[, 2] - means[, 1])[-1]
  alldf <- alldf %>% mutate(direction = ifelse(meanfc > 0, "high", "low")) %>%
    separate(Compound, into = c("subclass", "rest"), sep = "_", extra = "merge", remove = F)
  
}

# Data for association between CC status and GRS (return CC status and GRS only)
meta <- GRS.FA(metabs = F)

# Data for associations between GRS and metabolites for 1. Controls and 2. Cases
controls <- volcano.FAs(CaCo = 0)
cases    <- volcano.FAs(CaCo = 1)


# Plot results for metabolites vertically
library(ggplot2)
ggplot(controls, aes(y = reorder(displayname2, p.adj), x = log10(p.adj), shape = direction, colour = direction)) + 
  theme_minimal(base_size = 10) +
  geom_point(show.legend = F) + 
  geom_vline(xintercept = log10(0.05), linetype = "dashed") +
  xlab("-log10(FDR-adjusted p-value)") + ylab("") +
  facet_grid(sumgroup ~ ., scales = "free_y", space = "free_y", switch= "x") +
  theme(strip.text.y = element_blank())
#axis.text.x(angle = 0, size=7, hjust = 0.95, vjust = 0.5)) +
#ggtitle("Metabolite associations with WCRF score (cal)") +
#ggsave("WCRF score associations FAs.svg")

# or horizontally
library(ggplot2)
ggplot(cases, aes(x = reorder(displayname2, p.adj), y = p.adj, shape = direction, colour = direction)) + 
  theme_minimal(base_size = 10) +
  geom_point(show.legend = F) + 
  geom_hline(yintercept = 0.05, linetype = "dashed") +
  ylab("-log10(FDR-adjusted p-value)") + xlab("") +
  facet_grid(. ~ sumgroup, scales = "free_x", space = "free_x", switch= "y") +
  theme(strip.text.y = element_blank(),
        axis.text.x = element_text(angle = 90, size=7, hjust = 0.95, vjust = 0.5))


# CLR modelling case-control status from GRS (no metabolites)
par(mfrow= c(1,2))
boxplot(GRS ~ Cncr_Caco_Clrt, data = meta, ylab= "Genetic risk score", col = "limegreen")

# Do not adjust for country and study because they are already matched on this...
library(survival)
fit <- clogit(Cncr_Caco_Clrt ~ GRS + strata(Match_Caseset), data = meta)
summary(fit)

# unmatched
fit2 <- glm(Cncr_Caco_Clrt ~ GRS + Country + study + Sex, data = meta)
summary(fit2)
st <- stargazer(fit, type = "text", ci = T, apply.coef = exp)


# Get signature of PRS as a PLS object (formerly CRC_GRS_signature.R)
signatureGRS <- function() {
  # First adjust matrix for study, centre, sex, batch, BMI
  logmat <- select(controls, -(Match_Caseset:GRSgroup)) %>% as.matrix
  adj    <- function(x) residuals(lm(x ~ study + Country + Sex, data = controls))
  adjmat <- apply(logmat, 2, adj)
  
  # PLS model for metabolic signature of GRS
  # Subset GRS and bind to residuals-adjusted matrix
  grs <- data_frame(score = controls$GRS)
  plsdata <- cbind(grs, adjmat)
  
  library(pls)
  mod <- plsr(score ~ ., data = plsdata, validation = "CV")
  #summary(mod)
  
  # Find the number of dimensions with lowest cross validation error
  cv <- RMSEP(mod)
  plot(RMSEP(mod), legendpos = "topright")
  
  # Calculate optimal number of dimensions and rerun model
  best.dims <- which.min(cv$val[estimate = "adjCV", , ]) - 1
  mod <- plsr(score ~ ., data = plsdata, ncomp = best.dims)
  
  return(mod)
  
  # explained variances
  explvar(mod)
  
  # Coefficients and variable importance
  coefficients <- coef(mod)
  sum.coef <- sum(sapply(coefficients, abs))
  coefficients <- coefficients * 100 / sum.coef
  coefficients <- sort(coefficients[, 1 , 1])
  # plot(coefficients)
  
  # Get top and bottom deciles of compound coefficients
  df <- data.frame(as.list(coefficients)) %>% gather(Cmpd, VIP)
  
  qs <- quantile(coefficients, probs = seq(0, 1, 0.05))
  df1 <- df[df$VIP > qs[18], ]
  df2 <- df[df$VIP < qs[4], ]
  
  # if(modonly == T) return(list(df1, df2))
  # Vector of colours for plot points
  vec <- ifelse(df$VIP > qs[18] | df$VIP < qs[4], "black", "grey")
  
  # Now plot data, adding text
  plot(coefficients, pch = 17, col=vec, xlab = "", ylab = "Variable Importance",
       main = paste(nrow(plsdata), "subjects, optimal dimensions =", best.dims))
  text(nrow(df) - nrow(df1):1, df1$VIP, df1$Cmpd, pos=2, cex = 0.6)
  text(1:nrow(df2), df2$VIP, df2$Cmpd, pos=4, cex=0.6)
  abline(a=0, b=0, lty = "dotted")
}
mod <- signatureGRS()

