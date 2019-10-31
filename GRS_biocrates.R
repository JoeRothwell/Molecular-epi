# Test for overlap between GWAS data and other datasets
source("CRC_data_prep.R")
source("Metabolic_signatures.R")

# snps data with GRS calculated from Todd's STATA file
library(haven)
snps <- read_dta("clrt_gwas_gecco_snps_GRS.dta")

prepdata <- function(dat, group = 0, metabs = T) {

  library(haven)
  library(tidyverse)
  
  # Read in snps data with GRS calculated from Todd's STATA file
  dat <- read_dta("clrt_gwas_gecco_snps_GRS.dta")
  
  # Get small CRC obs that have SNP data
  # crcsmall <- readRDS("CRC_smallerCC.rds") 
  crcsmall <- crc1
  intersect(data$Idepic, crcsmall$Idepic) %>% length
  
  # Join small CRC to SNPs
  crc1 <- inner_join(crcsmall, dat, by="Idepic")
  # table(crc1$Cncr_Caco_Clrt)
  # 444 controls, 442 cases in small dataset
  
  # For either case-control model or associations with metabolites
  crc1 <- if(metabs == F) {
    crc1 %>% group_by(Match_Caseset) %>% filter(n() == 2) %>% mutate(study = "small") %>% ungroup()
  } else {
    crc1 %>% filter(Cncr_Caco_Clrt == group) %>% mutate(study = "small")
  }
    
  
  # Large CRC dataset
  crclarge <- read_csv("biocrates_p150.csv")
  intersect(dat$Idepic, crclarge$Idepic) %>% length
  
  crc2 <- inner_join(crclarge, dat, by="Idepic")
  # table(crc2$Cncr_Caco_Clrt)
  # 296 controls, 439 cases in dataset
  
  # For either case-control model or associations with metabolites
  crc2 <- if(metabs == F) {
    crc2 %>% group_by(Match_Caseset) %>% filter(n() == 2) %>% mutate(study = "large") %>% ungroup()
  } else {
    crc2 %>% filter(Cncr_Caco_Clrt == group) %>% mutate(study = "large")
  }
  
  # subset Biocrates data from the two datasets and arrange by the common columns
  # mat1 <- crc1 %>% select(Acylcarn_C0 : Sugars_H1, -Batch_MetBio)
  mat1 <- crc1 %>% select(matches("Acylcarn_|Aminoacid_|Biogenic_|Glyceroph_|Sphingo_|Sugars_"), 
                 -starts_with("Outdq"))
  # mat2 <- crc2 %>% select(Acylcarn_C0 : Sugars_H1, -starts_with("Outdq_"))
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
controls <- prepdata(snps, group = 0)
cases    <- prepdata(snps, group = 1)
meta     <- prepdata(snps, metabs = F)


# Univariate ----
# logistic regression method
df <- controls %>% filter(!is.na(GRS1))
table(controls$GRSgroup)

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
alldf <- alldf %>% mutate(direction = ifelse(meanfc > 0, "high", "low")) %>%
  separate(Compound, into = c("subclass", "rest"), sep = "_", extra = "merge", remove = F)

# Plot results for metabolites vertically
library(ggplot2)
ggplot(alldf, aes(y = reorder(displayname, p.adj), x = log10(p.adj), shape = direction, colour = direction)) + 
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

# Model case-control status from GRS metabolite signature







