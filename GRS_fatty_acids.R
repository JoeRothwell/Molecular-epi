# Association between genetic risk score and fatty acids

GRS.FA <- function(CaCo = 0, metabs = T) {
  
  library(haven)
  library(tidyverse)
  
  # Read in snps data with GRS calculated from Todd's STATA file
  snps <- read_dta("clrt_gwas_gecco_snps_GRS.dta")
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
controls <- GRS.FA(CaCo = 0)
cases    <- GRS.FA(CaCo = 1)


# CLR modelling case-control status from GRS
par(mfrow= c(1,2))
boxplot(GRS ~ Cncr_Caco_Clrt, data = meta, ylab= "Genetic risk score", col = "limegreen")

# Do not adjust for country and study because they are already matched on this...
library(survival)
fit <- clogit(Cncr_Caco_Clrt ~ GRS + strata(Match_Caseset), data = meta)
summary(fit)

# unmatched
fit2 <- glm(Cncr_Caco_Clrt ~ GRS + Country + study, data = meta)
summary(fit2)

st <- stargazer(fit, type = "text", ci = T, apply.coef = exp)



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