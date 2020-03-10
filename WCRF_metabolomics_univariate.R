# Performs logistic regression between low and high WCRF score groups for EPIC controls dataset
# Volcano plot of p-value vs fold change
# To get data prep from CRC_prep.data.R
source("CRC_prep_data.R")

library(tidyverse)

# Model scores from Biocrates data. Code 0 and 1 as categories 1 and 2 vs 4 and 5
ob <- ctrl %>% filter(Wcrf_C_Cal != 3) %>% mutate(score_cat = ifelse(Wcrf_C_Cal %in% 1:2, 0, 1))

# Subset metabolite matrix, replace zero with NA

controls <- ob %>% 
  select(matches("Acylcarn_|Aminoacid_|Biogenic_|Glyceroph_|Sphingo_|Sugars_"), -starts_with("Outdq"))

zerocols <- apply(controls, 2, function(x) sum(x, na.rm = T)) != 0
controls <- controls[, zerocols]
colnames(controls) %>% length # 159 variables

concs <- as.matrix(controls)
concs[concs == 0] <- NA
library(zoo)
metabo <- na.aggregate(concs, FUN = function(x) min(x)/2)
logconcs <- log2(metabo)

# logistic regression on categories 1 and 2 vs 4 and 5

glm.ob <- function(x) glm(score_cat ~ x + Sex + Center + Smoke_Intensity + Study, data = ob, family = "binomial") 

# apply function across metabolite matrix
multifit <- apply(logconcs, 2, glm.ob)

library(broom)
p <- map_df(multifit, tidy) %>% filter(term == "x") %>% select(p.value) %>% pull
p.adj <- p.adjust(p, method = "BH")

# calculation of fold changes for volcano plot
df       <- data.frame(ob$score_cat, metabo)
means    <- aggregate(. ~ ob.score_cat, data = df, mean) %>% t
# Get fold change of high score over low score
meanfc   <- means[, 2]/means[, 1]

# make final data frame adding extra columns for point labelling
cmpd.meta <- read.csv("Biocrates_cmpd_metadata.csv")
cal <- 
  data_frame(Compound = colnames(metabo), Fold_change = meanfc[-1], p.value = p.adj) %>%
  separate(Compound, into = c("subclass", "rest"), sep = "_", extra = "merge", remove = F) %>%
  mutate(
         # variables for Manhattan
         subclass1   = ifelse(str_detect(Compound, "Lysopc"), "LysoPC", subclass),
         direction   = ifelse(Fold_change > 1, "Increased with score", "Decreased with score"),
         Association = ifelse(p.value < 0.001, direction, "Not significant"),
         
         # For Volcano, specify points to label
         tolabel     = ifelse(p.adj < 0.001 & abs(1 - Fold_change) > 0.05, Compound, NA),
         tocolour    = ifelse(p.adj < 0.001 & abs(1 - Fold_change) > 0.05, direction, "col3")

         ) %>% 
  left_join(cmpd.meta, by = "Compound")

volcano.fa <- function() {
  
  #library(haven)
  library(tidyverse)
  
  # Read in fatty acids data and WCRF scores
  #fa <- read_sas("controls_fas1.sas7bdat") 
 # wcrf <- read_dta("Wcrf_Score.dta") %>% select(Idepic, Wcrf_C_Cal)
  
  # Join scores to fatty acids data
  #fa.scores <- fa %>% left_join(wcrf, by = "Idepic")
  #saveRDS(fa.scores, "FA_WCRF_scores1.rds")
  
  # Note: FA_WCRF_scores1 is the new file from Carine with technical covariates (16 Nov 2018)
  fa.scores <- readRDS("FA_WCRF_scores1.rds")
  
  fa <- fa.scores %>% filter(Wcrf_C_Cal != 3) %>% mutate(score_cat = ifelse(Wcrf_C_Cal %in% 1:2, 0, 1))
  # check numbers of low and high
  fa %>% group_by(score_cat) %>% summarise(subjects = n())
  
  # Prepare metabolite matrix. Set zeros to NA then impute with half min value
  concs <- fa %>% select(P14_0 : PCLA_10t_12c) %>% as.matrix
  concs[concs == 0] <- NA
  
  library(zoo)
  concs <- na.aggregate(concs, FUN = function(x) min(x)/2)
  logconcs <- log2(concs)
  
  #library(corrplot)
  #corrplot(cor(logconcs), method = "square", tl.col = "black")
  
  #logistic regression scores 1 and 2 vs 4 and 5, apply across metabo matrix
  glm.fa <- function(x) glm(score_cat ~ x + Center, data = fa, family = "binomial")
  multifit <- apply(logconcs, 2, glm.fa)
  
  library(broom)
  p <- map_df(multifit, tidy) %>% filter(term == "x") #%>% select(p.value) %>% pull
  p.adj <- data_frame(p.adj = p.adjust(p$p.value, method = "BH"))
  
  # calculation of fold changes for volcano plot. Make sure to use non-log data
  cats  <- data.frame(scorecat = fa$score_cat, concs)
  
  cmpd.meta <- read.csv("FA_compound_data.csv")
  
  # Data frame for results. FAs are roughly categorised by string count. Use factor_key to keep order
  output <- cats %>% group_by(scorecat) %>% summarise_all(funs(mean)) %>%
    gather(Compound, conc, -scorecat, factor_key = T) %>%
    spread(scorecat, conc) %>%
    bind_cols(p, p.adj) %>%
    mutate(meanfc = `1`/`0`, associated = ifelse(meanfc > 1, "incr_score", "decr_score"),
           Cmpd.length = as.factor(str_count(Compound)) ) %>%
    inner_join(cmpd.meta, by = "Compound")
  
}

# Biocrates metabolites
cal <- volcano(Cal = T, adj = F)
raw <- volcano(Cal = F)

# Fatty acids
alldf <- volcano.fa()

# Bind datasets, exclude outlier serotonin, add four new columns
all <- bind_rows(list("Raw" = raw, "Cal" = cal), .id = "id") %>% filter(rest != "Serotonin")

# saveRDS(all, "Wcrf_biocrates1")

# Manhattan ----
# For supplemental data


# Horizontal, by subclass, in order of p-value (cal or raw)
ggplot(cal, aes(x = reorder(displayname, p.value), y = -log10(p.value), shape = direction, colour = direction)) + 
  theme_minimal(base_size = 10) +
  geom_point(show.legend = F) + 
  geom_hline(yintercept = 3, linetype = "dashed") +
  ylab("log10(p-value)") +
  facet_grid(. ~ subclass1, scales = "free_x", space = "free_x", switch= "y") +
  theme(axis.title.x = element_blank(), strip.text.x = element_blank(),
        axis.text.x = element_text(angle = 90, size=7, hjust = 0.95, vjust = 0.5)) #+
  #ggtitle("Metabolite associations with WCRF score (cal)") +
  #ggsave("WCRF score associations for slide.svg")


# Horizontal
library(ggplot2)
ggplot(alldf, aes(x = reorder(displayname, -p.adj), y = -log10(p.adj), shape = associated, colour = associated)) + 
  theme_minimal(base_size = 10) +
  geom_point(show.legend = F) + 
  geom_hline(yintercept = 3, linetype = "dashed") +
  ylab("log10(FDR-adjusted p-value)") + xlab("") +
  facet_grid(. ~ sumgroup, scales = "free_x", space = "free_x", switch= "y") +
  theme(strip.text.x = element_blank(), 
        axis.text.x = element_text(angle = 90, size= 8, hjust = 0.95, vjust = 0.5)) #+
#ggsave("WCRF score associations FAs.svg")


# Old ---------------------


# Vertical, by subclass, in order of p-value (cal or raw)
ggplot(cal, aes(y = reorder(displayname, p.value), x = log10(p.value), shape = direction, colour = direction)) + 
  theme_minimal(base_size = 10) +
  geom_point() + geom_vline(xintercept = -3, linetype = "dashed") +
  xlab("log10(p-value)") +
  facet_grid(subclass1 ~ ., scales = "free_y", space = "free_y", switch= "x") +
  theme(axis.title.y = element_blank(),
        axis.text.y = element_text(size=6),
        legend.position = c(0.25, 0.4),
        legend.box.background = element_rect(colour="grey")) +
  #ggtitle("Metabolite associations with WCRF score (cal)")
  

  
# Horizontal, significant associations coloured
library(ggplot2) 

ggplot(raw, aes(x = reorder(Compound, subclass), y = -log10(p.value), colour = Association)) + 
#theme_bw() +
geom_point() + geom_hline(yintercept = 3, linetype = "dashed") +
scale_color_manual(values = c("blue", "red", "grey")) +
ylab("-log10(p-value) for association with WCRF score") + xlab("Compound") +
theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust=0.5, size = 7),
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
      legend.position = c(0.03,0.75),
      legend.justification = c(0, 0))  
  
# Volcano plot (facetted by raw/calibrated) ----
# Needs extra variables
library(ggplot2)
ggplot(cal, aes(x=Fold_change, y = -log10(p.value), colour = tocolour, shape = tocolour)) + 
  #geom_point(shape=1, colour = "dodgerblue", alpha = 0.7) + 
  geom_point(show.legend = F) +
  scale_colour_manual(values = c("limegreen", "red", "grey")) +
  scale_shape_manual(values = c(17,19,17)) +
  theme_bw() +
  geom_text(aes(label = tolabel), hjust = -0.05, vjust = 0, size = 2, colour = "black") +
  facet_grid(. ~ id) +
  xlab("Fold change (relative concentration in higher scoring subjects)") +
  ggtitle("Metabolite association with WCRF score", 
          subtitle = "Adjusted for sex, EPIC centre, study and smoking intensity")

# plot(cal$p.value)

raw %>% group_by(subclass) %>% 
  separate(Compound, into = c("subclass", "rest"), sep = "_", extra="merge", remove = F)

table(ctrl$Wcrf_C)
table(ctrl$Wcrf_C_Cal)

# Fatty acids

# Volcano plot
ggplot(alldf, aes(x=meanfc, y = -log10(p.adj))) + geom_point(show.legend = F) + 
  theme_minimal()
#geom_point(shape=1, colour = "dodgerblue", alpha = 0.7) + 
#xlab("Fold change (relative concentration in higher scoring subjects)")
#scale_colour_manual(values = c("grey", "limegreen", "red")) +
#scale_shape_manual(values = c(19,17,19)) +
#xlim(0.91, 1.2) +
#geom_text(aes(label = labname), hjust = -0.05, 
#vjust = 0, size = 2, colour = "black") 


# Vertical Manhattan
library(ggplot2)
ggplot(output, 
       aes(y = reorder(displayname, -p.adj), x = -log10(p.adj), shape = associated, colour = associated)) + 
  theme_minimal(base_size = 10) +
  geom_point(show.legend = F) + 
  geom_vline(xintercept = 3, linetype = "dashed") +
  xlab("-log10(FDR-adjusted p-value)") + ylab("") +
  facet_grid(sumgroup ~ ., scales = "free_y", space = "free_y", switch= "x") +
  theme(strip.text.y = element_blank())
#axis.text.x(angle = 0, size=7, hjust = 0.95, vjust = 0.5)) +
#ggsave("WCRF score associations FAs.svg")

# ----

# Case-control dataset
library(haven)
data <- read_sas("clrt_caco_metabo.sas7bdat") 
concs <- data %>% select(Acylcarn_C0 : Sugars_H1) %>% as.matrix

# Converted to dta
# write_dta(data, "clrt_caco_metabo.dta")

# subset only observations with Biocrates data
bioc <- apply(concs, 1, function(x) sum(!is.na(x)) > 0)
sum(bioc) # 988 observations with biocrates data
df <- data[bioc, ]

# Remove odd cases or controls
df1 <- df %>% group_by(Match_Caseset) %>% filter(n() == 2) %>% ungroup
write.csv(df1, file = "CRC biocrates data.csv")

# Serotonin analysis

boxplot(log10(Biogenic_Serotonin) ~ Cncr_Caco_Clrt, data = df1)
sum(!is.na(df1$Biogenic_Serotonin))
df1$Biogenic_Serotonin


