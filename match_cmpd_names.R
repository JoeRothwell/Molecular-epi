# Need to first match names to standard names in CRC study
cmpd_meta <- read_csv("Biocrates_cmpd_metadata.csv")

# Make data frame of Biocrates compound names
#allnames <- dat %>% select(lyso_pc_a_c16_0:c9)
matnames <- data.frame(cmpd.low = colnames(mat2), ord = 1:length(colnames(mat2)))
#df <- data.frame(Compound.low = colnames(allnames), ord = 1:length(colnames(allnames)))

# Make new names that match those in adenoma dataset
cmpd_meta2 <- cmpd_meta %>% separate(Compound, into = c("Compound.cl", "Compound.str"),
  remove = F, extra = "merge", fill = "right") %>%
  mutate(cmpd.low = str_to_lower(Compound.str)) %>% 
  mutate(cmpd.low = str_replace(cmpd.low, "lyso", "lyso_")) %>%
  mutate(cmpd.low = str_replace(cmpd.low, "c4_oh_", "c4_oh"))

library(fuzzyjoin)
df1 <- stringdist_right_join(cmpd_meta2, matnames, by = "cmpd.low", max_dist = 0.5)
