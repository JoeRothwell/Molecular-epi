# Extracts code names and label names and converts to a dataframe

  library(tidyverse)
  library(haven)
  
# Biocrates  
ctrl <- read_dta("obes_metabo.dta") %>%
        select(matches("Acylcarn_|Aminoacid_|Biogenic_|Glyceroph_|Sphingo_|Sugars_"), 
                              -starts_with("Outdq"))
labels <- sapply(controls, attr, "label")
tibble(name = names(labels), label = labels)


# Fatty acids
CRCfa1 <- read_dta("Database_Fatty acids.dta")
  
CRCfa <- CRCfa1 %>% select(P14_0 : PCLA_9t_11c) 
df <- tibble(Compound = colnames(CRCfa))



