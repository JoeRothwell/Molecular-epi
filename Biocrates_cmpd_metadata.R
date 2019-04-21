# Extracts code names and label names and converts to a dataframe

getnames <- function() {
  library(tidyverse)
  library(haven)
  ctrl <- read_dta("obes_metabo.dta") %>%
          select(matches("Acylcarn_|Aminoacid_|Biogenic_|Glyceroph_|Sphingo_|Sugars_"), 
                              -starts_with("Outdq"))
  labels <- sapply(controls, attr, "label")
  tibble(name = names(labels), label = labels)
}


