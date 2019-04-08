# LIFEPATH data exploratory

library(tidyverse)
library(haven)

path <- "Y:/RepPerso/Fabienne WILM/02_Demandes_ponctuelles/10_LIFEPATH/TABLES"
list.files(path)
dat1 <- read_csv("D01_20171031_LIFEPATH.csv")
dat2 <- read_csv("D01_20161018_LIFEPATH.csv")
dat3 <- read_csv("D01_20150917_LIFEPATH.csv")

write_dta(dat1, "D01_20171031_LIFEPATH.dta")
write_dta(dat2, "D01_20161018_LIFEPATH.dta")
write_dta(dat3, "D01_20150917_LIFEPATH.dta")


