library(qwraps2)
library(haven)
library(tidyverse)

# Descriptives for CRC case-control studies. 
# OUtput markdown to be copied and pasted into tables_crc_ms.Rmd and knitted

crc1a <- crc1 %>% 
  select(Sex, Age_Blood, Height_C, Weight_C, Bmi_C, Qe_Energy, Country, Pa_Mets, Smoke_Stat, 
         Wcrf_C_Cal, Cncr_Caco_Clrt) %>% mutate(Study = "CRC1")
crc2a <- crc2 %>% 
  select(Sex, Age_Blood, Height_C, Weight_C, Bmi_C, Qe_Energy, Country, Pa_Mets, Smoke_Stat, 
         Wcrf_C_Cal, Cncr_Caco_Clrt) %>% mutate(Study = "CRC2")
crcmetab <- bind_rows(crc1a, crc2a)

crc1_summaries <-
  list("Sex" = 
         list("Female"   =  ~ n_perc0(Sex == 1),
              "Male"     =  ~ n_perc0(Sex == 2)),
       "Age at blood collection" = 
         list("mean"     =  ~ mean_sd(Age_Blood)),
       "Height" =
         list("mean (sd)" = ~ mean_sd(Height_C)),
       "Weight" =
         list("mean (sd)" = ~ mean_sd(Weight_C)),
       "BMI" =
         list("mean (sd)" = ~ mean_sd(Bmi_C)),
       "Total energy intake" =
         list("mean (sd)" = ~ mean_sd(Qe_Energy)),
       "Country" = 
         list("France"      =  ~ n_perc0(Country == 1),
              "Italy"       =  ~ n_perc0(Country == 2),
              "Germany"     =  ~ n_perc0(Country == 3),
              "Spain"       =  ~ n_perc0(Country == 4),
              "Netherlands" =  ~ n_perc0(Country == 5),
              "UK"          =  ~ n_perc0(Country == 6),
              "Greece"      =  ~ n_perc0(Country == 7)),
       "Physical activity" =
         list("mean (sd)" = ~ mean_sd(Pa_Mets)),
       "Smoking status" =
         list("Never smoker" =  ~ n_perc0(Smoke_Stat == 1),
              "Non-smoker"   =  ~ n_perc0(Smoke_Stat == 2),
              "Smoker"       =  ~ n_perc0(Smoke_Stat == 3),
              "Unknown"      =  ~ n_perc0(Smoke_Stat == 4)),
       "WCRF score" = 
         list("mean (sd)" = ~ mean_sd(Wcrf_C_Cal))
  )

st1 <- summary_table(group_by(crc1a, Cncr_Caco_Clrt), crc1_summaries)
st2 <- summary_table(group_by(crc2a, Cncr_Caco_Clrt), crc1_summaries)
both <- cbind(st1, st2)
both

# Descriptives for all EPIC

library(haven)
library(tidyverse)
wcrf <- read_dta("wcrf_score.dta")
fullepic <- read_dta("D:/full epic.dta")
fullepic <- fullepic %>% left_join(wcrf, by = "Idepic")

fullepic_summaries <-
  list("Sex" = 
         list("Female"   =  ~ n_perc0(Sex == 1),
              "Male"     =  ~ n_perc0(Sex == 2)),
       "Height" =
         list("mean (sd)" = ~ mean_sd(Height_C)),
       "Weight" =
         list("mean (sd)" = ~ mean_sd(Weight_C)),
       "BMI" =
         list("mean (sd)" = ~ mean_sd(Bmi_C)),
       "Total energy intake" =
         list("mean (sd)" = ~ mean_sd(QE_ENERGY)),
       "Country" = 
         list("France"     =  ~ n_perc0(Country.x == 1),
              "Italy"      =  ~ n_perc0(Country.x == 2),
              "Germany"    =  ~ n_perc0(Country.x == 3),
              "Spain"      =  ~ n_perc0(Country.x == 4),
              "Netherlands" =  ~ n_perc0(Country.x == 5),
              "UK"         =  ~ n_perc0(Country.x == 6),
              "Greece"     =  ~ n_perc0(Country.x == 7)),
       "Physical activity" =
         list("mean (sd)" = ~ mean_sd(Pa_Mets)),
       "Smoking status" =
         list("Never smoker" =  ~ n_perc0(Smoke_Stat == 1),
              "Non-smoker"   =  ~ n_perc0(Smoke_Stat == 2),
              "Smoker"       =  ~ n_perc0(Smoke_Stat == 3),
              "Unknown"      =  ~ n_perc0(Smoke_Stat == 4)),
       "WCRF score" = 
         list("mean (sd)" = ~ mean_sd(Wcrf_C_Cal))
  )

st3 <- summary_table(fullepic, fullepic_summaries)




# From https://www.r-bloggers.com/baseline-characteristics-tables-with-qwraps2/
# Example on mtcars data. First set type to markdown

options(qwraps2_markup = "markdown")

mtcar_summaries <-
  list("Miles Per Gallon" =
         list("min:"         = ~ min(mpg),
              "mean (sd)"    = ~ qwraps2::mean_sd(mpg, denote_sd = "paren"),
              "median (IQR)" = ~ qwraps2::median_iqr(mpg),
              "max:"         = ~ max(mpg)),
       "Cylinders:" = 
         list("mean"             = ~ mean(cyl),
              "mean (formatted)" = ~ qwraps2::frmt(mean(cyl)),
              "4 cyl, n (%)"     = ~ qwraps2::n_perc0(cyl == 4),
              "6 cyl, n (%)"     = ~ qwraps2::n_perc0(cyl == 6),
              "8 cyl, n (%)"     = ~ qwraps2::n_perc0(cyl == 8)),
       "Weight" =
         list("Range" = ~ paste(range(wt), collapse = ", "))
  )

st <- summary_table(mtcars, mtcar_summaries)
# Now copy and paste markdown output into an rmd file
