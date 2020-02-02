# Baseline characteristics for the 2 CRC metabolomics studies
# 8 unpaired samples have been removed for CRC1

source("CRC_prep_data.R")

crc1a <- crc1 %>% group_by(Match_Caseset) %>% filter(n() == 2) %>%
  select(Sex, Age_Blood, Height_C, Weight_C, Bmi_C, Qe_Energy, Country, Pa_Mets, Smoke_Stat, Qe_Alc,
         Wcrf_C_Cal, Cncr_Caco_Clrt) %>% mutate(Study = "CRC1")

crc2a <- crc2 %>% #left_join(wcrf, by = "Idepic") %>%
  select(Sex, Age_Blood, Height_C, Weight_C, Bmi_C, Qe_Energy, Country, Pa_Mets, Smoke_Stat, Qe_Alc,
         Wcrf_C_Cal, Cncr_Caco_Clrt) %>% mutate(Study = "CRC2")
crcmetab <- bind_rows(crc1a, crc2a)

library(qwraps2)
options(qwraps2_markup = "markdown")

# Generate table automatically
#crcsum <- bind_rows(crc1a, crc2a) %>% qsummary(., numeric_summaries = list("Mean (SD)" = "~ mean_sd(%s)"),
#           n_perc_args = list(digits = 1, show_symbol = TRUE))

# Manually specified
crc_sum <-
  list("Total subjects"   =
         list("N"   =   ~ n()),
       "Sex" = 
         list("Female"   =  ~ n_perc0(Sex == 1, digits = 1),
              "Male"     =  ~ n_perc0(Sex == 2, digits = 1)),
       "Age at blood collection (years)" = 
         list("Mean"     =  ~ mean_sd(Age_Blood)),
       "Height (cm)" =
         list("Mean (SD)" = ~ mean_sd(Height_C)),
       "Weight (kg)" =
         list("Mean (SD)" = ~ mean_sd(Weight_C)),
       "BMI (kg/m2)" =
         list("Mean (SD)" = ~ mean_sd(Bmi_C)),
       "Total energy intake (kCal)" =
         list("Mean (SD)" = ~ mean_sd(Qe_Energy, na_rm = T, show_n = "never")),
       "Country" = 
         list("France"      =  ~ n_perc0(Country == 1, digits = 1),
              "Italy"       =  ~ n_perc0(Country == 2, digits = 1),
              "Spain"       =  ~ n_perc0(Country == 3, digits = 1),
              "United Kingdom"       =  ~ n_perc0(Country == 4, digits = 1),
              "Netherlands" =  ~ n_perc0(Country == 5, digits = 1),
              "Greece"          =  ~ n_perc0(Country == 6, digits = 1),
              "Germany"      =  ~ n_perc0(Country == 7, digits = 1),
              "Denmark"      = ~ n_perc0(Country == 9, digits = 1)),
       "Physical activity" =
         list("Mean (SD)" = ~ mean_sd(Pa_Mets, na_rm = T, show_n = "never")),
       "Alcohol intake (g/day)" =
         list("Mean (SD)" = ~ mean_sd(Qe_Alc, na_rm = T, show_n = "never")),
       "Smoking status" =
         list("Non smoker" =  ~ n_perc0(Smoke_Stat == 1, digits = 1),
              "Never smoker"   =  ~ n_perc0(Smoke_Stat == 2, digits = 1),
              "Smoker"       =  ~ n_perc0(Smoke_Stat == 3, digits = 1)),
       #"Unknown"      =  ~ n_perc0(Smoke_Stat == 4, digits = 1)),
       "WCRF score" = 
         list("Mean (SD)" = ~ mean_sd(Wcrf_C_Cal, na_rm = T, show_n = "never"))
  )
st1 <- summary_table(group_by(crc1a, Cncr_Caco_Clrt), crc_sum)
st2 <- summary_table(group_by(crc2a, Cncr_Caco_Clrt), crc_sum)
both <- cbind(st1, st2)
print(both, cnames = c("CRC1 cases", "CRC1 controls", "CRC2 cases", "CRC2 controls"))
# Copy and paste output into an Rmarkdown file and render to word/pdf etc


# Tests for p-values for table -----------------------------
# McNemar and Wilcoxon signed rank test?

# Case/control models
ll1 <- list(
  chisq.test(crc1a$Cncr_Caco_Clrt, crc1a$Height_C),
  chisq.test(crc1a$Cncr_Caco_Clrt, crc1a$Weight_C),
  chisq.test(crc1a$Cncr_Caco_Clrt, crc1a$Age_Blood),
  chisq.test(crc1a$Cncr_Caco_Clrt, crc1a$Bmi_C),
  chisq.test(crc1a$Cncr_Caco_Clrt, crc1a$Qe_Energy),
  chisq.test(crc1a$Cncr_Caco_Clrt, crc1a$Qe_Alc),
  chisq.test(crc1a$Cncr_Caco_Clrt, crc1a$Smoke_Stat)
)

# Case/control models
ll2 <- list(
  chisq.test(crc2a$Cncr_Caco_Clrt, crc2a$Height_C),
  chisq.test(crc2a$Cncr_Caco_Clrt, crc2a$Weight_C),
  chisq.test(crc2a$Cncr_Caco_Clrt, crc2a$Age_Blood),
  chisq.test(crc2a$Cncr_Caco_Clrt, crc2a$Bmi_C),
  chisq.test(crc2a$Cncr_Caco_Clrt, crc2a$Qe_Energy),
  chisq.test(crc2a$Cncr_Caco_Clrt, crc2a$Qe_Alc),
  chisq.test(crc2a$Cncr_Caco_Clrt, crc2a$Smoke_Stat)
)


t.test(meta$DURTHSDIAG ~ meta$CT)$p.value

# Pre/post menopausal models
l  <- list(
  chisq.test(meta$MENOPAUSE, meta$DIAGSAMPLINGCat1),
  chisq.test(meta$MENOPAUSE, meta$BEHAVIOUR),
  fisher.test(meta$MENOPAUSE, meta$SUBTYPE),
  chisq.test(meta$MENOPAUSE, meta$ER),
  chisq.test(meta$MENOPAUSE, meta$PR),
  chisq.test(meta$MENOPAUSE, meta$GRADE),
  fisher.test(meta$MENOPAUSE, meta$STADE)
)


t.test(meta$DURTHSDIAG ~ meta$CT)$names

# Extract data from models
library(purrr)
library(broom)
map_df(ll, tidy)
map_df(l, tidy)

# Alcohol intake