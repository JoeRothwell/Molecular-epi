# Baseline characteristics for the 2 CRC metabolomics studies
# 8 unpaired samples have been removed for CRC1
source("CRC_prep_data.R")
library(tidyverse)

#Get predicted WCRF score
load("predscore_df_subsite1.Rdata")
crc.ph0 <- crc.ph %>% select(Idepic, comp1)
crc.ph1 <- crc3.ph %>% select(Idepic, comp2)

# Generate time to diagnosis and tumour site variables for both studies
crc1a <- crc1 %>% group_by(Match_Caseset) %>% filter(n() == 2) %>% ungroup() %>%
  select(Idepic, Match_Caseset, Sex, location, Age_Blood, Tfollowup, Height_C, Bmi_C, Waist_C, Qe_Energy, Bdg1clrt,
         Country, Pa_Mets, Smoke_Stat, Qe_Alc, Wcrf_C_Cal, Cncr_Caco_Clrt) %>% 
  mutate(Study = "CRC1", Bdg1clrt_hist = ifelse(Bdg1clrt %in% c(53,54,55,56,60,70), 1, 0)) %>%
  left_join(crc.ph0, by = "Idepic") %>%
  left_join(crc.ph1, by = "Idepic")

crc2a <- crc2 %>%
  select(Idepic, Match_Caseset, Sex, location, Age_Blood, Tfollowup, Height_C, Bmi_C, Waist_C, Qe_Energy, Bdg1clrt,
         Country, Pa_Mets, Smoke_Stat, Qe_Alc, Wcrf_C_Cal, Cncr_Caco_Clrt) %>% 
  mutate(Study = "CRC2", Bdg1clrt_hist = ifelse(Bdg1clrt %in% c(53,54,55,56,60,70), 1, 0)) %>%
  left_join(crc.ph0, by = "Idepic")

# Merge for revised submission and add predicted WCRF score
crc.sig <- bind_rows(crc1a, crc2a)
crc.sig$Tfollowup <- as.numeric(crc.sig$Tfollowup)

library(qwraps2)
options(qwraps2_markup = "markdown")

# Manually specified summary
crc.sum <-
  list("Sex" = 
         list("Male"       =  ~ n_perc0(Sex == 1, digits = 1),
              "Female"     =  ~ n_perc0(Sex == 2, digits = 1)),
       "Age at blood collection (years)" = 
         list("Mean"     =  ~ mean_sd(Age_Blood, digits = 1)),
       
       "Follow-up time to diagnosis (years)" =
         list("Mean"     =  ~ mean_sd(Tfollowup, digits = 1)),
       
       "Country" = 
         list("France"         = ~ n_perc0(Country == 1, digits = 1),
              "Italy"          = ~ n_perc0(Country == 2, digits = 1),
              "Spain"          = ~ n_perc0(Country == 3, digits = 1),
              "United Kingdom" = ~ n_perc0(Country == 4, digits = 1),
              "Netherlands"   = ~ n_perc0(Country == 5, digits = 1),
              #"Greece"        = ~ n_perc0(Country == 6, digits = 1),
              "Germany"        = ~ n_perc0(Country == 7, digits = 1),
              "Denmark"        = ~ n_perc0(Country == 9, digits = 1)),
       
       "Tumor site" =
         list("Proximal colon"  = ~ n_perc0(location == 1, na_rm = T, digits = 1),
              "Distal colon"    = ~ n_perc0(location == 2, na_rm = T, digits = 1),
              "Rectum"          = ~ n_perc0(location == 3, na_rm = T, digits = 1),
              "Other"           = ~ n_perc0(location == 4, na_rm = T, digits = 1),
              "Unknown"         = ~ n_perc0(is.na(location), digits = 1)),
       "Diagnosis with histological verification" =
         list("Yes"    = ~ n_perc0(Bdg1clrt_hist == 1, na_rm = T, digits = 1),
              "No"     = ~ n_perc0(Bdg1clrt_hist == 0, na_rm = T, digits = 1)),
       
       "Smoking status" =
         list("Non smoker"     = ~ n_perc0(Smoke_Stat == 1, digits = 1),
              "Never smoker"   = ~ n_perc0(Smoke_Stat == 2, digits = 1),
              "Smoker"         = ~ n_perc0(Smoke_Stat == 3, digits = 1)),
       #"Unknown"      =  ~ n_perc0(Smoke_Stat == 4, digits = 1)),

       "Height (cm)" =
         list("Mean (SD)" = ~ mean_sd(Height_C, digits = 1)),
       "BMI (kg/m2)" =
         list("Mean (SD)" = ~ mean_sd(Bmi_C, digits = 1)),
       "Waist circumference (cm)" =
         list("Mean (SD)" = ~ mean_sd(Waist_C, digits = 1)),
       "Total energy intake (kCal)" =
         list("Mean (SD)" = ~ mean_sd(Qe_Energy, na_rm = T, show_n = "never", digits = 0)),

       "Physical activity" =
         list("Mean (SD)" = ~ mean_sd(Pa_Mets, na_rm = T, show_n = "never", digits = 1)),
       "Alcohol intake (g/day)" =
         list("Mean (SD)" = ~ mean_sd(Qe_Alc, na_rm = T, show_n = "never", digits = 1)),

       "WCRF score" = 
         list("Mean (SD)" = ~ mean_sd(Wcrf_C_Cal, na_rm = T, show_n = "never")),
       "Fatty acid metabolic signature" = 
         list("Mean (SD)" = ~ mean_sd(comp1, na_rm = T, show_n = "never")),
       "Endogenous metabolic signature" = 
         list("Mean (SD)" = ~ mean_sd(comp2, na_rm = T, show_n = "never"))
  )

st1 <- summary_table(group_by(crc.sig, Cncr_Caco_Clrt), crc.sum)
print(st1, cnames = c("Controls", "Cases"))
# Copy and paste output into an Rmarkdown file and render to word/pdf etc
# Note: for some reason Tfollowup didn't work. Calculated mean and SD manually

# Tests for p-values for table -----------------------------
# McNemar and Wilcoxon signed rank test?

# Arrange df to pair cases and controls
df1 <- crc.sig %>% arrange(Cncr_Caco_Clrt, Match_Caseset.x)

# Case/control models: CRC A
t.test(df1$Height_C ~ df1$Cncr_Caco_Clrt, paired = T)$p.value
t.test(df1$Age_Blood ~ df1$Cncr_Caco_Clrt, paired = T)$p.value
t.test(df1$Bmi_C ~ df1$Cncr_Caco_Clrt, paired = T)$p.value
t.test(df1$Waist_C ~ df1$Cncr_Caco_Clrt, paired = T)$p.value
chisq.test(df1$Cncr_Caco_Clrt, df1$Smoke_Stat)$p.value


# Contain NAs
df1a <- df1 %>% group_by(Match_Caseset) %>% filter(!anyNA(Qe_Energy))
t.test(df1a$Qe_Energy ~ df1a$Cncr_Caco_Clrt, paired = T)$p.value

df1b <- df1 %>% group_by(Match_Caseset) %>% filter(!anyNA(Qe_Alc))
wilcox.test(df1b$Qe_Alc ~ df1b$Cncr_Caco_Clrt, paired = T)$p.value

df1c <- df1 %>% group_by(Match_Caseset) %>% filter(!anyNA(Pa_Mets))
t.test(df1c$Pa_Mets ~ df1c$Cncr_Caco_Clrt, paired = T)$p.value

df1d <- df1 %>% group_by(Match_Caseset) %>% filter(!anyNA(Wcrf_C_Cal))
t.test(df1d$Wcrf_C_Cal ~ df1d$Cncr_Caco_Clrt, paired = T)$p.value


df2 <- crc2a %>% arrange(Cncr_Caco_Clrt, Match_Caseset)

# Case/control models: CRC B
t.test(df2$Height_C ~ df2$Cncr_Caco_Clrt, paired = T)$p.value
t.test(df2$Age_Blood ~ df2$Cncr_Caco_Clrt, paired = T)$p.value
t.test(df2$Bmi_C ~ df2$Cncr_Caco_Clrt, paired = T)$p.value
t.test(df2$Waist_C ~ df2$Cncr_Caco_Clrt, paired = T)$p.value
chisq.test(df2$Cncr_Caco_Clrt, df2$Smoke_Stat)$p.value

# Contain NAs
df2a <- df2 %>% group_by(Match_Caseset) %>% filter(!anyNA(Qe_Energy))
t.test(df2a$Qe_Energy ~ df2a$Cncr_Caco_Clrt, paired = T)$p.value

df2b <- df2 %>% group_by(Match_Caseset) %>% filter(!anyNA(Qe_Alc))
wilcox.test(df2b$Qe_Alc ~ df2b$Cncr_Caco_Clrt, paired = T)$p.value

df2c <- df2 %>% group_by(Match_Caseset) %>% filter(!anyNA(Pa_Mets))
t.test(df2c$Pa_Mets ~ df2c$Cncr_Caco_Clrt, paired = T)$p.value

df2d <- df2 %>% group_by(Match_Caseset) %>% filter(!anyNA(Wcrf_C_Cal))
t.test(df2d$Wcrf_C_Cal ~ df2d$Cncr_Caco_Clrt, paired = T)$p.value
