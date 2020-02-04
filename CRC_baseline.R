# Baseline characteristics for the 2 CRC metabolomics studies
# 8 unpaired samples have been removed for CRC1

source("CRC_prep_data.R")

crc1a <- crc1 %>% group_by(Match_Caseset) %>% filter(n() == 2) %>%
  select(Sex, Age_Blood, Height_C, Weight_C, Bmi_C, Qe_Energy, Country, Pa_Mets, Smoke_Stat, Qe_Alc,
         Wcrf_C_Cal, Cncr_Caco_Clrt) %>% mutate(Study = "CRC1")

crc2a <- crc2 %>% #left_join(wcrf, by = "Idepic") %>%
  select(Match_Caseset, Sex, Age_Blood, Height_C, Weight_C, Bmi_C, Qe_Energy, Country, Pa_Mets, Smoke_Stat, Qe_Alc,
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

# Arrange df to pair cases and controls
df1 <- crc1a %>% arrange(Cncr_Caco_Clrt, Match_Caseset)

# Case/control models: CRC A
t.test(df1$Height_C ~ df1$Cncr_Caco_Clrt, paired = T)$p.value
t.test(df1$Weight_C ~ df1$Cncr_Caco_Clrt, paired = T)$p.value
t.test(df1$Age_Blood ~ df1$Cncr_Caco_Clrt, paired = T)$p.value
t.test(df1$Bmi_C ~ df1$Cncr_Caco_Clrt, paired = T)$p.value
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
t.test(df2$Weight_C ~ df2$Cncr_Caco_Clrt, paired = T)$p.value
t.test(df2$Age_Blood ~ df2$Cncr_Caco_Clrt, paired = T)$p.value
t.test(df2$Bmi_C ~ df2$Cncr_Caco_Clrt, paired = T)$p.value
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
