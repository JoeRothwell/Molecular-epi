# Generate R markdown for baseline characteristics table

library(tidyverse)
meta <- read_csv("Lifepath_meta.csv", na = "9999") %>%
  select(CT,                    # case-control status
         AGE,                   
         BMICat1, 
         RTHCat1,               # Waist-hip ratio categorical
         MENOPAUSE,             # menopausal status at blood collection
         SMK, 
         DIABETE, 
         Life_Alcohol_Pattern_1, 
         BP, 
         CO,                    # previous oral contraceptive use
         Trait_Horm,            # menopausal treatment therapy taken 24h before blood collection
         DURTHSDIAG,            # Duration of use of therapy at date of diagnosis
         #CENTTIMECat1,          # time before centrifugation (?) OMIT not of interest
         FASTING, 
         STOCKTIME,             # Storage time (years)
         BEHAVIOUR,             # Tumour behaviour
         SUBTYPE, 
         #CERB2,                 # HER2 receptor OMIT mostly unknown
         ER,                    # Estrogen receptor
         PR,                    # Progesterone receptor
         SBR, 
         GRADE, 
         STADE, 
         DIAGSAMPLINGCat1) %>% 
  mutate_at(vars(-AGE, -STOCKTIME, -DURTHSDIAG), as.factor)

library(kableExtra)
library(qwraps2)
options(qwraps2_markup = "markdown")

# Manually generated ----
summ <-
  list(#"Total subjects"                  = list("N"   =   ~ n()),
       "Age at blood collection (years)" = list("Mean"     =  ~ mean_sd(AGE)),
       "BMI" =                          
          list("Underweight or normal"   =  ~ n_perc0(BMICat1 == 1, na_rm = T),
               "Overweight"              =  ~ n_perc0(BMICat1 == 2, na_rm = T),
               "Obese"                   =  ~ n_perc0(BMICat1 == 3, na_rm = T),
               "Unknown"                 =  ~ n_perc0(is.na(BMICat1))),
       "Waist to hip ratio" =             
                                           list("< 0.8"  =  ~ n_perc0(RTHCat1 == 0, na_rm = T),
                                                "> 0.8"  =  ~ n_perc0(RTHCat1 == 1, na_rm = T),
                                                "Unknown"   =  ~ n_perc0(is.na(RTHCat1))),
       "Menopausal status at blood collection" =
                                           list("Pre-menopausal"   =  ~ n_perc0(MENOPAUSE == 0),
                                                "Post-menopausal"  =  ~ n_perc0(MENOPAUSE == 1)),
       "Smoking status"                  = list("Yes"  =  ~ n_perc0(SMK == 1),
                                                "No"   =  ~ n_perc0(SMK == 0)),
       "Diabetic status"                 = list("Yes"   =  ~ n_perc0(DIABETE == 1),
                                                "No"  =  ~ n_perc0(DIABETE == 0)),
       "Lifetime alcohol drinking pattern" =
         list("Non-consumers (0 g/day)"        =  ~ n_perc0(Life_Alcohol_Pattern_1 == 0, na_rm = T),
              "Light consumers (1-10 g/day)"   =  ~ n_perc0(Life_Alcohol_Pattern_1 == 1, na_rm = T),
              "Drinkers (>10 g/day)"           =  ~ n_perc0(Life_Alcohol_Pattern_1 == 2, na_rm = T),
              "Unknown"                        =  ~ n_perc0(is.na(Life_Alcohol_Pattern_1))),
        "Blood pressure" =
         list("Normal tension"       =  ~ n_perc0(BP == 0, na_rm = T),
              "Hypertension"         =  ~ n_perc0(BP == 1, na_rm = T),
              "No information"       =  ~ n_perc0(is.na(BP))),
       "Previous oral contraceptive use" =
         list("Yes"  =  ~ n_perc0(CO == 1),
              "No"   =  ~ n_perc0(CO == 0)),
       "Menopause hormone therapy at blood collection" =
         list("Yes"  =  ~ n_perc0(Trait_Horm == 1),
              "No"   =  ~ n_perc0(Trait_Horm == 0)),
       "Duration of use of menopause hormonal treatment" =
         list("Mean (SD)"     =  ~ mean_sd(DURTHSDIAG)),
       #"Time before centrifugation" =
       #  list("< 12h"   =  ~ n_perc0(CENTTIMECat1 == 1, na_rm = T),
       #      "> 12-24h" =  ~ n_perc0(CENTTIMECat1 == 2, na_rm = T),
       #      "> 24h"    =  ~ n_perc0(CENTTIMECat1 == 3, na_rm = T),
       #      "Unknown"  =  ~ n_perc0(is.na(CENTTIMECat1)),
       "Fasting status" =
         list("Fasting"       =  ~ n_perc0(FASTING == 1),
              "Non-fasting"   =  ~ n_perc0(FASTING == 0)),
       "Biobank storage time (years)" =
         list("Mean (SD)" = ~ mean_sd(STOCKTIME)),
       "Time between sampling and diagnosis" =
         list("5 years or less"     =  ~ n_perc0(DIAGSAMPLINGCat1 == 1, na_rm = T),
              "More than 5 years"   =  ~ n_perc0(DIAGSAMPLINGCat1 == 2, na_rm = T)),
              #"No information"      =  ~ n_perc0(is.na(DIAGSAMPLINGCat1))),
       "Tumor behavior" =
         list("In situ"   =  ~ n_perc0(BEHAVIOUR == 2, na_rm = T),
              "Invasive"  =  ~ n_perc0(BEHAVIOUR == 3, na_rm = T),
              "Unknown"   =  ~ n_perc0(is.na(BEHAVIOUR))),
       "Subtype" =
         list("Lobular"  =  ~ n_perc0(SUBTYPE == 1, na_rm = T),
              "Ductal"   =  ~ n_perc0(SUBTYPE == 2, na_rm = T),
              "Tubular"  =  ~ n_perc0(SUBTYPE == 3, na_rm = T),
              "Mixed"    =  ~ n_perc0(SUBTYPE == 4, na_rm = T),
              "Others"   =  ~ n_perc0(SUBTYPE == 5, na_rm = T),
              "Unknown"  =  ~ n_perc0(is.na(SUBTYPE))),
       #"HER2" =
       # list("Negative"  =  ~ n_perc0(CERB2 == 1, na_rm = T),
        #     "Positive"  =  ~ n_perc0(CERB2 == 2, na_rm = T),
        #     "Unknown"   =  ~ n_perc0(is.na(CERB2)),
       "Estrogen receptor" =
         list("Negative"  =  ~ n_perc0(ER == 0, na_rm = T),
              "Positive"  =  ~ n_perc0(ER == 1, na_rm = T),
              "Unknown"   =  ~ n_perc0(is.na(ER))),
       "Progesterone receptor" =
         list("Negative" =  ~ n_perc0(PR == 0, na_rm = T),
              "Positive" =  ~ n_perc0(PR == 1, na_rm = T),
              "Unknown"  =  ~ n_perc0(is.na(PR))),
       "SBR Grade" =
         list("Favorable prognosis"     =  ~ n_perc0(SBR == 1, na_rm = T),
              "Intermediate prognosis"  =  ~ n_perc0(SBR == 2, na_rm = T),
              "Unfavorable prognosis"   =  ~ n_perc0(SBR == 3, na_rm = T),
              "No information"          =  ~ n_perc0(is.na(SBR))),
       "Grade" =
         list("1"       =  ~ n_perc0(GRADE == 1, na_rm = T),
              "2"       =  ~ n_perc0(GRADE == 2, na_rm = T),
              "3"       =  ~ n_perc0(GRADE == 3, na_rm = T),
              "Unknown" =  ~ n_perc0(is.na(GRADE))),
       "Stade" =
         list("1"    =  ~ n_perc0(STADE == 1, na_rm = T),
              "2"    =  ~ n_perc0(STADE == 2, na_rm = T),
              "3"    =  ~ n_perc0(STADE == 3, na_rm = T),
              "4"    =  ~ n_perc0(STADE == 4, na_rm = T),
              "No information"  =  ~ n_perc0(is.na(STADE)))
  )

# Automatic
#summ <- meta0 %>% qsummary(., numeric_summaries = list("Mean (SD)" = "~ mean_sd(%s)"),
#                                   n_perc_args = list(digits = 1, show_symbol = T))

# ----
st <- summary_table(group_by(meta, CT), summ)
# Copy and paste output into an Rmarkdown file and render to word/pdf etc
