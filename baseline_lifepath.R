# Generate R markdown for baseline characteristics table

library(tidyverse)
meta <- read.csv("Lifepath_meta.csv")
meta0 <- meta %>% 
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
         CENTTIMECat1,          # time before centrifugation (?)
         FASTING, 
         STOCKTIME,             # Storage time (years)
         BEHAVIOUR,             # Tumour behaviour
         SUBTYPE, 
         CERB2,                 # HER2 receptor
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

# Automatic
#summ <- meta0 %>% qsummary(., numeric_summaries = list("Mean (SD)" = "~ mean_sd(%s)"),
#                                   n_perc_args = list(digits = 1, show_symbol = TRUE))

# Manually generated ----
summ <-
  list("Total subjects"                  = list("N"   =   ~ n()),
       "Age at blood collection (years)" = list("Mean"     =  ~ mean_sd(AGE)),
       "BMI" =                          
          list("Underweight or normal"   =  ~ n_perc0(BMICat1 == 1),
               "Overweight"              =  ~ n_perc0(BMICat1 == 2),
               "Obese"                   =  ~ n_perc0(BMICat1 == 3),
               "Unknown"                 =  ~ n_perc0(BMICat1 == 9999)),
       "Waist to hip ratio" =             
                                           list("< 0.8"  =  ~ n_perc0(RTHCat1 == 0),
                                                "> 0.8"  =  ~ n_perc0(RTHCat1 == 1),
                                                "Unknown"   =  ~ n_perc0(RTHCat1 == 9999)),
       "Menopausal status at blood collection" =
                                           list("Pre-menopausal"   =  ~ n_perc0(MENOPAUSE == 0),
                                                "Post-menopausal"  =  ~ n_perc0(MENOPAUSE == 1)),
       "Smoking status"                  = list("Yes"  =  ~ n_perc0(SMK == 1),
                                                "No"   =  ~ n_perc0(SMK == 0)),
       "Diabetic status"                 = list("Yes"   =  ~ n_perc0(DIABETE == 1),
                                                "No"  =  ~ n_perc0(DIABETE == 0)),
       "Lifetime alcohol drinking pattern" =
         list("Non-consumers (0 g/day)"        =  ~ n_perc0(Life_Alcohol_Pattern_1 == 0),
              "Light consumers (1-10 g/day)"   =  ~ n_perc0(Life_Alcohol_Pattern_1 == 1),
              "Drinkers (>10 g/day)"           =  ~ n_perc0(Life_Alcohol_Pattern_1 == 2),
              "Unknown"                        =  ~ n_perc0(Life_Alcohol_Pattern_1 == 9999)),
        "Blood pressure" =
         list("Normal tension"       =  ~ n_perc0(BP == 0),
              "Hypertension"         =  ~ n_perc0(BP == 1),
              "No information"       =  ~ n_perc0(BP == 9999)),
       "Previous oral contraceptive use" =
         list("Yes"  =  ~ n_perc0(CO == 1),
              "No"   =  ~ n_perc0(CO == 0)),
       "Menopause hormone therapy at blood collection" =
         list("Yes"  =  ~ n_perc0(Trait_Horm == 1),
              "No"   =  ~ n_perc0(Trait_Horm == 0)),
       "Duration of use of menopause hormonal treatments" =
         list("Mean"     =  ~ mean_sd(DURTHSDIAG)),
       "Time before centrifugation" =
          list("< 12h"   =  ~ n_perc0(CENTTIMECat1 == 1),
              "> 12-24h" =  ~ n_perc0(CENTTIMECat1 == 2),
              "> 24h"    =  ~ n_perc0(CENTTIMECat1 == 3),
              "Unknown"  =  ~ n_perc0(CENTTIMECat1 == 9999)),
       "Fasting status" =
         list("Fasting" =  ~ n_perc0(FASTING == 1),
              "Non-fasting"     =  ~ n_perc0(FASTING == 0)),
       "Storage time (years)" =
         list("Mean (SD)" = ~ mean_sd(STOCKTIME)),
       "Behavior" =
         list("In situ"   =  ~ n_perc0(BEHAVIOUR == 2),
              "Invasive"  =  ~ n_perc0(BEHAVIOUR == 3),
              "Unknown"   =  ~ n_perc0(BEHAVIOUR == 9999)),
       "Subtype" =
         list("Lobular"  =  ~ n_perc0(SUBTYPE == 1),
              "Ductal"  =  ~ n_perc0(SUBTYPE == 2),
              "Tubular"  =  ~ n_perc0(SUBTYPE == 3),
              "Mixed"  =  ~ n_perc0(SUBTYPE == 4),
              "Others"  =  ~ n_perc0(SUBTYPE == 5),
              "Unknown"   =  ~ n_perc0(SUBTYPE == 9999)),
       "HER2" =
         list("Negative"  =  ~ n_perc0(CERB2 == 1),
              "Positive"  =  ~ n_perc0(CERB2 == 2),
              "Unknown"   =  ~ n_perc0(CERB2 == 9999)),
       "Estrogen receptor" =
         list("Negative"  =  ~ n_perc0(ER == 0),
              "Positive"  =  ~ n_perc0(ER == 1),
              "Unknown"   =  ~ n_perc0(ER == 9999)),
       "Progesterone receptor" =
         list("Negative" =  ~ n_perc0(PR == 0),
              "Positive" =  ~ n_perc0(PR == 1),
              "Unknown"  =  ~ n_perc0(PR == 9999)),
       "SBR Grade" =
         list("Favorable prognosis"     =  ~ n_perc0(SBR == 1),
              "Intermediate prognosis"  =  ~ n_perc0(SBR == 2),
              "Unfavorable prognosis"   =  ~ n_perc0(SBR == 3),
              "No information"          =  ~ n_perc0(SBR == 9999)),
       "Grade" =
         list("1"       =  ~ n_perc0(GRADE == 1),
              "2"       =  ~ n_perc0(GRADE == 2),
              "3"       =  ~ n_perc0(GRADE == 3),
              "Unknown" =  ~ n_perc0(GRADE == 9999)),
       "Stade" =
         list("1"    =  ~ n_perc0(STADE == 1),
              "2"    =  ~ n_perc0(STADE == 2),
              "3"    =  ~ n_perc0(STADE == 3),
              "4"    =  ~ n_perc0(STADE == 4),
              "No information"  =  ~ n_perc0(STADE == 9999)),
       "Time between sampling and diagnosis" =
         list("5 years or less"     =  ~ n_perc0(DIAGSAMPLINGCat1 == 1),
              "More than 5 years"   =  ~ n_perc0(DIAGSAMPLINGCat1 == 2),
              "No information"      =  ~ n_perc0(DIAGSAMPLINGCat1 == 9999))
  )

# ----
st <- summary_table(group_by(meta0, CT), summ)
# Copy and paste output into an Rmarkdown file and render to word/pdf etc
