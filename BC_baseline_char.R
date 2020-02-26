# Generate R markdown for baseline characteristics table
# To add something about triple negative?
# Source prep data to remove the 10 subjects not matched on menopausal status
source("BC_prep_data.R")

library(tidyverse)
meta <- read_csv("Lifepath_meta.csv", na = "9999") %>% filter(!(MATCH %in% unmatch_pairs)) %>%
  select(CT,                    # case-control status
         AGE,                   
         BMICat1, 
         RTHCat1,               # Waist-hip ratio categorical
         MENOPAUSE,             # menopausal status at blood collection
         SMK, 
         DIABETE, 
         Life_Alcohol_Pattern_1,
         ALCOHOL,
         BP, 
         CO,                    # previous oral contraceptive use
         Trait_Horm,            # menopausal treatment therapy taken 24h before blood collection
         #DURTHSDIAG,            # Duration of use of therapy at date of diagnosis
         DURTHSBMB,             # Duration of use of therapy at blood collection
         #CENTTIMECat1,          # time before centrifugation (?) OMIT not of interest
         FASTING, 
         STOCKTIME,             # Storage time (years)
         BEHAVIOUR,             # Tumour behaviour
         SUBTYPE, 
         #CERB2,                # HER2 receptor OMIT mostly unknown
         ER,                    # Estrogen receptor
         PR,                    # Progesterone receptor
         SBR, 
         GRADE, 
         STADE, 
         DIAGSAMPLINGCat1) %>% 
  mutate_at(vars(-AGE, -STOCKTIME, -ALCOHOL, -DURTHSBMB), as.factor)


library(kableExtra)
library(qwraps2)
options(qwraps2_markup = "markdown")

# Manually generated ----
summ <-
  list(#"Total subjects"                  = list("N"   =   ~ n()),
    
       "Age at blood collection (years)" = list("Mean"     =  ~ mean_sd(AGE)),
       
     "Menopausal status at blood collection" =
         
         list("Pre-menopausal"   =  ~ n_perc0(MENOPAUSE == 0, digits = 1),
              "Post-menopausal"  =  ~ n_perc0(MENOPAUSE == 1, digits = 1)),
       
     "BMI" =    
       
        list("Underweight or normal"   =  ~ n_perc0(BMICat1 == 1, na_rm = T, digits = 1),
             "Overweight"              =  ~ n_perc0(BMICat1 == 2, na_rm = T, digits = 1),
             "Obese"                   =  ~ n_perc0(BMICat1 == 3, na_rm = T, digits = 1),
             "Unknown"                 =  ~ n_perc0(is.na(BMICat1), digits = 1)),
       
     "Waist to hip ratio" =           
       
        list("< 0.8" = ~ n_perc0(RTHCat1 == 0, na_rm = T, digits = 1),
             "> 0.8" = ~ n_perc0(RTHCat1 == 1, na_rm = T, digits = 1),
             "Unknown" = ~ n_perc0(is.na(RTHCat1), digits = 1)),
       
     "Smoking status" = 
       
        list("Yes" = ~ n_perc0(SMK == 1, digits = 1),
              "No"  = ~ n_perc0(SMK == 0, digits = 1)),
       
     "Diabetic status"  = 
       
        list("Yes" = ~ n_perc0(DIABETE == 1, digits = 1),
              "No"  = ~ n_perc0(DIABETE == 0, digits = 1)),
       
     "Lifetime alcohol drinking pattern" =
       
        list("Non-consumers (0 g/day)"      =  ~ n_perc0(Life_Alcohol_Pattern_1 == 0, na_rm = T, digits = 1),
              "Light consumers (1-10 g/day)" =  ~ n_perc0(Life_Alcohol_Pattern_1 == 1, na_rm = T, digits = 1),
              "Drinkers (>10 g/day)"         =  ~ n_perc0(Life_Alcohol_Pattern_1 == 2, na_rm = T, digits = 1),
              "Unknown"                      =  ~ n_perc0(is.na(Life_Alcohol_Pattern_1), digits = 1)),
     
     "Alcohol intake (g/day)" =
       
       list("Mean (SD)" = ~ mean_sd(ALCOHOL)),
       
     "Blood pressure" =
       
        list("Normal tension"  =  ~ n_perc0(BP == 0, na_rm = T, digits = 1),
              "Hypertension"    =  ~ n_perc0(BP == 1, na_rm = T, digits = 1),
              "No information"  =  ~ n_perc0(is.na(BP), digits = 1)),
       
     "Previous oral contraceptive use" =
       
        list("Yes" = ~ n_perc0(CO == 1, digits = 1),
              "No" = ~ n_perc0(CO == 0, digits = 1)),
       
     "Menopause hormone therapy at blood collection" =
       
        list("Yes" =  ~ n_perc0(Trait_Horm == 1, digits = 1),
             "No"  =  ~ n_perc0(Trait_Horm == 0, digits = 1)),
       
     "Duration of use of menopause hormonal treatment at baseline" =
       
        list("Mean (SD)" =  ~ mean_sd(DURTHSBMB)),
       
     "Fasting status" =
       
        list("Fasting"     = ~ n_perc0(FASTING == 1, digits = 1),
              "Non-fasting" = ~ n_perc0(FASTING == 0, digits = 1)),
       
     "Biobank storage time (years)" =
       
        list("Mean (SD)" = ~ mean_sd(STOCKTIME)),
       
     "Time between sampling and diagnosis" =
         
        list("5 years or less"   =  ~ n_perc0(DIAGSAMPLINGCat1 == 1, na_rm = T, digits = 1),
              "More than 5 years" =  ~ n_perc0(DIAGSAMPLINGCat1 == 2, na_rm = T, digits = 1)),
              #"No information"  =  ~ n_perc0(is.na(DIAGSAMPLINGCat1))),
       
     "Tumor behavior" =
       
       list("In situ"  = ~ n_perc0(BEHAVIOUR == 2, na_rm = T, digits = 1),
            "Invasive" = ~ n_perc0(BEHAVIOUR == 3, na_rm = T, digits = 1),
            "Unknown"  = ~ n_perc0(is.na(BEHAVIOUR), digits = 1)),
       
     "Subtype" =
       
       list("Lobular"  =  ~ n_perc0(SUBTYPE == 1, na_rm = T, digits = 1),
            "Ductal"   =  ~ n_perc0(SUBTYPE == 2, na_rm = T, digits = 1),
            "Tubular"  =  ~ n_perc0(SUBTYPE == 3, na_rm = T, digits = 1),
            "Mixed"    =  ~ n_perc0(SUBTYPE == 4, na_rm = T, digits = 1),
            "Others"   =  ~ n_perc0(SUBTYPE == 5, na_rm = T, digits = 1),
            "Unknown"  =  ~ n_perc0(is.na(SUBTYPE), digits = 1)),
       
     #"HER2" =
     
     # list("Negative"  =  ~ n_perc0(CERB2 == 1, na_rm = T),
     #      "Positive"  =  ~ n_perc0(CERB2 == 2, na_rm = T),
     #      "Unknown"   =  ~ n_perc0(is.na(CERB2)),
       
     "Estrogen receptor" =
       
      list("Negative"   =  ~ n_perc0(ER == 0, na_rm = T, digits = 1),
            "Positive"  =  ~ n_perc0(ER == 1, na_rm = T, digits = 1),
            "Unknown"   =  ~ n_perc0(is.na(ER), digits = 1)),
       
     "Progesterone receptor" =
       
      list("Negative"  =  ~ n_perc0(PR == 0, na_rm = T, digits = 1),
            "Positive" =  ~ n_perc0(PR == 1, na_rm = T, digits = 1),
            "Unknown"  =  ~ n_perc0(is.na(PR), digits = 1)),
       
     "SBR Grade" =
       
      list("Favorable prognosis"      =  ~ n_perc0(SBR == 1, na_rm = T, digits = 1),
            "Intermediate prognosis"  =  ~ n_perc0(SBR == 2, na_rm = T, digits = 1),
            "Unfavorable prognosis"   =  ~ n_perc0(SBR == 3, na_rm = T, digits = 1),
            "No information"          =  ~ n_perc0(is.na(SBR), digits = 1)),
       
     "Grade" =
       
       list("1"       =  ~ n_perc0(GRADE == 1, na_rm = T, digits = 1),
            "2"       =  ~ n_perc0(GRADE == 2, na_rm = T, digits = 1),
            "3"       =  ~ n_perc0(GRADE == 3, na_rm = T, digits = 1),
            "Unknown" =  ~ n_perc0(is.na(GRADE), digits = 1)),
       
     "Stade" =
       
       list("1"    =  ~ n_perc0(STADE == 1, na_rm = T, digits = 1),
            "2"    =  ~ n_perc0(STADE == 2, na_rm = T, digits = 1),
            "3"    =  ~ n_perc0(STADE == 3, na_rm = T, digits = 1),
            "4"    =  ~ n_perc0(STADE == 4, na_rm = T, digits = 1),
            "No information"  =  ~ n_perc0(is.na(STADE), digits = 1))
       
  )

# Automatic
#summ <- meta0 %>% qsummary(., numeric_summaries = list("Mean (SD)" = "~ mean_sd(%s)"),
#                                   n_perc_args = list(digits = 1, show_symbol = T))

# ----
st <- summary_table(group_by(meta, CT), summ)
print(st, cnames = c("Controls (N=786)", "Cases (N=786)"))
# Copy and paste output into an Rmarkdown file and render to word/pdf etc

# ---------------------------------------------------------------

# Models for p-values for table
# Replace t-test and chi-sq with McNemar and Wilcoxon signed rank?

# Case/control models
ll <- list(
  chisq.test(meta$CT, meta$RTHCat1),
  chisq.test(meta$CT, meta$BMICat1),
  chisq.test(meta$CT, meta$SMK),
  chisq.test(meta$CT, meta$DIABETE),
  chisq.test(meta$CT, meta$BP),
  chisq.test(meta$CT, meta$Life_Alcohol_Pattern_1),
  chisq.test(meta$CT, meta$Trait_Horm)
)

ll <- list(
  mcnemar.test(meta$CT, meta$RTHCat1),
  mcnemar.test(meta$CT, meta$BMICat1),
  mcnemar.test(meta$CT, meta$SMK),
  mcnemar.test(meta$CT, meta$DIABETE),
  mcnemar.test(meta$CT, meta$BP),
  mcnemar.test(meta$CT, meta$Life_Alcohol_Pattern_1),
  mcnemar.test(meta$CT, meta$Trait_Horm)
)




#t.test(meta$DURTHSDIAG ~ meta$CT)$p.value
test1 <- wilcox.test(meta$DURTHSDIAG ~ meta$CT)


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



