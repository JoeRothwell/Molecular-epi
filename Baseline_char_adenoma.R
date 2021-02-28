# Baseline characteristics for adenoma project. Split by country
source("adenoma_crc.R")

library(qwraps2)
options(qwraps2_markup = "markdown")

# Manually specified summary
summ <-
  list("Sex" = 
         list("Male"       =  ~ n_perc0(sex == "M", digits = 1),
              "Female"     =  ~ n_perc0(sex == "F", digits = 1)),
       "Age at blood collection (years)" = 
         list("Mean"     =  ~ mean_sd(age, digits = 1)),
       
       "Diagnosed pathology" =
         list("Colorectal cancer"  = ~ n_perc0(pathsum == 1, digits = 1),
              "High grade dysplasia" = ~ n_perc0(pathsum == 2, digits = 1),
              "Adenoma"          = ~ n_perc0(pathsum == 3, digits = 1),
              "Polyp"           = ~ n_perc0(pathsum == 4, digits = 1),
              "Normal (after colonoscopy)"  = ~ n_perc0(pathsum == 5, digits = 1),
              "Healthy blood donor"   = ~ n_perc0(pathsum == 6, digits = 1)),
       
       "Prevalent diabetes" =
         list("Yes"    = ~ n_perc0(diabetes == 1, na_rm = T, digits = 1),
              "No"     = ~ n_perc0(diabetes == 0, na_rm = T, digits = 1)),
       
       "Smoking status" =
         list("Non smoker"     = ~ n_perc0(smoke == 1, na_rm = T, digits = 1),
              "Never smoker"   = ~ n_perc0(smoke == 0, na_rm = T, digits = 1),
              "Smoker"         = ~ n_perc0(smoke == 2, na_rm = T, digits = 1)),
       
       "BMI (kg/m2)" =
         list("Mean (SD)" = ~ mean_sd(bmi, na_rm = T, digits = 1)),
       
       #"Physical activity" =
        # list("Mean (SD)" = ~ mean_sd(Pa_Mets, na_rm = T, show_n = "never", digits = 1)),
       
       "Alcohol drinker" =
         list("Yes"    = ~ n_perc0(alcohol == 1, na_rm = T, digits = 1),
              "No"     = ~ n_perc0(alcohol == 0, na_rm = T, digits = 1)),
       
       "Alcoholic drinks per week" =
         list("Mean (SD)" = ~ mean_sd(alcohol_drinks_week, na_rm = T, show_n = "never", digits = 1))
  )

sumtab <- summary_table(group_by(mat, country), summ)
print(sumtab, cnames = c("Ireland", "Czech Republic"))
