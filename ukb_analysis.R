# Read in ----
#UKB data: different ways
#only works for Stat 12 files
library(foreign)
read.dta("D://ukb_summanalytical2017.dta")

#works for Stata 13 files
library(readstata13)
dat <- read.dta13("D://ukb_summanalytical2017.dta")
write.csv(dat, "D://ukb_summanalytical2017.csv")

#also works for Stata 13 files
library(haven)
ukb <- read_dta("D://ukb_summanalytical2017.dta") #original data
ukb.stata <- read_dta("D://ukb_cox_analyses_2017.dta") #renamed and recoded for cox model

#or, exported to csv and read in with readr (could be faster)
library(readr)
dat1 <- read_csv("D://ukb_summanalytical2017.csv")

# Data preparation for Cox model ----

ukbprep <- function() {
  
  #saveRDS(ukb, file="D:/UKB Cox dataset.rds")
  #saveRDS(ukb, file="UKB Cox dataset.rds")
  ukb <- readRDS("UKB_Cox_dataset.rds") #472526 observations
  
  #get dropped IDs from file
  droppedIDs <- scan("w25897_20180503.csv")

  library(tidyverse)
  ukb <- ukb %>% filter(!(eID %in% droppedIDs))

  #drop diabetes unknown (9) cases (n = 1956)
  ukb <- ukb %>% filter(!(eID %in% droppedIDs) & diabet != 9)
  
  #assign 6 to PA category unknown
  ukb$pa_met_cats[is.na(ukb$pa_met_cats)] <- 6
  
  #ukb.stata <- ukb.stata %>% filter(!(eID %in% droppedIDs))
  
  #Lifestyle variables used:
  #pa_total_mets: total physical activity level (MET hr/week)
  #pa_met_cats: category of total physical activity level
  #colorectal_inc: incident cancers of the colon and rectum
  #pa5cat: median value of total physical activity level for each category
  
  #Time variables:
  #age_exit_first: age at exit
  #agebirth: 0 for all observations
  #age_recr: age at recruitment
  
  #list categorical variables to be converted to factors. colorectal_inc stays numeric or produces error
  var.list <- c("pa_met_cats", "diabet", "qualif", "alc_freq", "ever_horm", "smoke_intensity", "redprocmeat_cat", 
                "fh_crc", "aspibu_use", "sex", "agecat", "q5town", "region", "bmi")
  
  #remove observations missing derived pa_activity MET variable
  ukb <- ukb %>% 
    #filter(!is.na(pa_total_mets)) %>% 
    mutate_at(vars(var.list), as.factor)
        
  #Check number of observations and count female (0) and male (1)
  #nrow(ukb)
  #ukb$sex
  #table(ukb$sex) #254846 Women 217700 Men 
  #ukb %>% group_by(bmi, diabet, sex) %>% summarise(n = n())
  
  #BMI
  #table(ukb.stata$bmi)
  #table(ukb$bmi)
  #2443 underweight (1)
  #153803 normal (2)
  #200917 overweight (3)
  #115360 obese (4)
  #table(ukb$sex, ukb$bmi)
  #table(ukb$sex, ukb$diabet)

  #Diabetes
  #table(ukb.stata$diabet) #0 No 1 Yes 9 Missing. 21922 cases, 1956 missing
  
  #get median of total physical activity level for each category
  pameds <- tapply(ukb$pa_total_mets, ukb$pa_met_cats, median) %>% as.numeric
  
  #create variable pa5cat with medians from tapply above (also use modify() from purrr?)
  ukb$pa5cat <- pameds[as.factor(ukb$pa_met_cats)]
  return(ukb)
}
ukb <- ukbprep()

# Cox proportional hazards models ----

# Create the survival object. In stata, age at exit, age at birth, age at recruitment, and CRC yes/no are used
library(survival)
survobj <- Surv(time = ukb$age_recr, time2 = ukb$age_exit_frst, event = ukb$colorectal_inc)

base <- survobj ~ diabet + fh_crc + alc_freq + smoke_intensity + qualif + redprocmeat_cat + 
          aspibu_use + ever_horm + height_m + pa_met_cats + strata(sex, agecat, q5town, region)

# Run models. 3 observations have exit time < start time. Order roughly by significance of co-variates

# All subjects, female only, male only, normal, overweight, obese.

fit0 <- coxph(update(base, .~. + bmi), data = ukb)
fit1 <- coxph(update(base, .~. + bmi), data = ukb, subset = sex == 0)
fit2 <- coxph(update(base, .~. + bmi), data = ukb, subset = sex == 1)
fit3 <- coxph(base,                    data = ukb, subset = bmi == 2)
fit4 <- coxph(base,                    data = ukb, subset = bmi == 3)
fit5 <- coxph(base,                    data = ukb, subset = bmi == 4)

library(broom)
t1 <- map_df(list(fit0, fit1, fit2, fit3, fit4, fit5), tidy) %>% filter(term == "diabet1")

library(metafor)
par(mar=c(5,4,2,2))
forest(t1$estimate, ci.lb = t1$conf.low, ci.ub = t1$conf.high,
       refline = 1, xlab = "Hazard ratio (per unit increase in score)", pch = 18, 
       transf = exp, psize = 1.5, 
       slab = c("All", "Female", "Male", "Normal", "Overweight", "Obese")) 
#ilab = studies[, 2:3], 
#rows = c(1:3, 5:7, 9:11),
#ylim = c(0, 14),
#ilab.pos = 4, ilab.xpos = c(-0.8, -0.6), 
#xlim = c(-1.2, 2))

summary(fit0)

