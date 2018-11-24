# Preparation of CRC case-control datasets

# Read in small CRC case control dataset (1) and subset the 146 biocrates compounds
# stored at \\inti\NME\EPIC_Projects\Epic_Colonrectum\Nested_CaCo_Study\2016
# Missings are already imputed

prepcrc1 <- function(){

  library(tidyverse)
  library(haven)
  crc <- read_sas("clrt_caco_metabo.sas7bdat")
  
  # Metadata and WCRF scores (keep on local drive because big file)
  wcrf <- read_dta("Wcrf_Score.dta") %>% select(Idepic, Wcrf_C_Cal)
  
  # Remove duplicated Idepics (with dplyr or base)
  meta <- read_dta("D:/clrt_caco.dta") %>% select(-Match_Caseset, -Country, -Center, -Cncr_Caco_Clrt) %>%
    distinct(Idepic, .keep_all = T)
  
  # meta.dedup <- meta[!duplicated(meta["Idepic"]), ]
  
  # Converted to dta
  # write_dta(data, "clrt_caco_metabo.dta")
  
  # Subset biocrates data by taking !is.na > 0
  concs <- crc %>% select(Acylcarn_C0 : Sugars_H1, -Batch_MetBio)
  biocrates <- apply(concs, 1, function(x) sum(!is.na(x)) > 0)
  
  # Subset Biocrates subjects only
  crc1 <- crc[biocrates, ]
  crc1 <- crc1 %>% inner_join(meta, by = "Idepic")
  
  var.list <- c("Country", "Center", "Sex", "Match_Caseset", "Smoke_Stat", "L_School")
  crc1 <- crc1 %>% mutate_at(vars(var.list), as.factor)
  
  #class(crc1$Sex)

}
crc1 <- prepcrc1()

prepcrc2 <- function(){
  library(haven)
  wcrf <- read_dta("Wcrf_Score.dta") %>% select(ends_with("_Cal"), Idepic)
  crc2 <- read_csv("biocrates_p150.csv")
  crc2 <- crc2 %>% left_join(wcrf, by = "Idepic")
  
  var.list <- c("Country", "Center", "Sex", "Match_Caseset", "Smoke_Stat", "L_School")
  crc2 <- crc2 %>% mutate_at(vars(var.list), as.factor)
}
crc2 <- prepcrc2()

# glmnet penalised model -----------------------------------------------------------------------------------------------
# didn't really work but useful to know how to do

#data prep from plsdata. Convert high and low scores to 1 and 0
filtmat <- plsdata[ plsdata[,1] != 3, ]
filtmat[1] <- ifelse(filtmat[1] > 3, 1, 0)
scorehighlow <- filtmat[, 1]
metabo <- as.matrix(filtmat[, -1])

library(glmnet)
# response and predictor variables must be matrices
fit <- glmnet(metabo, scorehighlow, family = "binomial")
summary(fit)

cv.fit <- cv.glmnet(metabo, scorehighlow, family = "binomial")
plot(cv.fit)
lambda.min <- cv.fit$lambda.min
#First dotted vertical line indicates minimal mse; 2nd indicates one sd from mse

plot(fit)
coef(fit)
obs2predict <- as.matrix(obs2predict)

#Predict high or low score from model
tpred <- predict(fit, obs2predict, family = "binomial", s = lambda.min, type = "class") %>% as.numeric
names(tpred) <- "glmscore"

# Bind predicted scores from glmnet model to the CRC dataset
crc <- bind_cols(glmscore = tpred, df)

# Remove odd cases or controls
crc <- crc %>% group_by(Match_Caseset) %>% filter(n() == 2) %>% ungroup

# CLR Predicted scores
library(survival)
fit.clr <- clogit(Cncr_Caco_Clrt ~ glmscore + strata(Match_Caseset), data = crc)
summary(fit.clr)
table(status = crc$Cncr_Caco_Clrt, scorehighlow = crc$glmscore)


