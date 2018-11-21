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

prepcrc2 <- function() {
  library(haven)
  wcrf <- read_dta("Wcrf_Score.dta") %>% select(ends_with("_Cal"), Idepic)
  crc2 <- read_csv("biocrates_p150.csv")
  crc2 <- crc2 %>% left_join(wcrf, by = "Idepic")
  
  var.list <- c("Country", "Center", "Sex", "Match_Caseset", "Smoke_Stat", "L_School")
  crc2 <- crc2 %>% mutate_at(vars(var.list), as.factor)
}
crc2 <- prepcrc2()


# Cox models for score ----


# Small case-control, test score and other co-variates
library(survival)

# 1: total; body fatness; physical activity; etc (see below)
# Adjusted for ENDB energy intake, highest level of schooling and smoking status
fit1 <- clogit(Cncr_Caco_Clrt ~ Wcrf_Bmi + Qe_Energy + L_School + Smoke_Stat +  strata(Match_Caseset), data = crc1)
fit2 <- clogit(Cncr_Caco_Clrt ~ Wcrf_Pa + Qe_Energy + L_School + Smoke_Stat + strata(Match_Caseset), data = crc1)
fit3 <- clogit(Cncr_Caco_Clrt ~ Wcrf_Fwg_Cal + Qe_Energy + L_School + Smoke_Stat + strata(Match_Caseset), data = crc1)
fit4 <- clogit(Cncr_Caco_Clrt ~ Wcrf_Pf_Cal + Qe_Energy + L_School + Smoke_Stat + strata(Match_Caseset), data = crc1)
fit5 <- clogit(Cncr_Caco_Clrt ~ Wcrf_Fv_Cal + Qe_Energy + L_School + Smoke_Stat + strata(Match_Caseset), data = crc1)
fit6 <- clogit(Cncr_Caco_Clrt ~ Wcrf_Fibt_Cal + Qe_Energy + L_School + Smoke_Stat + strata(Match_Caseset), data = crc1)
fit7 <- clogit(Cncr_Caco_Clrt ~ Wcrf_Meat_Cal + Qe_Energy + L_School + Smoke_Stat + strata(Match_Caseset), data = crc1)
fit0 <- clogit(Cncr_Caco_Clrt ~ Wcrf_C_Cal + Qe_Energy + L_School + Smoke_Stat + strata(Match_Caseset), data = crc1)


library(broom)
t1 <- map_df(list(fit1, fit2, fit3, fit4, fit5, fit6, fit7, fit0), tidy) %>% filter(grepl("Wcrf_", term))
scorecomp <- c("1. Body fatness", "2. Physical activity", "3. Energy density/sugary drinks",
          "4. FV intake", "5. Foods of plant origin", "6. Fibre intake", "7. Meat intake", 
          "Overall WCRF score (cal.)")


library(metafor)
par(mar=c(5,4,1,2))
forest(t1$estimate, ci.lb = t1$conf.low, ci.ub = t1$conf.high,
       refline = 1, xlab = "Odds ratio (per unit increase in score)", pch = 18, 
       transf = exp, psize = 1.5, xlim = c(-1.8, 4),
       slab = scorecomp) 
       #ilab = studies[, 2:3], rows = c(1:3, 5:7, 9:11),
       #ylim = c(0, 14),ilab.pos = 4, ilab.xpos = c(-0.8, -0.6))
hh <- par("usr")
text(hh[1], nrow(t1) + 2, "Component of score", pos = 4)
text(hh[2], nrow(t1) + 2, "OR [95% CI]", pos = 2)


#write_rds(crcB, "prepdata/CRC_smallerCC.rds")


# 1: total; body fatness; physical activity; etc (see below)
fit1 <- clogit(Cncr_Caco_Clrt ~ Wcrf_Bmi + Qe_Energy + L_School + Smoke_Stat + strata(Match_Caseset), data = crc2)
fit2 <- clogit(Cncr_Caco_Clrt ~ Wcrf_Pa + Qe_Energy + L_School + Smoke_Stat + strata(Match_Caseset), data = crc2)
fit3 <- clogit(Cncr_Caco_Clrt ~ Wcrf_Fwg_Cal + Qe_Energy + L_School + Smoke_Stat + strata(Match_Caseset), data = crc2)
fit4 <- clogit(Cncr_Caco_Clrt ~ Wcrf_Pf_Cal + Qe_Energy + L_School + Smoke_Stat + strata(Match_Caseset), data = crc2)
fit5 <- clogit(Cncr_Caco_Clrt ~ Wcrf_Fv_Cal + Qe_Energy + L_School + Smoke_Stat + strata(Match_Caseset), data = crc2)
fit6 <- clogit(Cncr_Caco_Clrt ~ Wcrf_Fibt_Cal + Qe_Energy + L_School + Smoke_Stat + strata(Match_Caseset), data = crc2)
fit7 <- clogit(Cncr_Caco_Clrt ~ Wcrf_Meat_Cal + Qe_Energy + L_School + Smoke_Stat + strata(Match_Caseset), data = crc2)
fit0 <- clogit(Cncr_Caco_Clrt ~ Wcrf_C_Cal + Qe_Energy + L_School + Smoke_Stat + strata(Match_Caseset), data = crc2)


library(broom)
t2 <- map_df(list(fit1, fit2, fit3, fit4, fit5, fit6, fit7, fit0), tidy) %>% filter(grepl("Wcrf_", term))


par(mar=c(5,4,1,2))
forest(t2$estimate, ci.lb = t2$conf.low, ci.ub = t2$conf.high,
       refline = 1, xlab = "Odds ratio (per unit increase in score)", pch = 18, 
       transf = exp, psize = 1.5, alim = c(0,2.5), steps = 6, xlim = c(-1.8, 4),
       slab = scorecomp) 
#ilab = studies[, 2:3], rows = c(1:3, 5:7, 9:11), ylim = c(0, 14),
#ilab.pos = 4, ilab.xpos = c(-0.8, -0.6), xlim = c(-1.2, 2))

hh <- par("usr")
text(hh[1], nrow(t2) + 2, "Component of score", pos = 4)
text(hh[2], nrow(t2) + 2, "OR [95% CI]", pos = 2)

# Both plots together
par(mfrow = c(2,1))
par(mar=c(5,4,2,2))

t3 <- bind_rows(t1, t2, .id = "Study") %>% arrange(scorecomp)

par(mar=c(5,4,2,2))
forest(t3$estimate, ci.lb = t3$conf.low, ci.ub = t3$conf.high,
       refline = 1, xlab = "Odds ratio (per unit increase in score)", pch = 18, 
       transf = exp, psize = 1.5,
       slab = rep(scorecomp, 2)) 




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


