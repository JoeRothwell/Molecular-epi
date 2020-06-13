# Get two CRC case controls (code repeated from CRC_prep_data.R)
library(tidyverse)
library(haven)
source("adenoma_crc.R")

# Remove duplicated Idepics (with dplyr or base). Also get follow up time and colorectal site
var.list <- c("Country", "Center", "Sex", "Match_Caseset", "L_School", #"Smoke_Int", 
              "Smoke_Stat", "Smoke_Intensity", "Fasting_C", "Menopause", "Phase_Mnscycle")

meta <- read_dta("clrt_caco.dta") %>% 
  mutate(Tfollowup.days = D_Dgclrt - D_Bld_Coll, Tfollowup = Tfollowup.days/365.25, 
         location = case_when(
           Case_Mal_Colon_Prox == 1 ~ 1, Case_Mal_Colon_Dist == 1 ~ 2,
           Case_Mal_Colon_Nos  == 1 ~ 4, Case_Mal_Rectum     == 1 ~ 3)) %>%
  group_by(Match_Caseset) %>% fill(c(D_Dgclrt, location), .direction = "downup") %>% ungroup() %>%
  select(-Match_Caseset, -Cncr_Caco_Clrt) %>%
  distinct(Idepic, .keep_all = T)

crc1 <- read_sas("clrt_caco_metabo.sas7bdat") %>% filter(!is.na(Aminoacid_Glu)) %>%
  left_join(meta, by = "Idepic", suffix = c("_1", "")) %>%
  mutate_at(vars(var.list), as.factor) %>% 
  mutate(Smoke_Int = fct_collapse(Smoke_Intensity, Other = c("8", "9", "10"))) %>% 
  filter(Country != 6)

library(lubridate)
crc2 <- read_csv("biocrates_p150.csv") %>% 
  select(Match_Caseset, Cncr_Caco_Clrt, ends_with("Idepic"), 
         matches("(carn|oacid|genic|roph|ingo|Sugars)[_]"), -contains("tdq")) %>%
  inner_join(meta, by = "Idepic") %>% mutate_at(vars(var.list), as.factor) %>% 
  mutate(Smoke_Int = fct_collapse(Smoke_Intensity, Other = c("8", "9", "10"))) %>%
  filter(Country != 6)

# Get compound overlaps between adenoma (128) and crc1 and crc2
# First subset compounds only from whole data and remove zero cols
expr <- "(carn|oacid|genic|roph|ingo|Sugars)[_]"
crc1p <- crc1 %>% select(matches(expr), -contains("tdq")) %>% select_if(~ sum(., na.rm = T) != 0) 
crc2p <- crc2 %>% select(matches(expr), -contains("tdq")) %>% select_if(~ sum(., na.rm = T) != 0)

# Get crc1 and crc2 overlap datasets, replace zeros with half min value, log and scale
overlap1 <- intersect(colnames(adenoma), colnames(crc1p))
overlap2 <- intersect(colnames(adenoma), colnames(crc2p))
hm <- function(x) min(x)/2
crc1sort <- crc1p[, overlap1] %>% na_if(0) %>% na.aggregate(FUN = hm) %>% log2 %>% scale
crc2sort <- crc2p[, overlap2] %>% na_if(0) %>% na.aggregate(FUN = hm) %>% log2 %>% scale

# Refit PLS models with adenoma and crc overlap dataset
# Bind case-control status to  matrix
# Adenoma (only 1 compound less)
plsdat1a <- data.frame(adenoma[, overlap1])
plsdat1b <- data.frame(adenoma[, overlap2])
plsdat1a$path.group <- as.factor(adenoma.meta$path.group)
plsdat1b$path.group <- as.factor(adenoma.meta$path.group)

library(caret)
# Overlap 1 for crc1

set.seed(111)
folds <- createMultiFolds(y = plsdat1a$path.group, k = 5, times = 5)
control <- trainControl("repeatedcv", index = folds, selectionFunction = "oneSE")
print(sapply(folds, length))

# Train PLS model
mod2 <- train(path.group ~ ., data = plsdat1a, method = "pls", metric = "Accuracy", 
              trControl = control, tuneLength = 20)
plot(mod2, main = paste("Model", length(mod2$coefnames), "compounds", sep = " "))
confusionMatrix(mod2)

predict.crc1 <- predict(mod2, newdata = crc1sort)
table(predict.crc1)
# 140 predicted adenomas, 801 predicted normal

# Overlap 2 for crc2
set.seed(111)
folds <- createMultiFolds(y = plsdat1b$path.group, k = 5, times = 5)
control <- trainControl("repeatedcv", index = folds, selectionFunction = "oneSE")
print(sapply(folds, length))

# Train PLS model
mod3 <- train(path.group ~ ., data = plsdat1b, method = "pls", metric = "Accuracy", 
              trControl = control, tuneLength = 20)
plot(mod3, main = paste("Model", length(mod3$coefnames), "compounds", sep = " "))
confusionMatrix(mod3)

predict.crc2 <- predict(mod3, newdata = crc2sort)
table(predict.crc2)
# 368 predicted adenomas, 1888 predicted normal



# CRC (13 compounds less)
plsdat2a <- data.frame(crc[, overlap1])
plsdat2b <- data.frame(crc[, overlap2])
plsdat2a$path.group <- as.factor(crc.meta$path.group)
plsdat2b$path.group <- as.factor(crc.meta$path.group)

# Overlap 1 for crc1

set.seed(111)
folds <- createMultiFolds(y = plsdat2a$path.group, k = 5, times = 5)
control <- trainControl("repeatedcv", index = folds, selectionFunction = "oneSE")
print(sapply(folds, length))

# Train PLS model
mod4 <- train(path.group ~ ., data = plsdat2a, method = "pls", metric = "Accuracy", 
              trControl = control, tuneLength = 20)
plot(mod4, main = paste("Model", length(mod4$coefnames), "compounds", sep = " "))
confusionMatrix(mod4)

predict.crc3 <- predict(mod4, newdata = crc1sort)
table(predict.crc3)
# 651 predicted crc, 290 predicted normal

# Overlap 2 for crc2

set.seed(111)
folds <- createMultiFolds(y = plsdat2b$path.group, k = 5, times = 5)
control <- trainControl("repeatedcv", index = folds, selectionFunction = "oneSE")
print(sapply(folds, length))

# Train PLS model
mod5 <- train(path.group ~ ., data = plsdat2b, method = "pls", metric = "Accuracy", 
              trControl = control, tuneLength = 20)
plot(mod5, main = paste("Model", length(mod5$coefnames), "compounds", sep = " "))
confusionMatrix(mod5)

predict.crc4 <- predict(mod5, newdata = crc2sort)
table(predict.crc4)
# 1567 predicted crc, 715 predicted normal


# Compare predictions with CRC case-control status
s1 <- cbind(crc1, pred.adenoma = predict.crc1, pred.crc = predict.crc3) %>% filter(Tfollowup < 2)
s2 <- cbind(crc2, pred.adenoma = predict.crc2, pred.crc = predict.crc4) %>% filter(Tfollowup < 2)

table(s1$Cncr_Caco_Clrt, s1$pred.adenoma)
table(s1$Cncr_Caco_Clrt, s1$pred.crc)
table(s2$Cncr_Caco_Clrt, s2$pred.adenoma)
table(s2$Cncr_Caco_Clrt, s2$pred.crc)
