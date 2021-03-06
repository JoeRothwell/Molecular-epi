# Get two CRC case controls (code repeated from CRC_prep_data.R)
library(tidyverse)
library(haven)
source("adenoma_crc.R")
source("CRC_prep_data.R")

# Remove duplicated Idepics (with dplyr or base). Also get follow up time and colorectal site
var.list <- c("Country", "Center", "Sex", "Match_Caseset", "L_School", #"Smoke_Int", 
              "Smoke_Stat", "Smoke_Intensity", "Fasting_C", "Menopause", "Phase_Mnscycle")

# Metadata (removed location for this study)
meta <- read_dta("clrt_caco.dta") %>% 
  group_by(Match_Caseset) %>% fill(D_Dgclrt, .direction = "downup") %>% ungroup() %>%
  mutate(Tfollowup.days = D_Dgclrt - D_Bld_Coll, Tfollowup = Tfollowup.days/365.25, 
         location = case_when(
           Case_Mal_Colon_Prox == 1 ~ 1, Case_Mal_Colon_Dist == 1 ~ 2,
           Case_Mal_Colon_Nos  == 1 ~ 4, Case_Mal_Rectum     == 1 ~ 3)) %>%
  #group_by(Match_Caseset) %>% fill(c(D_Dgclrt, location), .direction = "downup") %>% ungroup() %>%
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

# First subset compounds only from whole data and remove zero cols
crc.both <- bind_rows(crc1, crc2, .id = "lab")
expr <- "(carn|oacid|genic|roph|ingo|Sugars)[_]"
crcp <- crc.both %>% select(matches(expr), -contains("tdq")) %>% select_if(~ sum(., na.rm = T) != 0) 


# Get compound overlaps between adenoma (128) and crc case-control
# Get overlap dataset, replace zeros with half min value
overlap <- intersect(colnames(mat2a), colnames(crcp))
hm <- function(x) min(x)/2
crc.sort <- crcp[, overlap] %>% na_if(0) %>% na.aggregate(FUN = hm)

# Put cross-sectional and case-control samples together and make group labels
allmat <- rbind(crc.sort, mat2a) %>% log2 %>% scale
grps <- as.factor(c(rep("case-ctrlA", nrow(crc1)), rep("case-ctrlB", nrow(crc2)), rep("Hospital", nrow(mat2))))
grps1 <- as.factor(c(rep("case-ctrlA", nrow(crc1)), rep("case-ctrlB", nrow(crc2)), mat$path.group))


# Plot PCA
pca <- prcomp(allmat, scale. = F)
library(pca3d)
pca2d(pca, group = fct_inorder(grps1), legend = "topright")
box(which = "plot", lty = "solid")

# Adjust matrix with residuals method and repeat PCA
adjmat <- apply(allmat, 2, function(x) residuals(lm(x ~ grps)))
pca1 <- prcomp(adjmat, scale. = F)
pca2d(pca1, group = fct_inorder(grps1), legend = "topleft")
box(which = "plot", lty = "solid")

# Refit PLS models with adenoma and crc overlap dataset
# Bind case-control status to  matrix


# Make PLS data
# Adenoma (only 1 compound less)
#plsdat1 <- data.frame(adenoma[, overlap])
#plsdat1$path.group <- as.factor(adenoma.meta$path.group)

# Update: make PLS data from all data scaled together (? to check)
plsdat1 <- adjmat[grps1 %in% c("adenoma", "normal"), ] %>% data.frame()
plsdat1$path.group <- as.factor(adenoma.meta$path.group)

library(caret)
# Overlapping metabolites

set.seed(111)
folds <- createMultiFolds(y = plsdat1$path.group, k = 5, times = 5)
control <- trainControl("repeatedcv", index = folds, selectionFunction = "oneSE")
print(sapply(folds, length))

# Train PLS model
mod1 <- train(path.group ~ ., data = plsdat1, method = "pls", metric = "Accuracy", 
              trControl = control, tuneLength = 20)
plot(mod1, main = paste("Model", length(mod1$coefnames), "compounds", sep = " "))
confusionMatrix(mod1)

crc.sort <- adjmat[grps1 %in% c("case-ctrlA", "case-ctrlB"), ] %>% data.frame
predict.crc1 <- predict(mod1, newdata = crc.sort)
table(predict.crc1)
# 364 predicted adenomas, 2859 predicted normal (scaled separately)




# CRC (13 compounds less)
plsdat2 <- adjmat[grps1 %in% c("crc", "normal"), ] %>% data.frame()
plsdat2$path.group <- as.factor(crc.meta$path.group)

# Overlapping metabolites
set.seed(111)
folds <- createMultiFolds(y = plsdat2$path.group, k = 5, times = 5)
control <- trainControl("repeatedcv", index = folds, selectionFunction = "oneSE")
print(sapply(folds, length))

# Train PLS model
mod2 <- train(path.group ~ ., data = plsdat2, method = "pls", metric = "Accuracy", 
              trControl = control, tuneLength = 20)
plot(mod2, main = paste("Model", length(mod2$coefnames), "compounds", sep = " "))
confusionMatrix(mod2)

predict.crc2 <- predict(mod2, newdata = crc.sort)
table(predict.crc2)
# 1372 predicted crc, 1851 predicted normal





# Polyp (13 compounds less)
plsdat3 <- adjmat[grps1 %in% c("polyp", "normal"), ] %>% data.frame()
plsdat3$path.group <- as.factor(polyp.meta$path.group)

# Overlapping metabolites
set.seed(111)
folds <- createMultiFolds(y = plsdat3$path.group, k = 5, times = 5)
control <- trainControl("repeatedcv", index = folds, selectionFunction = "oneSE")
print(sapply(folds, length))

# Train PLS model
mod3 <- train(path.group ~ ., data = plsdat3, method = "pls", metric = "Accuracy", 
              trControl = control, tuneLength = 20)
plot(mod3, main = paste("Model", length(mod3$coefnames), "compounds", sep = " "))
confusionMatrix(mod3)

predict.crc3 <- predict(mod3, newdata = crc.sort)
table(predict.crc3)
# 332 predicted polyp, 2891 predicted normal



# Compare predictions with CRC case-control status
s1 <- cbind(crc.both, pred.adenoma = predict.crc1, pred.crc = predict.crc2, 
            pred.polyp = predict.crc3)
s2 <- s1 %>% filter(Tfollowup < 5)
s3 <- s1 %>% filter(Tfollowup < 2)
s4 <- s1 %>% filter(Age_Blood < 60)
s5 <- s1 %>% filter(Age_Blood > 60)

# Adenoma, all case-control subjects and <5 and <2 y follow up only
table(s1$Cncr_Caco_Clrt, s1$pred.adenoma)
table(s2$Cncr_Caco_Clrt, s2$pred.adenoma)
table(s3$Cncr_Caco_Clrt, s3$pred.adenoma)

# CRC
table(s1$Cncr_Caco_Clrt, s1$pred.crc)
table(s2$Cncr_Caco_Clrt, s2$pred.crc)
table(s3$Cncr_Caco_Clrt, s3$pred.crc)

# Polyp
table(s1$Cncr_Caco_Clrt, s1$pred.polyp)
table(s2$Cncr_Caco_Clrt, s2$pred.polyp)
table(s3$Cncr_Caco_Clrt, s3$pred.polyp)

table(s4$Cncr_Caco_Clrt, s4$pred.polyp)
table(s5$Cncr_Caco_Clrt, s5$pred.polyp)

library(sjmisc)
s1$Age_cat <- dicho(s1$Age_Blood)
table(s1$Age_cat, col = s1$pred.polyp)



