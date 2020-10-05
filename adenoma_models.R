source("adenoma_crc.R")

# Logistic regression to find discriminants
mod <- function(x) glm(ct ~ x + country + sex + age, data = adenoma.meta)
mod2 <- function(x) glm(ct ~ x + country + sex + age, data = crc.meta)
mod1 <- function(x) glm(ct ~ x + country + sex + age, data = polyp.meta)

fits <- apply(adenoma, 2, mod)
fits2 <- apply(crc, 2, mod2)
fits1 <- apply(polyp, 2, mod1)

library(broom)
mods.adenoma <- map_df(fits, tidy) %>% filter(term == "x") %>%
  mutate(p.fdr = p.adjust(p.value, method = "fdr")) %>% add_column(compound = colnames(mat2))
mods.crc <- map_df(fits2, tidy) %>% filter(term == "x") %>% 
  mutate(p.fdr = p.adjust(p.value, method = "fdr")) %>% add_column(compound = colnames(mat2))
mods.polyp <- map_df(fits1, tidy) %>% filter(term == "x") %>% 
  mutate(p.fdr = p.adjust(p.value, method = "fdr")) %>% add_column(compound = colnames(mat2))

library(ggplot2)
ggplot() + geom_point(data=mods.adenoma, aes(compound, log10(p.value)), colour = "red") +
  geom_point(data=mods.crc, aes(compound, log10(p.value)), colour = "blue") +
  theme_grey() + scale_y_reverse() + ylab("-log10 raw p-value") + xlab("Metabolite") +
  geom_hline(yintercept = log10(0.00424), linetype = "dotted") +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5))

# Model for transformation to residuals
# Warning, better predictions without this adjustment for country!!!
adj   <- function(x) residuals(lm(x ~ country, data = adenoma.meta))
adj2   <- function(x) residuals(lm(x ~ country, data = crc.meta))
adjmat1 <- apply(adenoma, 2, adj)
adjmat2 <- apply(crc, 2, adj2)

# Bind case-control status to adjusted (or unadjusted) matrix
plsdat1 <- data.frame(adenoma)
plsdat2 <- data.frame(crc)
plsdat3 <- data.frame(polyp)
plsdat1$path.group <- as.factor(adenoma.meta$path.group)
plsdat2$path.group <- as.factor(crc.meta$path.group)
plsdat3$path.group <- as.factor(polyp.meta$path.group)

# Predictive model by PLS
# Start with a sensible number of components eg 10
library(caret)
# Adenoma
# Split into training and test sets
inTrain <- createDataPartition(y = plsdat1$path.group, p = 0.75, list = F)
training <- plsdat1[inTrain, ]
testing <- plsdat1[-inTrain, ]

set.seed(111)
folds <- createMultiFolds(y = training$path.group, k = 5, times = 5)
control <- trainControl("repeatedcv", index = folds, selectionFunction = "oneSE")
print(sapply(folds, length))

# Train PLS model
mod0 <- train(path.group ~ ., data = training, method = "pls", metric = "Accuracy", 
              trControl = control, tuneLength = 10)

plot(mod0)
confusionMatrix(mod0)
predictions0 <- predict(mod0, newdata = testing)
confusionMatrix(predictions0, reference = testing$path.group)


# CRC
# Split into training and test sets
inTrain <- createDataPartition(y = plsdat2$path.group, p = 0.75, list = F)
training <- plsdat2[inTrain, ]
testing <- plsdat2[-inTrain, ]

set.seed(111)
folds <- createMultiFolds(y = training$path.group, k = 5, times = 5)
control <- trainControl("repeatedcv", index = folds, selectionFunction = "oneSE")
print(sapply(folds, length))

# Train PLS model
mod1 <- train(path.group ~ ., data = training, method = "pls", metric = "Accuracy", 
              trControl = control, tuneLength = 20) 
plot(mod1)
confusionMatrix(mod1)
predictions1 <- predict(mod1, newdata = testing)
confusionMatrix(predictions1, reference = testing$path.group)


# Polyp
# Split into training and test sets
inTrain <- createDataPartition(y = plsdat3$path.group, p = 0.75, list = F)
training <- plsdat3[inTrain, ]
testing <- plsdat3[-inTrain, ]

set.seed(111)
folds <- createMultiFolds(y = training$path.group, k = 5, times = 5)
control <- trainControl("repeatedcv", index = folds, selectionFunction = "oneSE")
print(sapply(folds, length))

# Train PLS model
mod1 <- train(path.group ~ ., data = training, method = "pls", metric = "Accuracy", 
              trControl = control, tuneLength = 20) 
plot(mod1)
confusionMatrix(mod1)
predictions1 <- predict(mod1, newdata = testing)
confusionMatrix(predictions1, reference = testing$path.group)

