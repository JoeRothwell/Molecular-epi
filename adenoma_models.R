source("adenoma_crc.R")

# Logistic regression to find discriminants
mod1 <- function(x) glm(ct ~ x + country + sex + age + batch, data = adenoma.meta)
mod2 <- function(x) glm(ct ~ x + country + sex + age + batch, data = crc.meta)
mod3 <- function(x) glm(ct ~ x + country + sex + age + batch, data = polyp.meta)

fits1 <- apply(adenoma, 2, mod1)
fits2 <- apply(crc, 2, mod2)
fits3 <- apply(polyp, 2, mod3)

library(broom)
mods.adenoma <- map_df(fits1, tidy) %>% filter(term == "x") %>%
  mutate(p.fdr = p.adjust(p.value, method = "fdr")) %>% add_column(compound = colnames(mat2)) %>%
  select(compound, estimate:p.fdr) #%>% arrange(p.fdr) %>% filter(p.fdr < 0.05)
# Copy and paste console output into Excel, then into powerpoint

mods.crc <- map_df(fits2, tidy) %>% filter(term == "x") %>% 
  mutate(p.fdr = p.adjust(p.value, method = "fdr")) %>% add_column(compound = colnames(mat2)) %>%
  select(compound, estimate:p.fdr) #%>% arrange(p.fdr) %>% filter(p.fdr < 0.05)

mods.polyp <- map_df(fits3, tidy) %>% filter(term == "x") %>% 
  mutate(p.fdr = p.adjust(p.value, method = "fdr")) %>% add_column(compound = colnames(mat2)) %>%
  select(compound, estimate:p.fdr) #%>% arrange(p.fdr) %>% filter(p.fdr < 0.05)

library(ggplot2)
ggplot() + geom_point(data=mods.adenoma, aes(compound, log10(p.value)), colour = "red") +
  geom_point(data=mods.crc, aes(compound, log10(p.value)), colour = "blue") +
  theme_grey() + scale_y_reverse() + ylab("-log10 raw p-value") + xlab("Metabolite") +
  geom_hline(yintercept = log10(0.00424), linetype = "dotted") +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5))

# Model for transformation to residuals
# Warning, better predictions without this adjustment for country!!!
adj1   <- function(x) residuals(lm(x ~ batch, data = adenoma.meta))
adj2   <- function(x) residuals(lm(x ~ batch, data = crc.meta))
adj3   <- function(x) residuals(lm(x ~ batch, data = polyp.meta))
adjmat1 <- apply(adenoma, 2, adj1)
adjmat2 <- apply(crc, 2, adj2)
adjmat3 <- apply(polyp, 2, adj3)

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
set.seed(111)
inTrain <- createDataPartition(y = plsdat1$path.group, p = 0.75, list = F)
training <- plsdat1[inTrain, ]
testing <- plsdat1[-inTrain, ]

# Create folds and training parameters
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
set.seed(111)
inTrain <- createDataPartition(y = plsdat2$path.group, p = 0.75, list = F)
training <- plsdat2[inTrain, ]
testing <- plsdat2[-inTrain, ]

# Create folds and training parameters
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
set.seed(111)
inTrain <- createDataPartition(y = plsdat3$path.group, p = 0.75, list = F)
training <- plsdat3[inTrain, ]
testing <- plsdat3[-inTrain, ]

# Create folds and training parameters
folds <- createMultiFolds(y = training$path.group, k = 5, times = 5)
control <- trainControl("repeatedcv", index = folds, selectionFunction = "oneSE")
print(sapply(folds, length))

# Train PLS model
mod3 <- train(path.group ~ ., data = training, method = "pls", metric = "Accuracy", 
              trControl = control, tuneLength = 20) 
plot(mod3)
confusionMatrix(mod3)
predictions3 <- predict(mod3, newdata = testing)
confusionMatrix(predictions3, reference = testing$path.group)



