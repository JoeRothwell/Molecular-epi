source("adenoma_crc.R")

# Logistic regression to find discriminants. Make functions for models

mod1 <- function(x, dat) glm(ct ~ x + country + sex + age + batch, data = dat)
# CRC with only CR subjects adjusted for bmi
mod2 <- function(x, dat) glm(ct ~ x + bmi + diabetes + sex + age + smoke + 
                               alcohol_drinks_week + batch, data = dat)

library(broom)
# Adenoma both countries
fits1 <- apply(mat2a[adenoma, ], 2, mod1, dat1[adenoma, ]) %>% map_df(tidy) %>% 
  filter(term == "x") %>% mutate(p.fdr = p.adjust(p.value, method = "fdr")) %>% 
  add_column(compound = colnames(mat2)) %>% select(compound, estimate:p.fdr)

# CRC all samples (no BMI adjustment etc)
fits2 <- apply(mat2a[crc, ], 2, mod1, dat1[crc, ]) %>% map_df(tidy) %>% 
  filter(term == "x") %>% mutate(p.fdr = p.adjust(p.value, method = "fdr")) %>% 
  add_column(compound = colnames(mat2)) %>% select(compound, estimate:p.fdr)

# CRC CR samples only (adjusted BMI and other covariates)
fits2a <- apply(mat2a[crc.cr, ], 2, mod2, dat1[crc.cr, ]) %>% map_df(tidy) %>% 
  filter(term == "x") %>% mutate(p.fdr = p.adjust(p.value, method = "fdr")) %>% 
  add_column(compound = colnames(mat2)) %>% select(compound, estimate:p.fdr)

# Polyp both countries 
fits3 <- apply(mat2a[polyp, ], 2, mod1, dat1[polyp, ]) %>% map_df(tidy) %>% 
  filter(term == "x") %>% mutate(p.fdr = p.adjust(p.value, method = "fdr")) %>% 
  add_column(compound = colnames(mat2)) %>% select(compound, estimate:p.fdr)

# Write to clipboard for pasting results into Excel
library(clipr)
write_clip(fits2a)

# Manhattan plot overlaying adenoma and CRC results 
library(ggplot2)
ggplot() + geom_point(data=fits1, aes(compound, log10(p.value)), colour = "red") +
  geom_point(data=fits2a, aes(compound, log10(p.value)), colour = "blue") +
  theme_grey() + scale_y_reverse() + ylab("-log10 raw p-value") + xlab("Metabolite") +
  geom_hline(yintercept = log10(0.00424), linetype = "dotted") +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5))

# Model function for transformation to residuals
# Warning, better predictions without this adjustment for country!!!
adj1   <- function(x, dat) residuals(lm(x ~ batch, data = dat))

adjmat1 <- apply(mat2a[adenoma, ], 2, adj1, dat1[adenoma, ])
adjmat2 <- apply(mat2a[crc, ], 2, adj1, dat1[crc, ])
adjmat3 <- apply(mat2a[polyp, ], 2, adj1, dat1[polyp, ])

# Bind case-control status to adjusted (or unadjusted) matrix
plsdat1 <- bind_cols(path.group = as.factor(dat1$ct), mat2a)[adenoma, ]
plsdat2 <- bind_cols(path.group = as.factor(dat1$ct), mat2a)[crc, ]
plsdat3 <- bind_cols(path.group = as.factor(dat1$ct), mat2a)[polyp, ]

# Predictive model by PLS
# Start with a sensible number of components eg 10

# Adenoma. First split into training and test sets
set.seed(111)
library(caret)
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



