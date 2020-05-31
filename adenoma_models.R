source("adenoma_crc.R")

# Logistic regression to find discriminants
mod <- function(x) glm(ct ~ x + country + sex + age, data = adenoma.meta)
mod2 <- function(x) glm(ct ~ x + country + sex + age, data = crc.meta)

fits <- apply(adenoma, 2, mod)
fits2 <- apply(crc, 2, mod2)

library(broom)
mods.adenoma <- map_df(fits, tidy) %>% filter(term == "x") %>% add_column(compound = colnames(mat2))
mods.crc <- map_df(fits2, tidy) %>% filter(term == "x") %>% add_column(compound = colnames(mat2))

library(ggplot2)
ggplot() + geom_point(data=mods.adenoma, aes(compound, log10(p.value)), colour = "red") +
  geom_point(data=mods.crc, aes(compound, log10(p.value)), colour = "blue") +
  theme_minimal() + scale_y_reverse() +
  theme(axis.text.x = element_text(angle = 90, hjust = 1))

# Model for transformation to residuals
adj   <- function(x) residuals(lm(x ~ country, data = adenoma.meta))
adj2   <- function(x) residuals(lm(x ~ country, data = crc.meta))
adjmat1 <- apply(adenoma, 2, adj)
adjmat2 <- apply(crc, 2, adj2)

# Bind case-control status to adjusted matrix
plsdat1 <- cbind(path.group = adenoma.meta$path.group, adjmat1) %>% as.tibble
plsdat2 <- cbind(path.group = crc.meta$path.group, adjmat2) %>% as.tibble

# Start with a sensible number of components eg 10

library(caret)
set.seed(100)
myfolds <- createMultiFolds(plsdat2$path.group, k = 5, times = 5)
control <- trainControl("repeatedcv", index = myfolds, selectionFunction = "oneSE")
# Fit models over different tuning parameters
mod1 <- train(path.group ~ ., data = plsdat2, method = "pls", metric = "Accuracy", 
              tuneLength = 20, trControl = control)

# Get components with lowest RMSEP, "one SE" and "permutation" methods
print(paste("Lowest RMSEP from", which.min(cv$val[estimate = "adjCV", , ]) - 1, "comps"))
print(paste("one SE method suggests", selectNcomp(mod, method = "onesigma", plot = F), "comps"))
permut     <- selectNcomp(mod, method = "randomization", plot = T)

print(paste("Lowest RMSEP from", best.dims, "comp(s);", 
            "one SE method suggests", onesigma, "comp(s);",
            "permutation method suggests", permut, "comp(s)"))
