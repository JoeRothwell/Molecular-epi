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
adj   <- function(x) residuals(lm(x ~ country, data = mat1))