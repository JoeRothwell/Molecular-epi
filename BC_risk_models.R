# BC risk models

library(tidyverse)
library(readxl)

# Read 1623 observations of 44 intensity variables (appears to be final scaled data) and metadata
ints <- read_tsv("1507_XMetabolite_std_cpmg_E3N.txt")
meta <- read.csv("Lifepath_meta.csv")

# Lifestyle data only
# Breast cancer risk model. Subset variables needed
meta1 <- meta %>%
  select(CODBMB, CT, MATCH, PLACE, AGE, BMI, BP, RTH, ALCOHOL, MENOPAUSE, FASTING, SMK, DIABETE, CENTTIMECat1, 
         CENTTIME, SAMPYEAR, STOCKTIME, DURTHSDIAG) %>%
  mutate_at(vars(-CODBMB, -CT, -AGE, -BMI, -STOCKTIME, -RTH, -ALCOHOL, -CENTTIME, -DURTHSDIAG), as.factor)

# Replace 9999 with NA (for just numeric or all columns)
meta2 <- meta1 %>% mutate_if(is.numeric, list( ~ na_if(., 9999))) %>% mutate(BP = na_if(BP, 9999))

# Subset pre or post menopausal
meta2.pre  <- meta2 %>% filter(MENOPAUSE == 0)
meta2.post <- meta2 %>% filter(MENOPAUSE == 1)

# Conditional logistic regression to get odds ratios for lifestyle factors
# Same co-variates as in original manuscript
library(survival)
fit <- clogit(CT ~ scale(BMI) + SMK + DIABETE + #BP + 
                scale(RTH) + scale(ALCOHOL) + scale(DURTHSDIAG) + 
                scale(CENTTIME) + STOCKTIME + strata(MATCH), data = meta2.pre) 
# output <- cbind(exp(coef(fit)), exp(confint(fit)))

library(broom)
t1 <- tidy(fit) %>% select(-(std.error:p.value))

library(metafor)
dev.off()
par(mar=c(5,4,1,2))
forest(t1$estimate, ci.lb = t1$conf.low, ci.ub = t1$conf.high, refline = 1, 
       xlab = "Multivariable adjusted odds ratio",
       transf = exp, pch = 18, psize = 1.5, slab = t1$term)  #, 
#alim = c(0.5,2.5), )
#xlim = c(-1, 3)
hh <- par("usr")
text(hh[1], nrow(t1) + 2, "Variable", pos = 4)
text(hh[2], nrow(t1) + 2, "OR [95% CI]", pos = 2)
# matching factors removed!

# Metabolite risk models --------------------------------------------------

# CLR models to get odds ratios for metabolites

data <- left_join(meta2, ints, by = "CODBMB")

clr <- function(x) { 
  clogit(CT ~ x + BMI + SMK + DIABETE + #BP + 
           RTH + ALCOHOL + DURTHSDIAG + 
           CENTTIME + STOCKTIME + strata(MATCH), data = meta2)
}

metabs <- data %>% select(`3Hydroxybutyrate`:Succinate) %>% as.matrix
multifit <- apply(metabs, 2, clr)
t2 <- map_df(multifit, tidy) %>% filter(term == "x")

dev.off()
par(mar=c(5,4,1,2))
library(metafor)
forest(t2$estimate, ci.lb = t2$conf.low, ci.ub = t2$conf.high, refline = 1, #xlab = xtitle, 
       xlab = "Multivariable adjusted odds ratio",
       transf = exp, pch = 18, psize = 1, slab = names(multifit))
hh <- par("usr")
text(hh[1], nrow(t2) + 2, "Compound", pos = 4)
text(hh[2], nrow(t2) + 2, "OR [95% CI]", pos = 2)
#alim = c(0,2), xlim = c(-1, 3)) 


# Funnel plots for metabolites
funnel(x = t2$estimate, sei = t2$std.error)
funnel(x = t2$estimate, sei = t2$std.error, yaxis = "vi")

# Investigation of Ethanol
boxplot(data$Ethanol ~ data$CT + data$MENOPAUSE, varwidth = T, outline = F,
        names = c("Control, pre", "Case, pre", "Control, post", "Case, post"),
        col = "dodgerblue", ylab = "Plasma ethanol conc (scaled)")