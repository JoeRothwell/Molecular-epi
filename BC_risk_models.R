# BC risk models for metabolites
library(tidyverse)
library(readxl)

# Read 1623 observations of 44 intensity variables (appears to be final scaled data) and metadata
ints <- read_tsv("1507_XMetabolite_std_cpmg_E3N.txt")

# Lifestyle data. Subset variables needed
meta <- read_csv("Lifepath_meta.csv", na = "9999") %>%
  select(CODBMB, CT, MATCH, PLACE, AGE, BMI, BP, RTH, ALCOHOL, MENOPAUSE, FASTING, SMK, DIABETE, CENTTIMECat1, 
         CENTTIME, SAMPYEAR, STOCKTIME, DURTHSDIAG) %>%
  mutate_at(vars(-CODBMB, -CT, -AGE, -BMI, -STOCKTIME, -RTH, -ALCOHOL, -CENTTIME, -DURTHSDIAG), as.factor)

# Metabolite risk models --------------------------------------------------

# CLR models to get odds ratios for metabolites
dat <- left_join(meta, ints, by = "CODBMB")
  
# Run models for all, pre-menopausal only and post-menopausal only
#base <- CT ~ BMI + SMK + DIABETE + RTH + ALCOHOL + DURTHSDIAG + CENTTIME + STOCKTIME + strata(MATCH)

Ints <- dat %>% select(`3Hydroxybutyrate`:Succinate) %>% as.matrix

library(survival)
fits0 <- apply(Ints, 2, function(x) clogit(CT ~ BMI + SMK + DIABETE + RTH + ALCOHOL + 
          DURTHSDIAG + CENTTIME + STOCKTIME + strata(MATCH) + x, data = dat))
fits1 <- apply(Ints, 2, function(x) clogit(CT ~ BMI + SMK + DIABETE + RTH + ALCOHOL + 
          DURTHSDIAG + CENTTIME + STOCKTIME + strata(MATCH) + x, data = dat, subset = MENOPAUSE == 0))
fits2 <- apply(Ints, 2, function(x) clogit(CT ~ BMI + SMK + DIABETE + RTH + ALCOHOL + 
          DURTHSDIAG + CENTTIME + STOCKTIME + strata(MATCH) + x, data = dat, subset = MENOPAUSE == 1))

cmpd_meta <- read.csv("NMR_cmpd_metadata.csv")

# Function to tidy and present output table

tidy.output <- function(fits) {
  
  library(broom) 
  
  # Function to exponentiate and round
  round.exp <- function(x) round(exp(x), 2)
  df <- map_df(fits, tidy) %>% filter(term == "x") %>%
    mutate_at(.vars = c(OR = "estimate", "conf.low", "conf.high"), .funs = round.exp) %>%
    cbind(Compound = names(fits)) %>%
    left_join(cmpd_meta, by  = "Compound")
  
  # Make columns
  df$P.value <- round(df$p.value, 3)
  df$FDR <- round(p.adjust(df$p.value, method = "fdr"), 3)
  df$Bonferroni <- round(p.adjust(df$p.value, method = "bonferroni"), 3)
  df$CI.95 <- paste("(", df$conf.low, ", ", df$conf.high, ")", sep = "")

  # Select columns and order
  df1 <- df %>% 
    select(Compound = "display_name", "description", "OR", "CI.95", "P.value", "FDR", "Bonferroni") %>%
    arrange(description)
  
}

all <- tidy.output(fits0)
pre <- tidy.output(fits1)
post <- tidy.output(fits2)

# Retain only metabolite groups with at least one p-value < 0.05
tab <- bind_rows("All" = all, "Pre" = pre, "Post" = post, .id = "Analysis") %>%
  group_by(Compound) %>% filter(min(P.value) < 0.05) %>% 
  select(Compound, description, everything()) %>%
  arrange(description, Compound) %>% as.data.frame

library(broom)
t2 <- map_df(fits0, tidy) %>% filter(term == "x") %>% bind_cols(cmpd_meta) %>%
  arrange(description)

library(stargazer)
stargazer(tab, summary = F, type = "html", out = "metabolite_table_selected.html")

# Plot data with Metafor
rowvec <- rev(c(1, 3, 5:7, 9:17, 19, 21, 23:24, 26:28, 30:31, 33:47, 49:51, 53:55))
rowvec2 <- cumsum(as.numeric(t2$description))
tabulate(t2$description)

par(mar=c(5,4,1,2))
library(metafor)
forest(t2$estimate, ci.lb = t2$conf.low, ci.ub = t2$conf.high, refline = 1, #xlab = xtitle, 
       xlab = "Multivariable adjusted odds ratio", ylim = c(1, 60), 
       rows = rowvec, efac=0,
       transf = exp, #pch = 18, 
       col = "grey",
       psize = 1, slab = t2$display_name)
hh <- par("usr")
text(hh[1], max(rowvec) + 4, "Compound", pos = 4)
text(hh[2], max(rowvec) + 4, "OR [95% CI]", pos = 2)


# Funnel plots for metabolites
funnel(x = t2$estimate, sei = t2$std.error)
funnel(x = t2$estimate, sei = t2$std.error, yaxis = "vi")

# Investigation of Ethanol
data <- left_join(meta2, ints, by = "CODBMB")
boxplot(data$Ethanol ~ data$CT + data$MENOPAUSE, varwidth = T, outline = F,
        names = c("Control, pre", "Case, pre", "Control, post", "Case, post"),
        col = "dodgerblue", ylab = "Plasma ethanol conc (scaled)")


# ---------------------------------------------------------------------------------------------------

# Conditional logistic regression to get odds ratios for lifestyle factors
# Same co-variates as in original manuscript
library(survival)
fit <- clogit(CT ~ scale(BMI) + SMK + DIABETE + #BP + 
                scale(RTH) + scale(ALCOHOL) + scale(DURTHSDIAG) + 
                scale(CENTTIME) + STOCKTIME + strata(MATCH), data = meta, subset = MENOPAUSE == 1) 
# output <- cbind(exp(coef(fit)), exp(confint(fit)))

library(broom)
t1 <- tidy(fit) %>% select(-(std.error:p.value))

library(metafor)
#dev.off()
par(mar=c(5,4,2,2))
forest(t1$estimate, ci.lb = t1$conf.low, ci.ub = t1$conf.high, refline = 1, 
       xlab = "Multivariable adjusted odds ratio",
       alim = c(0, 5),
       xlim = c(-5, 10),
       transf = exp, pch = 18, psize = 1.5, slab = t1$term,
       main = paste(fit$nevent, "case-control pairs"))
#xlim = c(-1, 3)
hh <- par("usr")
text(hh[1], nrow(t1) + 2, "Variable", pos = 4)
text(hh[2], nrow(t1) + 2, "OR [95% CI]", pos = 2)
# matching factors removed!

# Non-metabolite model lme4
library(lme4)
fit1 <- glmer(CT ~ BMI + SMK + DIABETE + (1|MATCH), data = meta1, family = "binomial")
summary(fit1)
