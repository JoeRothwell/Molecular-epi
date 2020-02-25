# BC risk models for metabolites and also for lifestyle variables alone
# Also heatmap of differences
source("BC_prep_data.R")

# CLR models to get odds ratios for metabolites: all, pre-menopausal only and post-menopausal only

# All subjects -------------------

fits0 <- apply(ints, 2, function(x) clogit(CT ~ BMI + SMK + DIABETE + RTH + ALCOHOL + 
         DURTHSBMB + CENTTIME + STOCKTIME + strata(MATCH) + x, data = meta))

# Analysis by quartiles of metabolite concentration (resist outliers). 
quartiles <- ints0 %>% mutate_all(funs(cut_number(., n = 4, labels = 1:4))) 

# Apply across all compounds
fits0a <- apply(quartiles, 2, function(x) {
  Q1Q4 <- x == 1 | x == 4
  clogit(CT ~ x[Q1Q4] + BMI + SMK + DIABETE + RTH + ALCOHOL + DURTHSBMB + CENTTIME + #RACK +
           STOCKTIME + strata(MATCH), data = meta[Q1Q4, ])
} )

# Tidy and plot
library(broom)
t1 <- map_df(fits0, tidy) %>% filter(str_detect(term, "x")) %>% bind_cols(cmpd.meta) %>% arrange(description)

par(mfrow = c(1,2))
par(mar=c(5,4,1,2))
library(metafor)
forest(t1$estimate, ci.lb = t1$conf.low, ci.ub = t1$conf.high, refline = 1,
       ylim = c(1, max(rowvec) + 3), xlab = "Odds ratio per SD increase in conc", 
       transf = exp, rows = rowvec, efac = 0.5, pch = 18, cex = 0.8, psize = 1.5, 
       annosym = c("  (", " to ", ")"), slab = t1$display_name)

hh <- par("usr")
text(hh[1], max(rowvec) + 2, "Metabolite", pos = 4, cex = 0.8)
text(hh[2], max(rowvec) + 2, "OR [95% CI]", pos = 2, cex = 0.8)

t2 <- map_df(fits0a, tidy) %>% filter(str_detect(term, "x")) %>% bind_cols(cmpd.meta) %>% arrange(description)

par(mar=c(5,4,1,2))
forest(t2$estimate, ci.lb = t2$conf.low, ci.ub = t2$conf.high, refline = 1,
       ylim = c(1, max(rowvec) + 3), xlab = "Odds ratio Q4 vs Q1", 
       transf = exp, rows = rowvec, efac = 0.5, pch = 18, cex = 0.8, psize = 1.5, 
       annosym = c("  (", " to ", ")"), slab = t2$display_name)

hh <- par("usr")
text(hh[1], max(rowvec) + 2, "Metabolite", pos = 4, cex = 0.8)
text(hh[2], max(rowvec) + 2, "OR [95% CI]", pos = 2, cex = 0.8)

# Smile plot
#plot(t1$statistic, -log10(t1$p.value))
library(ggplot2)
library(ggrepel)
library(scales)
# For slides, dimension 415 x 377
ggplot(t1, aes(exp(estimate), log10(p.value))) + geom_point() + theme_bw() +
  xlim(0.8, 1.2) +
  scale_y_reverse(breaks = c(-2, -1, 0), labels = function(x) 10^x) +
  xlab("Odds ratio per SD increase concentration") + ylab("P-value") +
  geom_text_repel(aes(label = display_name), size = 3,
            data = t1[t1$p.value < 0.3, ]) +
  geom_hline(yintercept = log10(0.05), linetype = "dashed") +
  geom_vline(xintercept = 1, linetype = "dashed") +
  ggtitle("All subjects")


# Pre-menopausal -------------------

quartiles <- ints0[pre, ] %>% mutate_all(funs(cut_number(., n = 4, labels = 1:4))) 
meta1 <- meta[pre, ]

# Note: need to remove hormone treatment therapy variable
fits1 <- apply(ints[pre, ], 2, function(x) clogit(CT ~ BMI + SMK + DIABETE + RTH + ALCOHOL + 
         #DURTHSBMB + 
           CENTTIME + STOCKTIME + strata(MATCH) + x, data = meta[pre, ]))

fits1a <- apply(quartiles, 2, function(x) {
  Q1Q4 <- x == 1 | x == 4
  clogit(CT ~ x[Q1Q4] + BMI + SMK + DIABETE + RTH + ALCOHOL + #DURTHSBMB + 
           CENTTIME + STOCKTIME + strata(MATCH), data = meta1[Q1Q4, ])
} )

library(broom)
t1 <- map_df(fits1, tidy) %>% filter(str_detect(term, "x")) %>% bind_cols(cmpd.meta) %>% arrange(description)

par(mfrow = c(1,2))
par(mar=c(5,4,1,2))
library(metafor)
forest(t1$estimate, ci.lb = t1$conf.low, ci.ub = t1$conf.high, refline = 1,
       ylim = c(1, max(rowvec) + 3), xlab = "Odds ratio per SD increase in conc", 
       transf = exp, rows = rowvec, efac = 0.5, pch = 18, cex = 0.8, psize = 1.5, 
       annosym = c("  (", " to ", ")"), slab = t1$display_name)

hh <- par("usr")
text(hh[1], max(rowvec) + 2, "Metabolite", pos = 4, cex = 0.8)
text(hh[2], max(rowvec) + 2, "OR [95% CI]", pos = 2, cex = 0.8)

t2 <- map_df(fits1a, tidy) %>% filter(str_detect(term, "x")) %>% bind_cols(cmpd.meta) %>% arrange(description)

par(mar=c(5,4,1,2))
forest(t2$estimate, ci.lb = t2$conf.low, ci.ub = t2$conf.high, refline = 1,
       ylim = c(1, max(rowvec) + 3), xlab = "Odds ratio Q4 vs Q1", 
       transf = exp, rows = rowvec, efac = 0.5, pch = 18, cex = 0.8, psize = 1.5, 
       annosym = c("  (", " to ", ")"), slab = t2$display_name)

hh <- par("usr")
text(hh[1], max(rowvec) + 2, "Metabolite", pos = 4, cex = 0.8)
text(hh[2], max(rowvec) + 2, "OR [95% CI]", pos = 2, cex = 0.8)

library(ggplot2)
library(ggrepel)
ggplot(t1, aes(exp(estimate), log10(p.value))) + geom_point() + theme_bw() +
  xlim(c(0, 2)) + 
  scale_y_reverse(breaks = c(-3, -2, -1, 0), labels = function(x) 10^x) +
  xlab("Odds ratio per SD increase concentration") + ylab("P-value") +
  geom_text_repel(aes(label = display_name), size = 3,
                  data = t1[t1$p.value < 0.04, ] ) +
  geom_hline(yintercept = c(log10(0.05), log10(0.014)), linetype = "dashed") +
  geom_vline(xintercept = 1, linetype = "dashed") +
  ggtitle("Pre-menopausal")

# Post-menopausal --------------------

quartiles <- ints0[post, ] %>% mutate_all(funs(cut_number(., n = 4, labels = 1:4))) 
meta2 <- meta[post, ]

fits2 <- apply(ints[post, ], 2, function(x) clogit(CT ~ BMI + SMK + DIABETE + RTH + ALCOHOL + 
          DURTHSBMB + CENTTIME + STOCKTIME + strata(MATCH) + x, data = meta2))

fits2a <- apply(quartiles, 2, function(x) {
  Q1Q4 <- x == 1 | x == 4
  clogit(CT ~  x[Q1Q4] + BMI + SMK + DIABETE + RTH + ALCOHOL + DURTHSBMB + CENTTIME + 
           STOCKTIME + strata(MATCH), data = meta2[Q1Q4, ])
} )

library(broom)
t1 <- map_df(fits2, tidy) %>% filter(str_detect(term, "x")) %>% bind_cols(cmpd.meta) %>% arrange(description)

par(mfrow = c(1,2))
par(mar=c(5,4,1,2))
library(metafor)
forest(t1$estimate, ci.lb = t1$conf.low, ci.ub = t1$conf.high, refline = 1,
       ylim = c(1, max(rowvec) + 3), xlab = "Odds ratio per SD increase in conc", 
       transf = exp, rows = rowvec, efac = 0.5, pch = 18, cex = 0.8, psize = 1.5, 
       annosym = c("  (", " to ", ")"), slab = t1$display_name)

hh <- par("usr")
text(hh[1], max(rowvec) + 2, "Metabolite", pos = 4, cex = 0.8)
text(hh[2], max(rowvec) + 2, "OR [95% CI]", pos = 2, cex = 0.8)

t2 <- map_df(fits2a, tidy) %>% filter(str_detect(term, "x")) %>% bind_cols(cmpd.meta) %>% arrange(description)

par(mar=c(5,4,1,2))
forest(t2$estimate, ci.lb = t2$conf.low, ci.ub = t2$conf.high, refline = 1,
       ylim = c(1, max(rowvec) + 3), xlab = "Odds ratio Q4 vs Q1", 
       transf = exp, rows = rowvec, efac = 0.5, pch = 18, cex = 0.8, psize = 1.5, 
       annosym = c("  (", " to ", ")"), slab = t2$display_name)

hh <- par("usr")
text(hh[1], max(rowvec) + 2, "Metabolite", pos = 4, cex = 0.8)
text(hh[2], max(rowvec) + 2, "OR [95% CI]", pos = 2, cex = 0.8)

library(ggplot2)
library(ggrepel)
ggplot(t1, aes(exp(estimate), log10(p.value))) + geom_point() + theme_bw() +
  xlim(c(0.85, 1.15)) +  
  scale_y_reverse(breaks = c(-2, -1, 0), labels = function(x) 10^x) +
  xlab("Odds ratio per SD increase concentration") + ylab("P-value") +
  geom_text_repel(aes(label = display_name), size = 3,
            data = t1[t1$p.value < 0.15, ]) +
  geom_hline(yintercept = log10(0.05), linetype = "dashed") +
  geom_vline(xintercept = 1, linetype = "dashed") +
  ggtitle("c) Post-menopausal")


# Funnel plots for metabolites
funnel(x = t2$estimate, sei = t2$std.error)
funnel(x = t2$estimate, sei = t2$std.error, yaxis = "vi")

# ---------------------------------------------------------------

# Tables for manuscript
# Generate tidy output table from models
tidy.output <- function(mod) {
  
  library(broom) 
  # Function to exponentiate and round
  round.exp <- function(x) round(exp(x), 2)
  df <- map_df(mod, tidy) %>% filter(term == "x") %>%
    mutate_at(.vars = c(OR = "estimate", "conf.low", "conf.high"), .funs = round.exp) %>%
    cbind(Compound = names(mod)) %>%
    left_join(cmpd.meta, by  = "Compound")
  
  # Make columns
  df$P.value    <- round(df$p.value, 3)
  df$FDR        <- round(p.adjust(df$p.value, method = "fdr"), 3)
  df$Bonferroni <- round(p.adjust(df$p.value, method = "bonferroni"), 3)
  df$CI.95      <- paste("(", df$conf.low, ", ", df$conf.high, ")", sep = "")

  # Select columns and order
  output <- df %>% 
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

library(stargazer)
stargazer(tab, summary = F, type = "html", out = "metabolite_table_selected_new.html")

# Manhattan plot
df <- bind_rows("All" = all, "Pre" = pre, "Post" = post, .id = "Group")

library(ggplot2)
ggplot(pre, aes(y = reorder(Compound, P.value), x = log10(P.value))) + 
  theme_minimal() + geom_point() + 
  #geom_vline(xintercept = -3, linetype = "dashed") +
  xlab("-log10(p-value)") +
  facet_grid(description ~ ., scales = "free_y", space = "free_y", switch= "x") +
  theme(axis.title.y = element_blank(), #axis.text.y = element_text(size=9),
        legend.position = c(0.25, 0.4),
        legend.box.background = element_rect(colour="grey")) #+
  #ggtitle("Metabolite associations with WCRF score (cal)")

# Plot data with Metafor
# Get vectors for row spacings using groups (may add compound classes later)
cmpds_ordered <- cmpd_meta %>% arrange(description) %>% mutate(row = 1:n() + (as.numeric(description)-1))
rowvec <- cmpds_ordered$row

# Investigation of Ethanol
boxplot(dat$Ethanol ~ dat$CT + dat$MENOPAUSE, varwidth = T, outline = F,
        names = c("Control, pre", "Case, pre", "Control, post", "Case, post"),
        col = "dodgerblue", ylab = "Plasma ethanol conc (scaled)")

boxplot(log(dat$ALCOHOL) ~ dat$CT + dat$MENOPAUSE, varwidth = T, outline = F,
        names = c("Control, pre", "Case, pre", "Control, post", "Case, post"),
        col = "hotpink", ylab = "Plasma ethanol conc (scaled)")

mod1 <- wilcox.test(ALCOHOL ~ CT, data = dat, subset = MENOPAUSE == 0)
mod2 <- wilcox.test(ALCOHOL ~ CT, data = dat, subset = MENOPAUSE == 1)


hist(log2(dat$Ethanol))
plot(log(dat$Ethanol), dat$ALCOHOL)
fit.e <- lm(Ethanol ~ ALCOHOL, data = dat)
summary(fit.e)

# For publication
library(ggsignif)
ggplot(dat, aes(x = as.factor(CT), y = Ethanol)) + geom_boxplot() + ylim(0, 2) +
  geom_signif(comparisons = c(0, 1), map_signif_level = T)

# Heatmap of differences---------------------

# Make a lookup df for menopausal status
match.men <- dat %>% select(MATCH, MENOPAUSE) %>% unique()

# Wide and long datasets for gplots and ggplot2
diff.wide <- dat %>% select(MATCH, CT,`3Hydroxybutyrate`:Succinate) %>% arrange(CT) %>%
  group_by(MATCH) %>% summarise_all(list(d = diff)) #%>% left_join(match.men, by = "MATCH")

diff.long <- gather(diff.wide, compound, difference, -MATCH)

# Vectors of pre and post-menopausal 
pre <- diff.wide$MENOPAUSE == 0
post <- diff.wide$MENOPAUSE == 1

hist(dat1$difference, breaks = 50)
which.max(dat1$difference)
which.min(dat1$difference)
dat1[10291, ]
dat1[34429, ]

# With heatmap2
mat1 <- as.matrix(diff.wide[pre, -1])
mat2 <- diff.wide[post, -1]

quantile.range <- quantile(as.matrix(diff.wide[, -(1:2)]), probs = seq(0, 1, 0.01))
palette.breaks <- seq(quantile.range["1%"], quantile.range["99%"], 0.1)
colpalette1 <- colorRampPalette(c("#FC8D59", "#FFFFBF", "#91CF60"))(length(palette.breaks) - 1)
colpalette2 <- colorRampPalette(c('#ef8a62','#f7f7f7','#67a9cf'))(length(palette.breaks) - 1)

library(gplots)
par(mar=c(10, 4.5, 0, 0.5))
heatmap.2(as.matrix(diff.wide[, -(1:2)]), trace = "none", col = colpalette2,
          breaks = palette.breaks, dendrogram = "col", margins = c(9, 3),
          offsetCol = 0, srtCol = 60, cexRow = 0.05)

# With ggplot
dat1 <- inner_join(diff.long, dat, by = "MATCH") %>% filter(CT == 0)

ggplot(dat1, aes(x= MATCH, y = compound, fill = difference)) + geom_tile() +
  facet_grid(. ~ MENOPAUSE, scales = "free", space = "free") +
  scale_fill_gradientn(colours = colpalette2)

# ---------------------------------------------------------------------------------------------------

# Conditional logistic regression to get odds ratios for lifestyle factors
# Adjusting for same co-variates as in original manuscript
library(survival)
# Test for associations between alcohol and BC
fit0 <- clogit(CT ~ scale(ALCOHOL) + BMI + SMK + DIABETE + RTH + DURTHSDIAG + CENTTIME + STOCKTIME + strata(MATCH), data = meta)
fit1 <- clogit(CT ~ scale(ALCOHOL) + BMI + SMK + DIABETE + RTH + DURTHSDIAG + CENTTIME + STOCKTIME + strata(MATCH), data = meta, subset = MENOPAUSE == 0) 
fit2 <- clogit(CT ~ scale(ALCOHOL) + BMI + SMK + DIABETE + RTH + DURTHSDIAG + CENTTIME + STOCKTIME + strata(MATCH), data = meta, subset = MENOPAUSE == 1) 

# fit0: 1.0784 [0.9704, 1.199]
# fit1: 0.9569 [0.77228, 1.186]
# fit2: 1.1035 [0.9739 , 1.250]

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
