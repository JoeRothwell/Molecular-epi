library(survival)
library(broom)
library(ggrepel)
source("BC_prep_data.R")

fits0 <- apply(ints, 2, function(x) clogit(CT ~ BMI + SMK + DIABETE + RTH + ALCOHOL + 
  DURTHSBMB + CENTTIME + STOCKTIME + strata(MATCH) + x, data = meta))

t1 <- map_df(fits0, tidy) %>% filter(str_detect(term, "x")) %>% 
  mutate(p.valueFDR = p.adjust(p.value, method = "fdr")) %>% bind_cols(cmpd.meta)

# Create base plot to cut down code
base <- ggplot(t1, aes(exp(estimate), log10(p.value))) + geom_point() + theme_bw() +
  xlab("Odds ratio per SD increase concentration") + ylab("P-value") +
  geom_vline(xintercept = 1, linetype = "dashed")

p1 <- 
  base %+% xlim(0.8, 1.2) +
  scale_y_reverse(breaks = c(-2, -1, 0), labels = function(x) 10^x) +
  geom_text_repel(aes(label = display_name), size = 3, data = t1[t1$p.value < 0.3, ]) + 
  geom_hline(yintercept = log10(0.05), linetype = "dashed") + ggtitle("All subjects")

# Pre-menopausal
meta1 <- meta[pre, ]
# Note: need to remove hormone treatment therapy variable DURTHSBMB
fits1 <- apply(ints[pre, ], 2, function(x) clogit(CT ~ BMI + SMK + DIABETE + RTH + ALCOHOL + 
  CENTTIME + STOCKTIME + strata(MATCH) + x, data = meta[pre, ]))

t2 <- map_df(fits1, tidy) %>% filter(str_detect(term, "x")) %>% 
  mutate(p.valueFDR = p.adjust(p.value, method = "fdr")) %>% bind_cols(cmpd.meta)

p2 <- 
  base %+% t2 + xlim(0.3, 1.7) +
  scale_y_reverse(breaks = c(-3, -2, -1, 0), labels = function(x) 10^x) +
  geom_text_repel(aes(label = display_name), size = 3, data = t2[t2$p.value < 0.04, ] ) +
  geom_hline(yintercept = c(log10(0.05), log10(0.014)), linetype = c("dashed", "dotted")) +
  ggtitle("Pre-menopausal")


# Post menopausal
meta2 <- meta[post, ]
fits2 <- apply(ints[post, ], 2, function(x) clogit(CT ~ BMI + SMK + DIABETE + RTH + ALCOHOL + 
  DURTHSBMB + CENTTIME + STOCKTIME + strata(MATCH) + x, data = meta2))

t3 <- map_df(fits2, tidy) %>% filter(str_detect(term, "x")) %>% 
  mutate(p.valueFDR = p.adjust(p.value, method = "fdr")) %>% bind_cols(cmpd.meta)

p3 <-
  base %+% t3 + xlim(0.85, 1.15) + 
  scale_y_reverse(breaks = c(-2, -1, 0), labels = function(x) 10^x) +
  geom_text_repel(aes(label = display_name), size = 3, data = t3[t3$p.value < 0.15, ]) +
  geom_hline(yintercept = log10(0.05), linetype = "dashed") +
  ggtitle("Post-menopausal")


# Sensitivity analysis adjusted ethanol
# All subjects
eth <- ints[pre, 14]
meta.eth <- cbind(meta1, ETHANOL = eth)

fits1e <- apply(ints[pre, -14], 2, function(x) clogit(CT ~ BMI + SMK + DIABETE + RTH + ALCOHOL +
          ETHANOL + CENTTIME + STOCKTIME + strata(MATCH) + x, data = meta.eth))

t4 <- map_df(fits1e, tidy) %>% filter(str_detect(term, "x")) %>% 
  mutate(p.valueFDR = p.adjust(p.value, method = "fdr")) %>% bind_cols(cmpd.meta[-14, ])


p4 <- 
  base %+% t4 + xlim(c(0.4, 1.6)) +
  scale_y_reverse(breaks = c(-3, -2, -1, 0), labels = function(x) 10^x) +
  geom_text_repel(aes(label = display_name), size = 3, data = t4[t4$p.value < 0.04, ] ) +
  geom_hline(yintercept = c(log10(0.05), log10(0.001)), linetype = c("dashed", "dotted")) +
  ggtitle("Pre-menopausal, adjusted for plasma ethanol")
  
library(cowplot)
plot_grid(p1, p3, p2, p4, labels = c('A', 'B', 'C', 'D'), label_size = 12)  


# Pre-menopausal, exclude first 2 years of follow up
meta3 <- meta[pre0, ]
# Note: need to remove hormone treatment therapy variable DURTHSBMB
fits1 <- apply(ints[pre0, ], 2, function(x) clogit(CT ~ BMI + SMK + DIABETE + RTH + ALCOHOL + 
    CENTTIME + STOCKTIME + strata(MATCH) + x, data = meta3))

t5 <- map_df(fits1, tidy) %>% filter(str_detect(term, "x")) %>% 
  mutate(p.valueFDR = p.adjust(p.value, method = "fdr")) %>% bind_cols(cmpd.meta)

p5 <- 
  base %+% t5 + xlim(0.3, 1.7) + 
  scale_y_reverse(breaks = c(-3, -2, -1, 0), labels = function(x) 10^x) +
  geom_text_repel(aes(label = display_name), size = 3, data = t5[t5$p.value < 0.04, ] ) +
  geom_hline(yintercept = c(log10(0.05), log10(0.014)), linetype = c("dashed", "dotted")) +
  ggtitle("Pre-menopausal, follow-up < 2 y excluded")

# Stratification by age
table(Age = meta$Age1, Menopause = meta$MENOPAUSE)

# Under 55
meta4 <- meta[agelo, ]
fits0 <- apply(ints[agelo, ], 2, function(x) clogit(CT ~ BMI + SMK + DIABETE + RTH + ALCOHOL + 
  DURTHSBMB + CENTTIME + STOCKTIME + strata(MATCH) + x, data = meta4))

t6 <- map_df(fits0, tidy) %>% filter(str_detect(term, "x")) %>% 
  mutate(p.valueFDR = p.adjust(p.value, method = "fdr")) %>% bind_cols(cmpd.meta)

p6 <- 
  base %+% t6 + #xlim(0.8, 1.2) +
  scale_y_reverse(breaks = c(-2, -1, 0), labels = function(x) 10^x) +
  geom_text_repel(aes(label = display_name), size = 3, data = t6[t6$p.value < 0.3, ]) +
  geom_hline(yintercept = log10(0.05), linetype = "dashed") + ggtitle("Age < 55")

# Over 55
meta5 <- meta[agehi, ]
fits0 <- apply(ints[agehi, ], 2, function(x) clogit(CT ~ BMI + SMK + DIABETE + RTH + ALCOHOL + 
   DURTHSBMB + CENTTIME + STOCKTIME + strata(MATCH) + x, data = meta5))

t7 <- map_df(fits0, tidy) %>% filter(str_detect(term, "x")) %>% 
  mutate(p.valueFDR = p.adjust(p.value, method = "fdr")) %>% bind_cols(cmpd.meta)

p7 <- 
  base %+% t7 + #xlim(0.8, 1.2) +
  scale_y_reverse(breaks = c(-2, -1, 0), labels = function(x) 10^x) +
  geom_text_repel(aes(label = display_name), size = 3, data = t7[t7$p.value < 0.3, ]) +
  geom_hline(yintercept = log10(0.05), linetype = "dashed") + ggtitle("Age > 55")













