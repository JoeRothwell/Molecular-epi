library(survival)
library(broom)
library(ggrepel)
source("BC_prep_data.R")

fits0 <- apply(ints, 2, function(x) clogit(CT ~ BMI + SMK + DIABETE + RTH + ALCOHOL + 
  DURTHSBMB + CENTTIME + STOCKTIME + strata(MATCH) + x, data = meta))

t1 <- map_df(fits0, tidy) %>% filter(str_detect(term, "x")) %>% bind_cols(cmpd.meta) %>% 
  arrange(description)

p1 <- ggplot(t1, aes(exp(estimate), log10(p.value))) + geom_point() + theme_bw() +
  xlim(0.8, 1.2) +
  scale_y_reverse(breaks = c(-2, -1, 0), labels = function(x) 10^x) +
  xlab("Odds ratio per SD increase concentration") + ylab("P-value") +
  geom_text_repel(aes(label = display_name), size = 3,
                  data = t1[t1$p.value < 0.3, ]) +
  geom_hline(yintercept = log10(0.05), linetype = "dashed") +
  geom_vline(xintercept = 1, linetype = "dashed") +
  ggtitle("All subjects")

# Pre-menopausal
meta1 <- meta[pre, ]
# Note: need to remove hormone treatment therapy variable DURTHSBMB
fits1 <- apply(ints[pre, ], 2, function(x) clogit(CT ~ BMI + SMK + DIABETE + RTH + ALCOHOL + 
  CENTTIME + STOCKTIME + strata(MATCH) + x, data = meta[pre, ]))

t1 <- map_df(fits1, tidy) %>% filter(str_detect(term, "x")) %>% bind_cols(cmpd.meta) %>% 
  arrange(description)

p2 <- ggplot(t1, aes(exp(estimate), log10(p.value))) + geom_point() + theme_bw() +
  xlim(0.3, 1.7) + 
  scale_y_reverse(breaks = c(-3, -2, -1, 0), labels = function(x) 10^x) +
  xlab("Odds ratio per SD increase concentration") + ylab("P-value") +
  geom_text_repel(aes(label = display_name), size = 3,
                  data = t1[t1$p.value < 0.04, ] ) +
  geom_hline(yintercept = c(log10(0.05), log10(0.014)), linetype = c("dashed", "dotted")) +
  geom_vline(xintercept = 1, linetype = "dashed") +
  ggtitle("Pre-menopausal")

# Post menopausal
meta2 <- meta[post, ]
fits2 <- apply(ints[post, ], 2, function(x) clogit(CT ~ BMI + SMK + DIABETE + RTH + ALCOHOL + 
  DURTHSBMB + CENTTIME + STOCKTIME + strata(MATCH) + x, data = meta2))

t1 <- map_df(fits2, tidy) %>% filter(str_detect(term, "x")) %>% bind_cols(cmpd.meta) %>% 
  arrange(description)

p3 <- ggplot(t1, aes(exp(estimate), log10(p.value))) + geom_point() + theme_bw() +
  xlim(0.85, 1.15) +  
  scale_y_reverse(breaks = c(-2, -1, 0), labels = function(x) 10^x) +
  xlab("Odds ratio per SD increase concentration") + ylab("P-value") +
  geom_text_repel(aes(label = display_name), size = 3,
                  data = t1[t1$p.value < 0.15, ]) +
  geom_hline(yintercept = log10(0.05), linetype = "dashed") +
  geom_vline(xintercept = 1, linetype = "dashed") +
  ggtitle("Post-menopausal")


# Sensitivity analysis adjusted ethanol
# All subjects
eth <- ints[pre, 14]
meta.eth <- cbind(meta1, ETHANOL = eth)

fits1e <- apply(ints[pre, -14], 2, function(x) clogit(CT ~ BMI + SMK + DIABETE + RTH + ALCOHOL +
          ETHANOL +
          CENTTIME + STOCKTIME + strata(MATCH) + x, data = meta.eth))

t1 <- map_df(fits1e, tidy) %>% filter(str_detect(term, "x")) %>% bind_cols(cmpd.meta[-14, ]) %>% arrange(description)

p4 <- ggplot(t1, aes(exp(estimate), log10(p.value))) + geom_point() + theme_bw() +
  xlim(c(0.4, 1.6)) + 
  scale_y_reverse(breaks = c(-3, -2, -1, 0), labels = function(x) 10^x) +
  xlab("Odds ratio per SD increase concentration") + ylab("P-value") +
  geom_text_repel(aes(label = display_name), size = 3,
                  data = t1[t1$p.value < 0.04, ] ) +
  geom_hline(yintercept = c(log10(0.05), log10(0.001)), linetype = c("dashed", "dotted")) +
  geom_vline(xintercept = 1, linetype = "dashed") +
  ggtitle("Pre-menopausal, adjusted for plasma ethanol")


library(cowplot)
plot_grid(p1, p3, p2, p4, labels = c('A', 'B', 'C', 'D'), label_size = 12)  



# Pre-menopausal
meta1e <- meta.e[pre, ]
# Note: need to remove hormone treatment therapy variable DURTHSBMB
fits1 <- apply(ints[pre, ], 2, function(x) clogit(CT ~ BMI + SMK + DIABETE + RTH + ALCOHOL +
                                                    ethanol +
  CENTTIME + STOCKTIME + strata(MATCH) + x, data = meta1e[pre, ]))

t1 <- map_df(fits1, tidy) %>% filter(str_detect(term, "x")) %>% bind_cols(cmpd.meta) %>% 
  arrange(description)

p2 <- ggplot(t1, aes(exp(estimate), log10(p.value))) + geom_point() + theme_bw() +
  xlim(0.3, 1.7) + 
  scale_y_reverse(breaks = c(-3, -2, -1, 0), labels = function(x) 10^x) +
  xlab("Odds ratio per SD increase concentration") + ylab("P-value") +
  geom_text_repel(aes(label = display_name), size = 3,
                  data = t1[t1$p.value < 0.04, ] ) +
  geom_hline(yintercept = c(log10(0.05), log10(0.014)), linetype = "dashed") +
  geom_vline(xintercept = 1, linetype = "dashed") +
  ggtitle("Pre-menopausal")
