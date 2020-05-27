library(tidyverse)
load("predicted_score_tables2.Rdata")

# Maintain body weight within the normal range
# Be moderately physically active
# Avoid sugary drinks: Wcrf_Drinks_Cal
# Limit consumption of energy-dense foods and avoid sugary drinks
# Eat at least 5 portions of non-starchy vegetables/fruits every day: Wcrf_Fv_Cal
# Eat unprocessed cereals (grains) and/or pulses (legumes): Wcrf_Fibt_Cal
# Eat mostly foods of plant origin
# Energy dense foods: Wcrf_Ed
# Animal foods: Wcrf_Meat_Cal
# Avoid alcohol
# Breastfeed infants exclusively up to 6 months: Wcrf_Bf
# Overall score

varlist <- c("Wcrf_Bmi", "Wcrf_Pa", "Wcrf_Fwg_Cal", "Wcrf_Pf_Cal", "Wcrf_Meat_Cal", 
             "Wcrf_Alc", "Wcrf_C_Cal")
nlist <- c("Maintain normal\nbody weight", "Be moderately\nphysically active", 
           "Limit foods that\npromote weight gain", "Eat mostly\nplant foods", "Limit red and\nprocessed meat", "Avoid\nalcohol", 
           "Overall WCRF/AICR\nscore")

crc1f <- crc3.ph %>% ungroup %>% filter(Cncr_Caco_Clrt == 0) %>% select(one_of(varlist)) %>% as.matrix
crc1b <- crc1.ph %>% ungroup %>% filter(Cncr_Caco_Clrt == 0) %>% select(one_of(varlist)) %>% as.matrix
crc2b <- crc2.ph %>% ungroup %>% filter(Cncr_Caco_Clrt == 0) %>% select(one_of(varlist)) %>% as.matrix

comp2.ct1 <- crc3.ph$comp2[crc3.ph$Cncr_Caco_Clrt == 0]
comp1.ct1 <- crc1.ph$comp1[crc1.ph$Cncr_Caco_Clrt == 0]
comp1.ct2 <- crc2.ph$comp1[crc2.ph$Cncr_Caco_Clrt == 0]

library(psych)
fa1 <- corr.test(comp2.ct1, crc1f, use = "pairwise.complete.obs")
em1 <- corr.test(comp1.ct1, crc1b, use = "pairwise.complete.obs")
em2 <- corr.test(comp1.ct2, crc2b, use = "pairwise.complete.obs")

cor.ci <- bind_rows(fa1$ci, em1$ci, em2$ci, .id = "Model") %>% 
  bind_cols(tibble(component = rep(varlist, 3)))

# Partial correlations
# Subset controls only from prediction tables
dat1 <- crc3.ph[crc3.ph$Cncr_Caco_Clrt == 0, ]
dat2 <- crc1.ph[crc1.ph$Cncr_Caco_Clrt == 0, ]
dat3 <- crc2.ph[crc2.ph$Cncr_Caco_Clrt == 0, ]

# Function to get partial correlations, omitting NAs
get.pcor <- function(x, dat, predscore) {
  mod1 <- lm(x[!is.na(x)] ~         Qe_Energy + L_School + Smoke_Int + Smoke_Stat + Height_C, data = dat[!is.na(x), ] )
  mod2 <- lm(predscore[!is.na(x)] ~ Qe_Energy + L_School + Smoke_Int + Smoke_Stat + Height_C, data = dat[!is.na(x), ] )
  output <- cor.test(residuals(mod1), residuals(mod2))
}

# Run function and extract data 
pcor1 <- apply(crc1f, 2, function(x, dat, predscore) get.pcor(x, dat1, comp2.ct1)) 
pcor2 <- apply(crc1b, 2, function(x, dat, predscore) get.pcor(x, dat2, comp1.ct1))
pcor3 <- apply(crc2b, 2, function(x, dat, predscore) get.pcor(x, dat3, comp1.ct2))

pcor.ci <- bind_rows("Study A, fatty acids" = map_df(pcor1, tidy), 
                     "Study A, endogenous metabolites" = map_df(pcor2, tidy), 
                     "Study B, endogenous metabolites" = map_df(pcor3, tidy), .id = "Model") %>%
  bind_cols(tibble(component = rep(varlist, 3)))

# Plot data: partial correlation
ggplot(pcor.ci, aes(x=fct_inorder(component), y=estimate, colour = fct_inorder(Model))) + 
theme_bw() +
scale_colour_manual(values = c("red","blue","blue")) +
geom_hline(yintercept = 0, linetype = "dashed") +
geom_point(size = 2, position = position_dodge(width = 0.5)) +
ylab("Partial Pearson correlation") +
geom_errorbar(width=0.2, aes(ymin = conf.low, ymax = conf.high), 
              position= position_dodge(width = 0.5)) +
scale_x_discrete(labels = nlist) +
theme(axis.title.x = element_blank(), legend.position = "bottom",
      legend.title = element_blank())



# Raw correlations
library(psych)
fa1 <- corr.test(comp2.ct1, crc1f, use = "pairwise.complete.obs")
em1 <- corr.test(comp1.ct1, crc1b, use = "pairwise.complete.obs")
em2 <- corr.test(comp1.ct2, crc2b, use = "pairwise.complete.obs")

cor.ci <- bind_rows(fa1$ci, em1$ci, em2$ci, .id = "Model") %>% 
  bind_cols(tibble(component = rep(varlist, 3)))

# Plot data: raw correlation
library(ggplot2)
cor.ci %>%
  ggplot(aes(x=fct_inorder(component), y=r, shape = signature)) + 
  theme_bw() +
  geom_hline(yintercept = 0, linetype = "dashed") +
  geom_point(size = 2, position = position_dodge(width = 0.5)) +
  ylab("Pearson correlation") + xlab("Fatty acid or endogenous signature") +
  geom_errorbar(width=0.2, aes(ymin = lower, ymax = upper), 
                position= position_dodge(width = 0.5), colour="black") +
  theme(axis.title.x = element_blank())


# Old: facetted
ggplot(cor.ci, aes(x = Model, y = r)) + #geom_point() + 
  geom_hline(yintercept = 0, linetype = "dashed") +
  geom_errorbar(width=.2, aes(ymin=lower, ymax=upper), colour="black") +
  geom_point() + theme_minimal() +
  ylab("Pearson correlation") + xlab("Fatty acid or endogenous signature") +
  facet_grid(. ~ fct_inorder(component), scales = "free_x") +
  theme(panel.spacing=unit(1.5,"lines"))
