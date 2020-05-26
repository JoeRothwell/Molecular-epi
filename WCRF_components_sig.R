source("CRC_prep_data")
source("CRC_get_signatures")

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

crc1.wcrf <- crc1.ph %>% ungroup %>% filter(Cncr_Caco_Clrt == 0) %>% select(one_of(varlist)) %>% as.matrix
crc2.wcrf <- crc2.ph %>% ungroup %>% filter(Cncr_Caco_Clrt == 0) %>% select(one_of(varlist)) %>% as.matrix
crc3.wcrf <- crc3.ph %>% ungroup %>% filter(Cncr_Caco_Clrt == 0) %>% select(one_of(varlist)) %>% as.matrix

comp1.ct1 <- crc1.ph$comp1[crc1.ph$Cncr_Caco_Clrt == 0]
comp1.ct2 <- crc2.ph$comp1[crc2.ph$Cncr_Caco_Clrt == 0]
comp2.ct3 <- crc3.ph$comp2[crc3.ph$Cncr_Caco_Clrt == 0]

em1 <- cor(comp1.ct1, crc1.wcrf, use = "pairwise.complete.obs")
em2 <- cor(comp1.ct2, crc2.wcrf, use = "pairwise.complete.obs")
fa1 <- cor(comp2.ct3, crc3.wcrf, use = "pairwise.complete.obs")

cor.wcrf <- t(rbind(em1, em2, fa1))
colnames(cor.wcrf) <- c("Bioc1", "Bioc2", "FAs")
heatmap(cor.wcrf, Colv = NA, Rowv = NA)

library(psych)
em1 <- corr.test(comp1.ct1, crc1.wcrf, use = "pairwise.complete.obs")
em2 <- corr.test(comp1.ct2, crc2.wcrf, use = "pairwise.complete.obs")
fa1 <- corr.test(comp2.ct3, crc3.wcrf, use = "pairwise.complete.obs")

cor.ci <- bind_rows(em1$ci, em2$ci, fa1$ci, .id = "Signature") %>% 
  bind_cols(tibble(component = rep(varlist, 3)))
# Facetted plot?
ggplot(cor.ci, aes(x = Signature, y = r)) + #geom_point() + 
  geom_hline(yintercept = 0, linetype = "dashed") +
  geom_errorbar(width=.2, aes(ymin=lower, ymax=upper), colour="black") +
  geom_point() + theme_minimal() +
  ylab("Pearson correlation") +
  facet_grid(. ~ fct_inorder(component), scales = "free_x") +
  theme(panel.spacing=unit(2,"lines"))
