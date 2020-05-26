source("CRC_prep_data")
source("CRC_get_signatures")

# Maintain body weight within the normal range
# Be moderately physically active
# Avoid sugary drinks
# Limit consumption of energy-dense foods and avoid sugary drinks
# Eat at least 5 portions of non-starchy vegetables/fruits every day
# Eat unprocessed cereals (grains) and/or pulses (legumes)
# Eat mostly foods of plant origin
# Energy dense foods
# Animal foods
# Avoid alcohol
# Breastfeed infants exclusively up to 6 months
# Overall score

varlist <- c("Wcrf_Bmi", "Wcrf_Pa", "Wcrf_Drinks_Cal", "Wcrf_Fwg_Cal", "Wcrf_Fv_Cal", "Wcrf_Fibt_Cal",
             "Wcrf_Pf_Cal", "Wcrf_Ed", "Wcrf_Meat_Cal", "Wcrf_Alc", "Wcrf_Bf", "Wcrf_C_Cal")

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

# Facetted plot?
