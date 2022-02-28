# EPIC P180 and P150 studies now pooled instead of meta-analysed. 
# p180 (small) had 740 subjects, p150 (large) had 556. Note: no NAs and no zeros in these data
source("CRC_prep_data.R")
rm(list = ls(pattern = "colon|dist|prox|rect|crc1|crc2"))

library(survival)
# Model co-variates with and without matching factors
#base <- Cncr_Caco_Clrt ~ Bmi_Cat + Smoke_Stat + Alc_Drinker + Pa_Index
base <- Cncr_Caco_Clrt ~ Bmi_Cat + Smoke_Stat + Alc_Drinker + Pa_Index + Age_Blood +
  Center

# Filter amino acids with over 31% missings
# Update: keep all amino acids and give missings instead (to keep p180 AAs)
mat0 <- crc %>% select(contains("Aminoacid_")) #%>% 
#select_if(~ sum(is.na(.)) < (nrow(crc) * 0.31))
mat <- crc %>% select(colnames(mat0)) #1308 w/o Greece
scalemat <- scale(mat)

ala <- scalemat[, 1]
arg <- scalemat[, 2]
asn <- scalemat[, 3]
asp <- scalemat[, 4]
cit <- scalemat[, 5]
gln <- scalemat[, 6]
glu <- scalemat[, 7]
gly <- scalemat[, 8]
his <- scalemat[, 9]
ile <- scalemat[, 10]
leu <- scalemat[, 11]
lys <- scalemat[, 12]
met <- scalemat[, 13]
orn <- scalemat[, 14]
phe <- scalemat[, 15]
pro <- scalemat[, 16]
ser <- scalemat[, 17]
thr <- scalemat[, 18]
trp <- scalemat[, 19]
tyr <- scalemat[, 20]
val <- scalemat[, 21]

# Heterogeneity by sex
library(lmtest)
mod1 <- glm(update(base, ~. + ala + Sex), data = crc, family = "binomial")
mod2 <- glm(update(base, ~. + ala * Sex), data = crc, family = "binomial")
lrtest(mod1, mod2)
# pHET = 0.4376

mod1 <- glm(update(base, ~. + arg + Sex), data = crc, family = "binomial")
mod2 <- glm(update(base, ~. + arg * Sex), data = crc, family = "binomial")
lrtest(mod1, mod2)
# pHET = 0.2415

mod1 <- glm(update(base, ~. + asn + Sex), data = crc, family = "binomial")
mod2 <- glm(update(base, ~. + asn * Sex), data = crc, family = "binomial")
lrtest(mod1, mod2)
# pHET = 0.3003

mod1 <- glm(update(base, ~. + asp + Sex), data = crc, family = "binomial")
mod2 <- glm(update(base, ~. + asp * Sex), data = crc, family = "binomial")
lrtest(mod1, mod2)
# pHET = 0.5491

mod1 <- glm(update(base, ~. + cit + Sex), data = crc, family = "binomial")
mod2 <- glm(update(base, ~. + cit * Sex), data = crc, family = "binomial")
lrtest(mod1, mod2)
# pHET = 0.07931

mod1 <- glm(update(base, ~. + gln + Sex), data = crc, family = "binomial")
mod2 <- glm(update(base, ~. + gln * Sex), data = crc, family = "binomial")
lrtest(mod1, mod2)
# pHET = 0.4376

mod1 <- glm(update(base, ~. + glu + Sex), data = crc, family = "binomial")
mod2 <- glm(update(base, ~. + glu * Sex), data = crc, family = "binomial")
lrtest(mod1, mod2)
# pHET = 0.3944

mod1 <- glm(update(base, ~. + gly + Sex), data = crc, family = "binomial")
mod2 <- glm(update(base, ~. + gly * Sex), data = crc, family = "binomial")
lrtest(mod1, mod2)
# pHET = 0.3595

mod1 <- glm(update(base, ~. + his + Sex), data = crc, family = "binomial")
mod2 <- glm(update(base, ~. + his * Sex), data = crc, family = "binomial")
lrtest(mod1, mod2)
# pHET = 0.5426

mod1 <- glm(update(base, ~. + ile + Sex), data = crc, family = "binomial")
mod2 <- glm(update(base, ~. + ile * Sex), data = crc, family = "binomial")
lrtest(mod1, mod2)
# pHET = 0.2351

mod1 <- glm(update(base, ~. + leu + Sex), data = crc, family = "binomial")
mod2 <- glm(update(base, ~. + leu * Sex), data = crc, family = "binomial")
lrtest(mod1, mod2)
# pHET = 0.379

mod1 <- glm(update(base, ~. + lys + Sex), data = crc, family = "binomial")
mod2 <- glm(update(base, ~. + lys * Sex), data = crc, family = "binomial")
lrtest(mod1, mod2)
# pHET = 0.7577

mod1 <- glm(update(base, ~. + met + Sex), data = crc, family = "binomial")
mod2 <- glm(update(base, ~. + met * Sex), data = crc, family = "binomial")
lrtest(mod1, mod2)
# pHET = 0.6444

mod1 <- glm(update(base, ~. + orn + Sex), data = crc, family = "binomial")
mod2 <- glm(update(base, ~. + orn * Sex), data = crc, family = "binomial")
lrtest(mod1, mod2)
# pHET = 0.9053

mod1 <- glm(update(base, ~. + phe + Sex), data = crc, family = "binomial")
mod2 <- glm(update(base, ~. + phe * Sex), data = crc, family = "binomial")
lrtest(mod1, mod2)
# pHET = 0.1744

mod1 <- glm(update(base, ~. + pro + Sex), data = crc, family = "binomial")
mod2 <- glm(update(base, ~. + pro * Sex), data = crc, family = "binomial")
lrtest(mod1, mod2)
# pHET = 0.3485

mod1 <- glm(update(base, ~. + ser + Sex), data = crc, family = "binomial")
mod2 <- glm(update(base, ~. + ser * Sex), data = crc, family = "binomial")
lrtest(mod1, mod2)
# pHET = 0.4515

mod1 <- glm(update(base, ~. + thr + Sex), data = crc, family = "binomial")
mod2 <- glm(update(base, ~. + thr * Sex), data = crc, family = "binomial")
lrtest(mod1, mod2)
# pHET = 0.1919

mod1 <- glm(update(base, ~. + trp + Sex), data = crc, family = "binomial")
mod2 <- glm(update(base, ~. + trp * Sex), data = crc, family = "binomial")
lrtest(mod1, mod2)
# pHET = 0.6753

mod1 <- glm(update(base, ~. + tyr + Sex), data = crc, family = "binomial")
mod2 <- glm(update(base, ~. + tyr * Sex), data = crc, family = "binomial")
lrtest(mod1, mod2)
# pHET = 0.9737

mod1 <- glm(update(base, ~. + val + Sex), data = crc, family = "binomial")
mod2 <- glm(update(base, ~. + val * Sex), data = crc, family = "binomial")
lrtest(mod1, mod2)
# pHET = 0.8024

# Heterogeneity by subsite in UK Biobank, colon and rectal
library(readxl)
library(metafor)
results <- read_xlsx("ukb_aminoacids_het.xlsx")
gln <- rma(yi = estimate, sei = std.error, dat = results, #method="FE", 
              subset = c(11,38))
hist <- rma(yi = estimate, sei = std.error, dat = results, #method="FE", 
            subset = c(13,40))
malist <- list(gln, hist)
output <- lapply(malist, "[", c("QEp", "I2"))
phets <- do.call(rbind, output)
