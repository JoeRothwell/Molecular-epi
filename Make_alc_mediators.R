# Make mediator variables and bind to small case control (CRC1)
# Load the 3 CRC studies with all variables
load("alc_mediation.Rdata")
library(tidyverse)

# First mediator: from controls of CRC1

# Prepare metabolite matrix for PLS. Subset controls only, log2, scale
expr <- "(carn|oacid|genic|roph|ingo|Sugars)[_]"
ctlmat <- crc1 %>% filter(Cncr_Caco_Clrt == 0) %>%
  select(matches(expr), -contains("tdq")) %>% select_if(~ sum(., na.rm = T) != 0) %>% 
  log2 %>% scale

# Control subset all data
dat <- crc1[crc1$Cncr_Caco_Clrt == 0, ]

library(lme4)
adj <- function(x) residuals(lm(x ~ Center + Batch_MetBio, data = dat))
adjmat <- apply(ctlmat, 2, adj)

library(pls)
# Bind log transformed alcohol intakes to adjusted metabolite matrix
#plsdat <- bind_cols(alc = log2(dat$Qe_Alc + 0.1), adjmat) #%>% filter(!is.na(score))
plsdat <- bind_cols(alc = log(dat$Qe_Alc/12 + 0.01), adjmat)

set.seed(111)
mod <- plsr(alc ~ ., ncomp = 10, data = plsdat, validation = "CV")
# Find the number of dimensions with lowest cross validation error
cv <- RMSEP(mod)
plot(RMSEP(mod), legendpos = "topright")

# Get components with lowest RMSEP, min, "one SE" and "permutation" methods
which.min(cv$val[estimate = "adjCV", , ]) - 1
selectNcomp(mod, method = "onesigma", plot = F)
selectNcomp(mod, method = "randomization", plot = T)

# Use 2 components
mod1 <- plsr(alc ~ ., data = plsdat, ncomp = 2)
# Get coefficients
coeff <- data.frame(value = round(coef(mod)[ , 1, 1], 3))

# Get cases for prediction, log2 and scale as for discovery matrix
mat <- crc1 %>% #filter(Cncr_Caco_Clrt == 1) %>%
  select(matches(expr), -contains("tdq")) %>% select_if(~ sum(., na.rm = T) != 0) %>% 
  log2 %>% scale

crc1$M.alc <- exp(predict(mod1, mat)[,,1])
#fit3a <- clogit(update(base, ~. + M.alc), data = crc1) 
#OR for mediator 2.63 (1.43-4.81)


# Second mediator: from EPIC pooled controls

# Repeat using discovery set. Read in EPIC controls and subset metabolites
# Get cases for prediction, log2 and scale as for discovery matrix
expr <- "(carn|oacid|genic|roph|ingo|Sugars)[_]"
ctlmat <- crc1 %>% select(matches(expr), -contains("tdq")) %>% 
  select_if(~ sum(., na.rm = T) != 0) %>% log2 %>% scale

# Repeat using EPIC controls as the discovery set.
library(haven)
ctrl <- read_dta("obes_metabo.dta") %>% mutate(Study = 
  fct_recode(Study, Colonrectum_S = "Colonrectum", Colonrectum_L = "Colonrectum (bbmri)")) %>%
  separate(Batch_MetBio, into = c("batch", "rest")) %>%
  mutate(batch_no = as.numeric(flatten(str_extract_all(batch, "[0-9]+")))) %>% 
  filter(!(Study %in% c("Colonrectum_L", "Colonrectum_S")), Fasting_C == 2) %>%
  filter(Country != 6) #%>% droplevels()

# Prepare metabolite matrix for PLS. Subset controls only, log2, scale
# Note: need to remove NAs from this dataset unlike crc1
# Subset metabolites and remove ones with many missings
ctrl0 <- ctrl %>% select(matches(expr), -contains("tdq")) %>% select_if(~ sum(., na.rm = T) != 0)
ctrl0 <- ctrl0[, colSums(is.na(ctrl0)) < 697]

# Get CC matrix and discovery matrix and put columns in same order
predmat <- ctlmat[, intersect(colnames(ctrl0), colnames(ctlmat))]

library(zoo)
hm <- function(x) min(x)/2
discmat <- ctrl0[, intersect(colnames(ctrl0), colnames(ctlmat))] %>% as.matrix %>% na_if(0) %>%
  na.aggregate(FUN = hm) %>% log2 %>% scale


library(lme4)
adj   <- function(x) residuals(lm(x ~ Center + batch_no, data = ctrl))
#adj   <- function(x) residuals(lmer(x ~ Center + batch_no + Sex + (1|Study), data = ctrl))
adjmat <- apply(discmat, 2, adj)

library(pls)
# Bind log transformed alcohol intakes to adjusted metabolite matrix
#plsdat <- bind_cols(alc = log2(dat$Qe_Alc + 0.1), adjmat) #%>% filter(!is.na(score))
plsdat <- bind_cols(alc = log(ctrl$Qe_Alc/12 + 0.01), adjmat)

set.seed(111)
mod <- plsr(alc ~ ., ncomp = 10, data = plsdat, validation = "CV")
# Find the number of dimensions with lowest cross validation error
cv <- RMSEP(mod)
plot(RMSEP(mod), legendpos = "topright")

# Get components with lowest RMSEP, min, "one SE" and "permutation" methods
which.min(cv$val[estimate = "adjCV", , ]) - 1
selectNcomp(mod, method = "onesigma", plot = F)
selectNcomp(mod, method = "randomization", plot = T)

# Use 4 components
mod1 <- plsr(alc ~ ., data = plsdat, ncomp = 4)
# Get coefficients
coeff <- data.frame(value = round(coef(mod)[ , 1, 1], 3))

# Get cases for prediction, log2 and scale as for discovery matrix
crc1$M.alc1 <- exp(predict(mod1, data.frame(ctlmat))[,,1])
# OR for mediator