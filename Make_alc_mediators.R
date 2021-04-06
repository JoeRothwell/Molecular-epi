# Make mediator variables and bind to small case control (CRC1)
# Load the 3 CRC studies with all variables
load("alc_mediation.Rdata")

# Get small case-control study
#meta <- read_dta("clrt_caco.dta") 
#crc1 <- read_sas("clrt_caco_metabo.sas7bdat") %>% filter(!is.na(Batch_MetBio) & Country != 6)
#crc1 <- inner_join(crc1, meta, by = "Idepic", suffix = c("_1", "")) 

library(tidyverse)

# Mediators M1, M2, M3, derived from controls of the CRC case-control

# Subset case-control metabolite matrix, log, scale
expr1 <- "(carn|oacid|genic|roph|ingo|Sugars)[_]"
expr2 <- "oacid[_]"
expr3 <- "carn|roph|ingo[_]"

# Get matrices for all biocrates compounds, amino acids and lipids
mat1 <- crc1 %>% select(matches(expr1), -contains("tdq")) %>% log2 %>% scale
mat2 <- crc1 %>% select(matches(expr2), -contains("tdq")) %>% log2 %>% scale
mat3 <- crc1 %>% select(matches(expr3), -contains("tdq")) %>% log2 %>% scale

# Subset controls
vec <- crc1$Cncr_Caco_Clrt == 0 & crc1$Qe_Alc != 0
ctlmat1 <- mat1[vec, ]
ctlmat2 <- mat2[vec, ]
ctlmat3 <- mat3[vec, ]
dat <- crc1[vec, ]

# Adjust matrix and bind log transformed alcohol intakes
library(lme4)
adj <- function(x) residuals(lm(x ~ Center + Batch_MetBio, data = dat))
adjmat1 <- apply(ctlmat1, 2, adj)
adjmat2 <- apply(ctlmat2, 2, adj)
adjmat3 <- apply(ctlmat3, 2, adj)

plsdat1 <- bind_cols(alc = log(dat$Qe_Alc/12), adjmat1)
plsdat2 <- bind_cols(alc = log(dat$Qe_Alc/12), adjmat2)
plsdat3 <- bind_cols(alc = log(dat$Qe_Alc/12), adjmat3)

# Train PLS. Find the number of dimensions with lowest cross validation error
library(pls)
set.seed(111)
mod <- plsr(alc ~ ., ncomp = 10, data = plsdat1, validation = "CV")
cv <- RMSEP(mod)
plot(cv, legendpos = "topright")

# Get no. components by lowest RMSEP, min, "one SE" and "permutation" methods
which.min(cv$val[estimate = "adjCV", , ]) - 1
selectNcomp(mod, method = "onesigma", plot = F)
selectNcomp(mod, method = "randomization", plot = T)

# Run PLS with 2 components and get coefficients, predict for mat
mod1 <- plsr(alc ~ ., data = plsdat1, ncomp = 1)
coeff <- data.frame(value = round(coef(mod1)[ , 1, 1], 3))
crc1$M1 <- exp(predict(mod1, mat1)[,,1])


# Amino acids only
mod <- plsr(alc ~ ., ncomp = 10, data = plsdat2, validation = "CV")
cv <- RMSEP(mod)
plot(cv, legendpos = "topright")

# Get no. components by lowest RMSEP, min, "one SE" and "permutation" methods
which.min(cv$val[estimate = "adjCV", , ]) - 1
selectNcomp(mod, method = "onesigma", plot = F)
selectNcomp(mod, method = "randomization", plot = T)

# Run PLS with 2 components and get coefficients, predict for mat
mod1 <- plsr(alc ~ ., data = plsdat2, ncomp = 1)
coeff <- data.frame(value = round(coef(mod1)[ , 1, 1], 3))
crc1$M2 <- exp(predict(mod1, mat2)[,,1])


# Lipid metabolites only
mod <- plsr(alc ~ ., ncomp = 10, data = plsdat3, validation = "CV")
cv <- RMSEP(mod)
plot(cv, legendpos = "topright")

# Get no. components by lowest RMSEP, min, "one SE" and "permutation" methods
which.min(cv$val[estimate = "adjCV", , ]) - 1
selectNcomp(mod, method = "onesigma", plot = F)
selectNcomp(mod, method = "randomization", plot = T)

# Run PLS with 2 components and get coefficients, predict for mat
mod1 <- plsr(alc ~ ., data = plsdat3, ncomp = 1)
coeff <- data.frame(value = round(coef(mod1)[ , 1, 1], 3))
crc1$M3 <- exp(predict(mod1, mat3)[,,1])

#fit3a <- clogit(update(base, ~. + M.alc), data = crc1) 
#OR for mediator 2.63 (1.43-4.81)

# Mediator M4 derived from EPIC pooled controls

# Get case-control for prediction, log2 and scale as for discovery matrix
expr <- "(carn|oacid|genic|roph|ingo|Sugars)[_]"
mat <- crc1 %>% select(matches(expr), -contains("tdq")) %>% log2 %>% scale

library(haven)
ctrl <- read_dta("obes_metabo.dta") %>% mutate(Study = 
  fct_recode(Study, Colonrectum_S = "Colonrectum", Colonrectum_L = "Colonrectum (bbmri)")) %>%
  separate(Batch_MetBio, into = c("batch", "rest")) %>%
  mutate(batch_no = as.numeric(flatten(str_extract_all(batch, "[0-9]+")))) %>% 
  filter(!(Study %in% c("Colonrectum_L", "Colonrectum_S")), Fasting_C == 2) %>%
  filter(Country != 6) #%>% droplevels()

# Prepare metabolite matrix for PLS. Subset controls only, log2, scale
# Note: need to remove NAs from this dataset because not imputed like crc1
# Subset metabolites and remove ones with many missings
ctrl0 <- ctrl %>% select(matches(expr), -contains("tdq")) %>% select_if(~ sum(., na.rm = T) != 0)
ctrl0 <- ctrl0[, colSums(is.na(ctrl0)) < 697]

# Get CC matrix and discovery matrix and put columns in same order
#predmat <- mat[, intersect(colnames(ctrl0), colnames(mat))]

library(zoo)
hm <- function(x) min(x)/2
discmat <- ctrl0[, intersect(colnames(ctrl0), colnames(mat))] %>% as.matrix %>% na_if(0) %>%
  na.aggregate(FUN = hm) %>% log2 %>% scale

# Exclude zero alcohol consumers N = 1741->1635
discmat <- discmat[ctrl$Qe_Alc != 0, ]
dat1 <- ctrl[ctrl$Qe_Alc != 0, ]

# Adjust matrix and bind log transformed alcohol intakes
library(lme4)
adj   <- function(x) residuals(lm(x ~ Center + batch_no, data = dat1))
adjmat <- apply(discmat, 2, adj)
plsdat <- bind_cols(alc = log(dat1$Qe_Alc/12), adjmat)

# Train PLS. Find the number of dimensions with lowest cross validation error
library(pls)
set.seed(111)
mod <- plsr(alc ~ ., ncomp = 10, data = plsdat, validation = "CV")
cv <- RMSEP(mod)
plot(cv, legendpos = "topright")

# Get components with lowest RMSEP, min, "one SE" and "permutation" methods
which.min(cv$val[estimate = "adjCV", , ]) - 1
selectNcomp(mod, method = "onesigma", plot = F)
selectNcomp(mod, method = "randomization", plot = T)

# Fit PLSR with 4 components and get coefficients
mod1 <- plsr(alc ~ ., data = plsdat, ncomp = 3)
coeff <- data.frame(value = round(coef(mod)[ , 1, 1], 3))

# Get cases for prediction, log2 and scale as for discovery matrix
crc1$M4 <- exp(predict(mod1, data.frame(mat))[,,1])

# Fatty acids mediator from controls of CRC1
crc3 <- crc3.ph %>% ungroup()

# Remove problem FAs and ones with > 20% CV
concs <- crc3 %>% select(P14_0 : PCLA_9t_11c, -P24_0, -P20_0)
concs <- concs %>% select(-P14_1n_5, -P18_1n_7t, -P22_0, -P18_2n_6tt)

# Subset controls
vec <- crc3$Cncr_Caco_Clrt == 0 & crc3$Qe_Alc != 0
ctlmat <- concs[vec, ]
dat <- crc3[vec, ]

# Adjust matrix and bind log transformed alcohol intakes
library(lme4)
adj <- function(x) residuals(lm(x ~ Center + Batch_Fa, data = dat))
adjmat <- apply(ctlmat, 2, adj)
plsdat <- bind_cols(alc = log(dat$Qe_Alc/12), adjmat)

# Train PLS. Find the number of dimensions with lowest cross validation error
library(pls)
set.seed(111)
mod <- plsr(alc ~ ., ncomp = 10, data = plsdat, validation = "CV")
cv <- RMSEP(mod)
plot(cv, legendpos = "topright")

# Get components with lowest RMSEP, min, "one SE" and "permutation" methods
which.min(cv$val[estimate = "adjCV", , ]) - 1
selectNcomp(mod, method = "onesigma", plot = F)
selectNcomp(mod, method = "randomization", plot = T)

# Fit PLSR with 4 components and get coefficients
mod1 <- plsr(alc ~ ., data = plsdat, ncomp = 2)
coeff <- data.frame(value = round(coef(mod)[ , 1, 1], 3))

# Get cases for prediction, log2 and scale as for discovery matrix
crc3$M5 <- exp(predict(mod1, data.frame(concs))[,,1])
