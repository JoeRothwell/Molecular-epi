# Make mediator variables and bind to small case control (CRC1)
# Load the 3 CRC studies with all variables
load("alc_mediation.Rdata")
library(tidyverse)

# First mediator derived from controls of CRC1

# Subset case-control metabolite matrix, log, scale
expr <- "(carn|oacid|genic|roph|ingo|Sugars)[_]"
mat <- crc1 %>% select(matches(expr), -contains("tdq")) %>% log2 %>% scale

# Subset controls 
ctlmat <- mat[crc1$Cncr_Caco_Clrt == 0 & crc1$Qe_Alc != 0, ]
dat <- crc1[crc1$Cncr_Caco_Clrt == 0 & crc1$Qe_Alc != 0, ]

# Adjust matrix and bind log transformed alcohol intakes
library(lme4)
adj <- function(x) residuals(lm(x ~ Center + Batch_MetBio, data = dat))
adjmat <- apply(ctlmat, 2, adj)
plsdat <- bind_cols(alc = log(dat$Qe_Alc/12 + 0.01), adjmat)

# Train PLS. Find the number of dimensions with lowest cross validation error
library(pls)
set.seed(111)
mod <- plsr(alc ~ ., ncomp = 10, data = plsdat, validation = "CV")
cv <- RMSEP(mod)
plot(cv, legendpos = "topright")

# Get no. components by lowest RMSEP, min, "one SE" and "permutation" methods
which.min(cv$val[estimate = "adjCV", , ]) - 1
selectNcomp(mod, method = "onesigma", plot = F)
selectNcomp(mod, method = "randomization", plot = T)

# Run PLS with 2 components and get coefficients, predict for mat
mod1 <- plsr(alc ~ ., data = plsdat, ncomp = 2)
coeff <- data.frame(value = round(coef(mod1)[ , 1, 1], 3))
crc1$M.alc <- exp(predict(mod1, mat)[,,1])

#fit3a <- clogit(update(base, ~. + M.alc), data = crc1) 
#OR for mediator 2.63 (1.43-4.81)


# Second mediator derived from EPIC pooled controls

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
plsdat <- bind_cols(alc = log(dat1$Qe_Alc/12 + 0.01), adjmat)

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
mod1 <- plsr(alc ~ ., data = plsdat, ncomp = 4)
coeff <- data.frame(value = round(coef(mod)[ , 1, 1], 3))

# Get cases for prediction, log2 and scale as for discovery matrix
crc1$M.alc1 <- exp(predict(mod1, data.frame(mat))[,,1])
