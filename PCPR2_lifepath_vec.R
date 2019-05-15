#PCPR2_lifepath_vec

library(tidyverse)
library(MetabolAnalyze)
library(car)

# 1623 observations of 44 intensity variables. Looks scaled version of dat 4 and final prepared data
Xdata <- read_tsv("1507_XMetabolite_std_cpmg_E3N.txt") %>% filter(CODBMB != 11094738)
meta <- read_csv("Lifepath_meta.csv") %>% filter(CODBMB != 11094738)
alldata <- inner_join(meta, Xdata, by = "CODBMB")
# 1109473 is not in the metadata anyway

X_DataMatrixScaled <- alldata %>% select(`3Hydroxybutyrate`:Succinate) %>% as.matrix
Z_Meta <- alldata %>%
  select(WEEKS, PLACE, AGE, BMI, MENOPAUSE, FASTING, SMK, DIABETE, CENTTIMECat1, SAMPYEAR, STOCKTIME) %>%
  mutate_at(vars(-AGE, -BMI, -WEEKS, -STOCKTIME), as.factor)

# Center or scale the data, edit the parameter "pareto" or "unit"
#X_DataMatrixCentered = scale(X_DataMatrix, center = TRUE, scale = FALSE)
#X_DataMatrixScaled = scaling(X_DataMatrixCentered, type = "pareto")

Z_MetaRowN <- nrow(Z_Meta)
Z_MetaColN <- ncol(Z_Meta)
ColNames <- names(Z_Meta)

# get number of PCs for threshold
# Amount of variability desired to be explained
pct_threshold <- 0.8

# Get eigenvectors and eigenvalues
X_DataMatrixScaled_t <- t(X_DataMatrixScaled)
symMat <- X_DataMatrixScaled %*% X_DataMatrixScaled_t
eigenData      <- eigen(symMat)
eigenValues    <- eigenData$values
eigenVecMatrix <- eigenData$vectors
percents_PCs   <- eigenValues/sum(eigenValues)

# Get number of PCs required for threshold
my_counter_2 <- sum(1 - cumsum(rev(percents_PCs)) <= 0.8)
if(my_counter_2 > 3) pc_n <- my_counter_2 else pc_n <- 3

pc_data_matrix <- eigenVecMatrix[, 1:pc_n ]

DataCol <- Z_MetaColN + 1

#TotSumSq1 <- var(pc_data_matrix1) * (Z_MetaRowN - 1)
#fit1 <- lm(pc_data_matrix1 ~ ., data = Z_Meta)
#AnovaTab1   <- Anova(fit1, type=3, singular.ok = T)

# # Run a linear model with each eigenvector as the response
TotSumSq <- apply(pc_data_matrix, 2, var) * (Z_MetaRowN - 1)
fit      <- lm(pc_data_matrix ~ ., data = Z_Meta)
AnovaTab <- Anova(fit, type=3, singular.ok = T)
SSP      <- AnovaTab$SSP

# Extract sum of squares for each factor, removing intercept column
#type3mat0 <- sapply(SSP, "[", 1:pc_n)[, -1] NO!!!
# Need to take the diagonal of each factor matrix to get sums of squares
Residuals     <- diag(AnovaTab$SSPE)
RR            <- Residuals/TotSumSq

type3mat0 <- sapply(SSP, diag)[, -1]
type3mat  <- cbind(type3mat0, "SumSqResiduals" = Residuals)
ST_ResidualR2 <- cbind("ST_R2" = 1-RR, "ST_Residuals" = RR)

# Make partial R2 matrix and populate from typeIII matrix
partialR2mat <- type3mat[, -DataCol] / (type3mat[, -DataCol] + type3mat[, DataCol])

# Apply eigenvalues as weights
eigenValues <- eigenData$values[1: pc_n]
weight     <- eigenValues/sum(eigenValues) 
partialR2MatWtProp <- cbind(partialR2mat, ST_ResidualR2[, 1]) * weight
colnames(partialR2MatWtProp) <- NULL
pR2Sums <- colSums(partialR2MatWtProp) * 100

# Plot data
#par(mfrow = c(1,2))
bp <- barplot(pR2Sums, ylab = "Weighted Rpartial2", ylim = c(0, 25), xlab = "",
              col = "red", las=2, cex.main = 1, 
              main = paste("Sources of variation in", Z_MetaRowN, "BC cases and controls"),
              font.main = 1)
axis(1, at = bp, labels = c(ColNames, "R2"), cex.axis = 0.8, las = 2) 
rounded <- round(pR2Sums, 3)
text(bp, pR2Sums, labels = rounded, pos = 3, cex = 0.8)