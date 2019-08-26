#Liver cancer/liver disease mortality case control. 
#Subsetting of samples, QCs, etc
library(readr)
library(dplyr)
library(stringr)
library(tidyr)
library(forcats)

#get distribution of features in samples
rawdata <- read.delim("Results table.txt", skip=4) 
pt      <- rawdata %>% select(contains("raw")) %>% t # ends_with doesn't work for some reason
allfeat <- rawdata$Compound
hist(colSums(pt > 1), breaks=100)

#write the peak table order to a csv
write.csv(colnames(rawdata), "peak table order.csv")

#make a factor for the type of observation: sample, QC or blank
type <- str_length(rownames(pt)) %>% as.factor %>%
  recode("29"="QC", "36"="QC", "40"="Sample", "39"="Sample", "38"="Sample", "30"="Blank")

#make a factor for batch: batch 1 or batch 2
batch <- str_sub(rownames(pt), 14L, 15L) %>% as.factor

logmat <- (log2(pt))

#make different matrix subsets (log or non log transformed)
sampmat <- logmat[type == "Sample", ]
noQCmat <- logmat[type != "QC", ]
QCmat   <- pt[type == "QC", ]
QCmatB1 <- pt[type == "QC" & batch == "04", ]
QCmatB2 <- pt[type == "QC" & batch == "14", ]
IDsamp  <- rownames(sampmat)
IDnoQC  <- rownames(noQCmat)
s.batch <- str_sub(IDsamp, 14L, 15L)
q.batch <- str_sub(IDnoQC, 14L, 15L)

#----------------------------------------------------------------------------------------------
#filtering of batch specific features
#find any features in at least 10% of batch 1 OR 10% of batch 2. First make a logical vector for batch
d1 <- colSums(sampmat[ s.batch == "04", ] != 0)
d2 <- colSums(sampmat[ s.batch == "14", ] != 0)

sum(d1 > 51 & d2 > 51) #4477 features satisfy criterion
fss <- which(d1 > 51 & d2 > 51) #vector of features to be retained

#filter original peak table, convert to matrix and transpose, add feature names
ptfilt <- t(as.matrix(pt[, fss]))
min(colSums(ptfilt > 1))
rownames(ptfilt) <- allfeat[fss]

#-----------------------------------------------------------------------------------------------
#filtering of blank specific features (background)
logmat1 <- t(log(ptfilt))
noQCfilt <- logmat1[, type != "QC" ]
medblanks <- apply(noQCfilt[, c(1:6, 518:523) ], 1, median)
medsample <- apply(noQCfilt[, c(7:517, 518:ncol(noQCfilt)) ], 1, median)
sum(medblanks > medsample) #1578 features satisfy criterion
fss1 <- which(medblanks > medsample)

#filter batch filtered peak table to give final peak table
ptfilt2 <- ptfilt[, -fss1]
hist(colSums(ptfilt2 > 1), breaks=100)
colnames(ptfilt2) <- colnames(ptfilt)[-fss1]

#------------------------------------------------------------------------------------------------
#CV histogram for features found in at least 90% of the QCs using QCmat
rsd.matrix <- function(x) {
  cm   <- colMeans(x, na.rm=T)
  sds  <- apply(x, 2, sd, na.rm=T)
  rsds <- (sds/cm)*100
  return(rsds)
}

ind <- colSums(QCmatB2 != 1) > 10
QCmat2 <- QCmatB2[ , ind ]

hist(rsd.matrix(QCmatB2), breaks=100)
plot(ecdf(rsd.matrix(QCmat2)))
#-----------------------------------------------------------------------------------------------

#Calculation of targeted RSDs
library(readr)
library(stringr)
batch1 <- read_csv("005_QC37.csv") %>% mutate(Batch = "Batch 1", Sequence = rank(File))
batch2 <- read_csv("010_QC38.csv") %>% mutate(Batch = "Batch 2", Sequence = rank(File))
#Need to group by batch to find ranks separately for each batch
#Dense rank gets the sequence of files with ties
all    <- bind_rows(batch1, batch2) %>% group_by(Batch) %>% 
          mutate("Area2" = Area/1000000, rank = dense_rank(File)) %>% 
          filter(!str_detect(Name, "Background"))

library(ggplot2)
library(grid)
ggplot(all, aes(x=rank, y=Area2, colour=Name, group=Name)) + geom_line() + 
  facet_grid(. ~ Batch) + theme_bw(base_size = 13) + 
  xlab("Pooled QC order") + 
  ylab("Peak area (millions of units)") +
  theme(legend.title=element_blank(), legend.key.height=unit(4, "mm"))

#-----------------------------------------------------------------------------------------------
#Scaling and PCA
#create metadata from text labels
df <- tbl_df(colnames(pt)) %>% separate(value, into=c("Cof", "Date", "type", "id"), sep="_") %>%
  select(Date, type) %>% separate("type", into=c("F1", "F2", "F3"), sep="-")

library(MetabolAnalyze)
logmat <- log2(ptfilt2)
scalemat <- scaling(logmat, type="pareto")
matpc    <- prcomp(scalemat, scale. = F)

pcs <- matpc$x
plot(pcs[, 1], pcs[, 2], pch=17, col="darkorange")
plot(pcs[, 2], pcs[, 3], pch=17, col="dodgerblue")

plot(pcs[, 1], pcs[, 2], col=type)
legend("topleft", c("sample", "QC", "blank"), col="type")
plot(pcs[, 2], pcs[, 3], col=type)

library(car)
library(rgl)
scatter3d(pcs[, 1], pcs[, 2], pcs[, 3], surface=F, id.n=nrow(matpc$x), point.col = "red")

