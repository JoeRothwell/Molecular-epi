## Liver cancer/liver disease mortality case control ##
## Filtering of redundant spectral features

#load packages (you may need to install them first)
library(dplyr)
library(tidyr)
library(stringr)

#read in metabolomics data and subset intensity matrix (pt = peak table)
rawdata <- read.delim("Unfiltered peak table.txt", skip=4) 
pt      <- rawdata %>% select(contains("raw")) %>% t

#store the feature names in a variable
allfeat <- rawdata$Compound

#write the peak table sample order to a csv
#write.csv(colnames(rawdata), "peak table order.csv")

#make a factor for the type of observation from row names: sample, QC or blank
type <- str_length(rownames(pt)) %>% as.factor %>%
  recode("29"="QC", "36"="QC", "40"="Sample", "39"="Sample", "38"="Sample", "30"="Blank")

#make a factor for batch: 1 or 2
batch <- str_sub(rownames(pt), 14L, 15L) %>% as.factor

#log2 transform intensities
logmat <- log2(pt)

#make different matrix subsets. Samples only:
sampmat <- logmat[type == "Sample", ]  
#samples and blanks only:
noQCmat <- logmat[type != "QC", ]  

#get sample names for two subsets
IDsamp  <- rownames(sampmat) #Get names for sample only matrix
IDnoQC  <- rownames(noQCmat) #Get names for sample and blank only matrix

s.batch <- str_sub(IDsamp, 14L, 15L) #Get batch numbers for sample only matrix
q.batch <- str_sub(IDnoQC, 14L, 15L) #Get batch numbers for sample and blank only matrix

#----------------------------------------------------------------------------------------------
#filtering of batch specific features
#find any features in at least 10% of batch 1 OR 10% of batch 2 (51 samples) (or change as desired)
#Make a logical vector for batch
d1 <- colSums(sampmat[ s.batch == "04", ] != 0)
d2 <- colSums(sampmat[ s.batch == "14", ] != 0)

sum(d1 > 51 & d2 > 51)          #4477 features satisfy criterion
fss <- which(d1 > 51 & d2 > 51) #vector of features to be retained

#filter original peak table by retained features and NCI samples only, convert to matrix and transpose
ptfilt <- as.matrix(pt[type=="Sample", fss])
#add feature names
colnames(ptfilt) <- allfeat[fss]

#check dimensions of final dataset, should be 1023 obs of 4477 variables
dim(ptfilt)

#save as an R object or write to .csv
save(ptfilt, file="Filtered peak table.Rdata")
write.csv(ptfilt, file="Filtered peak table.csv")

#to load it:
load("Filtered peak table no QCs.Rdata")

#-----------------------------------------------------------------------------------------------
#optional filtering of blank specific features (background)

#get filtered table with samples and blanks only
ptfilt2 <- as.matrix(pt[, fss])
noQCfilt <- ptfilt2[ type != "QC",  ]

#get filtered table with samples and blanks only
blanks   <- noQCfilt[ c(1:6, 518:523), ]
samples  <- noQCfilt[ c(7:517, 524:nrow(noQCfilt)), ]

#Get median intensity vectors of blanks (n=12) and samples (n=1023)
medblanks <- apply(blanks, 2, median)
medsample <- apply(samples, 2, median)

#Set criteria for blank filtering (change as desired)
sum(medblanks > medsample) #1598 features satisfy criterion
fss1 <- which(medblanks > medsample)

#filter batch filtered peak table to give final peak table
ptfilt2 <- ptfilt[, -fss1]
colnames(ptfilt2) <- colnames(ptfilt)[-fss1]

#export to file as above




