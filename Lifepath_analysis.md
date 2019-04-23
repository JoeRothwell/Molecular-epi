Lifepath NMR metabolomics analysis
================

Exploratory analysis
--------------------

Read in metabolomics data from text file. There are 44 compounds measured for 1623 subjects, and subject identifiers in the first column, `CODBMB`.

``` r
library(tidyverse)
ints <- read_tsv("1507_XMetabolite_std_cpmg_E3N.txt")
dim(ints)
```

    ## [1] 1623   45

``` r
colnames(ints)
```

    ##  [1] "CODBMB"                   "3Hydroxybutyrate"        
    ##  [3] "Acetate"                  "Acetoacetate"            
    ##  [5] "Acetone"                  "Albumin"                 
    ##  [7] "Cholesterol"              "Choline"                 
    ##  [9] "cis-Aconitate"            "Creatine"                
    ## [11] "Creatinine"               "Glucose"                 
    ## [13] "Mannose"                  "Dimethylamine"           
    ## [15] "Ethanol"                  "Fatty acid"              
    ## [17] "Fatty acid (mainly LDL)"  "Fatty acid (mainly VLDL)"
    ## [19] "Formate"                  "Glutamine"               
    ## [21] "Glycerol"                 "Glycerophosphocholine"   
    ## [23] "Glycine"                  "Hypoxanthine"            
    ## [25] "Inosine"                  "Alanine"                 
    ## [27] "Aspartate"                "Glutamate"               
    ## [29] "Histidine"                "Isoleucine"              
    ## [31] "Lactate"                  "Lysine"                  
    ## [33] "Methionine"               "Ornithine"               
    ## [35] "Phenylalanine"            "Proline"                 
    ## [37] "Tyrosine"                 "Valine"                  
    ## [39] "Leucine"                  "Malonate"                
    ## [41] "Methanol"                 "NAC 1"                   
    ## [43] "NAC 2"                    "Pyruvate"                
    ## [45] "Succinate"

Intensity data are both positive and negative and most values in between -2 and 2. Data appear to be scaled and centered. We can plot the total intensity of each metabolite:

``` r
plot(colSums(ints[ , -1]), xlab = "Compound number", ylab = "Scaled intensity",
     pch = 19, col = "dodgerblue")
```

![](Lifepath_analysis_files/figure-markdown_github/totalintensity-1.png)

Total intensities are similar between compounds and normally distributed.We can now look at correlations between compounds:

``` r
cormat <- cor(ints)
colnames(cormat) <- NULL
library(corrplot)
```

    ## corrplot 0.84 loaded

``` r
corrplot(cormat, method = "square", tl.col = "black", tl.cex = 0.8)
```

![](Lifepath_analysis_files/figure-markdown_github/unnamed-chunk-2-1.png)

Of note is that the three fatty acid variables are correlated, as well as valine and leucine and NAC1 and NAC2. Next we can plot a PCA of the compound profiles:

``` r
ints <- ints[ , -1]
pca <- prcomp(ints, scale.=F)
library(pca3d)
pca2d(pca, title = "Metabolite profiles of 1623 subjects", xlab = "Score on PC1", ylab = "Score on PC2")
```

![](Lifepath_analysis_files/figure-markdown_github/unnamed-chunk-3-1.png)

There is one notable outlying profile in the bottom right-hand corner of the plot. This could be caused by a technical issue with the sample.
