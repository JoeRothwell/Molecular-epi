# Graphical abstract for gastrotenterology submission

estimate <- c(0.68, 0.53, 0.93)
conf.low <- c(0.56, 0.3, 0.86)
conf.high <- c(0.86, 0.95, 1)

labs <- c("      Endogenous metabolites", 
          "      Fatty acids", 
          "      WCRF/AICR score")

par(mar=c(3,4,0,2))
library(metafor)
forest(estimate, ci.lb = conf.low, ci.ub = conf.high, refline = 1, 
       slab = labs, #alabs = c(0.2,1.2),
       efac = c(0,0), psize = 1.5,
       bg = c("dodgerblue","grey", "grey"),
       header =c("Measure of Healthy Lifestyle", "Odds ratio per unit increase"),
       rows = c(1,2,5), ylim = c(0,9), alim = c(0.25, 1.25),
       xlim = c(-0.54, 1.6),
       xlab = NA, annosym = c(" (", "-", ")"),
       pch = 23)

text(-0.54, c(3,6), c("Metabolic signature of WCRF/AICR score", "Questionnaire-based assessment"), pos=4)

par("usr")


# Cut-down version
par(mar=c(3,4,0,2))
library(metafor)
forest(estimate, ci.lb = conf.low, ci.ub = conf.high, refline = 1, 
       efac = c(0,0), psize = 1.5,
       rows = c(1,2,5), ylim = c(0,6), alim = c(0.2, 1.2),
       steps = 6,
       slab = NA,
       header = F,
       top = 0,
       xlab = NA, annosym = c(" (", "-", ")"),
       pch = 18)

