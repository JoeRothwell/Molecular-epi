# output tables and format for paper. The appropriate script should be run first
library(stargazer)
stargazer(mod, type = "html", out="C:/Your/Path/Name.html")

# Association calculated scores whole dataset (from list of models)
# Need to find a way to transpose

stargazer(modlist, type = "text", ci = T, apply.coef = exp, keep = "x", omit.stat = "all", flip = F,
          column.labels = scorecomp, out = "calc_scores.html", no.space = T, 
          digits = 2, single.row = T)
# notes = "Additionally adjusted for total energy intake, highest education level and 
          #smoking status.")

# Get most important compounds in PLS model for manuscript table 3
table3 <- bind_rows(table3a, table3b)
stargazer(table3, type = "text", summary = F, digits = 3, out = "Cmpd_importance.html")


# Association calculated and predicted scores 2 case-controls
ll2 <- list(fit7, fit8, fit1, fit2, fit5, fit6)
stargazer(ll2, type = "html", ci = T, apply.coef = exp, keep = c("score.2.comps", "Wcrf_C_Cal"), 
          omit.stat = "all", digits = 2,
          #column.labels = scorecomp, 
          out = "assoc_signature.html", no.space = T, single.row = T)

# Association metabolites and GRS
controls <- controls %>% select(cmpd, subclass, p.adj, meanfc) %>% 
               mutate(mean.fc = exp(meanfc)) %>% filter(p.adj < 0.05)

