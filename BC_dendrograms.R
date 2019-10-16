library(corrplot)
cormat <- cor(ints1)
colnames(cormat) <- NULL
corrplot(cormat, method = "square", tl.col = "black", tl.cex = 0.8)

dend1 <- cormat %>% dist %>% hclust %>% as.dendrogram

plot(dend0)
plot(dend1)

library(dendextend)
dend1 %>% set("labels_col", "blue") %>% hang.dendrogram() %>% 
  plot(main = "Change label's color")

dend1 %>% set("labels_cex", 1) %>% set("labels_col", value = c(3,4), k=2) %>% 
  plot(main = "Color labels \nper cluster")

dend1 %>% set("branches_k_color", k = 3) %>% plot(main = "Nice defaults")

library(circlize)


rownames(cormat) <- abbreviate(rownames(cormat), named = F)
dend <- cormat %>% dist %>% hclust %>% as.dendrogram %>%
  set("branches_k_color", k=3) %>% set("labels_cex", c(1)) #%>%
  #hang.dendrogram(hang_height = 0.3)

circos.par(gap.after = 0) #gap.degree for gaps between sectors

circlize_dendrogram(dend, dend_track_height = 0.85)
