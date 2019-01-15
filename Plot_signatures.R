# Get coefficients object from Metabolic_signature_WCRF.R
coefficients
df <- data.frame(as.list(coefficients)) %>% gather(Cmpd, VIP)
plot(coefficients)

# Get top and bottom deciles of compound coefficients
qs <- quantile(coefficients, probs = seq(0, 1, 0.05))
df1 <- df[df$VIP > qs[18], ]
df2 <- df[df$VIP < qs[4], ]
# Vector of colours for  plot points
vec <- ifelse(df$VIP > qs[18] | df$VIP < qs[4], "black", "grey")

#colvec <- cut(df$VIP, breaks = c(min(df$VIP), -1.5, 1.5, max(df$VIP)), 
    #include.lowest = T, labels = c("grey", "red", "blue"))

#col = I(brewer.pal(nlevels(colvec), name = "Dark2"))

# Now plot data, adding text
plot(coefficients, pch = 17, col=vec, xlab = "")
text(nrow(df) - nrow(df1):1, df1$VIP, df1$Cmpd, pos=2, cex = 0.6)
text(1:nrow(df2), df2$VIP, df2$Cmpd, pos=4, cex=0.6)
abline(a=0, b=0, lty = "dotted")


ggplot(df, aes(y = ind, x = 1, fill = values)) + geom_tile()
heatmap(as.matrix((coefficients)))

library(RColorBrewer)

        