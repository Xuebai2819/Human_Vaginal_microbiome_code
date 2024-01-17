######## Figure1B Traits Correlation heatmap ##############
# import data
Traits_table <- read.csv("Traits_table.csv", header = T, row.names =1)
#
library(corrplot)
# Calculate correlation
cor_matrix <- cor(Traits_table, method = "spearman")
# 
write.csv(cor_matrix, "Traits_spearman_cor.csv")
# matrix of the p-value of the correlation
p.mat <- cor.mtest(Traits_table, method = "spearman")$p
# 
write.csv(p.mat, "Traits_spearman_cor_P.csv")
# 
pdf(file = "Traits_cor_heatmap.pdf", height = 5, width = 5)
col <- colorRampPalette(c("#BB4444", "#EE9988", "#FFFFFF", "#77AADD", "#4477AA"))
corrplot(cor_matrix, method = "circle", col = col(200),
         type = "lower", order = "original", number.cex = .7,
         addCoef.col = NULL, # Add coefficient of correlation
         tl.col = "black", tl.srt = 45, # Text label color and rotation
         # Combine with significance
         p.mat = p.mat, sig.level = c(.001, .01, .05), insig = "label_sig",
         pch.cex = .9, pch.col = "black",
         # hide correlation coefficient on the principal diagonal
         diag = FALSE)
dev.off()

# 
rm(list = ls())
dev.off()
#
#
######## Figure1C genus Correlation heatmap ##############
# import data
abundance_table <- read.csv("genus_relative_abundance_filter.csv", header = T, row.names = 1)
# 
group_table = read.table("CST_group.csv", header = T, row.names = 1)
# 
abundance_table = abundance_table[rownames(group_table), ]
# 加载包
library(corrplot)
# 
cor_matrix <- cor(abundance_table, method = "spearman")
#
write.csv(cor_matrix, "Genus_spearman_cor.csv", row.names = T)
# matrix of the p-value of the correlation
p.mat <- cor.mtest(abundance_table, method = "spearman")$p
write.csv(p.mat, "Genus_spearman_cor_P.csv", row.names = T)
# 
pdf(file = "Genus_cor_heatmap.pdf", height = 5, width = 5)
col <- colorRampPalette(c("#BB4444", "#EE9988", "#FFFFFF", "#77AADD", "#4477AA"))
corrplot(cor_matrix, method = "circle", col = col(200),
         type = "lower", order = "original", number.cex = 0.6,
         addCoef.col = NULL, # Add coefficient of correlation，
         tl.col = "black", tl.srt = 45, tl.cex = 0.5, # Text label color and rotation, size
         # Combine with significance
         p.mat = p.mat, sig.level = c(0.001, 0.01, 0.05), insig = "label_sig",
         pch.cex = 0.5, pch.col = "black",
         # hide correlation coefficient on the principal diagonal
         diag = FALSE)
dev.off()

# 
rm(list = ls())
dev.off()
