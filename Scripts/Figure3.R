######## Figure3 combination chart ##############
# Correlation between genera relative abundance and trait data
library(correlation)
library(tidyr)
library(openxlsx)
library(readxl)
library(psych)
# import data
table_OTU = read.csv("genus_relative_abundance_filter.csv", header = T, row.names = 1)
#
group_table = read.csv("CST_group.csv", header = T, row.names = 1)
#
table_traits = read.csv("Traits_table.csv", header = T, row.names = 1)
table_traits = table_traits[rownames(group_table), ]
#
OTU_Traits_cor = corr.test(x = table_OTU, y = table_traits, method = "spearman", adjust = "fdr", normal=FALSE)
OTU_Traits_cor_value = as.data.frame(OTU_Traits_cor$r)
OTU_Traits_cor_P_value = as.data.frame(OTU_Traits_cor$p)
colnames(OTU_Traits_cor_P_value) = paste0("P_value_", colnames(OTU_Traits_cor_P_value))
#
merge_table = cbind(OTU_Traits_cor_value, OTU_Traits_cor_P_value)
write.csv(merge_table, "genus_merge_table_cor_P_value.csv")
##
#
table = merge_table[,c(12:22)]  # Extract data on the significant degree of correlation between phenotype and bacteria
list_Bacteria_all = list() # Create a list of all phenotypically related bacteria
list_Bacteria_0.05 = list() # Create a list of bacteria with significant (P<0.05) phenotypes
list_Bacteria_0.01 = list() # Create a list of bacteria with significant (P<0.01) phenotypes
list_Bacteria_0.001 = list() # Create a list of bacteria associated with a phenotype that is significantly (P<0.001)

#Cycling to extract significant bacteria
for (i in colnames(table)) {
  list_Bacteria_all[[i]] = rownames(table)
  list_Bacteria_0.05[[i]] = rownames(table[table[,i] < 0.05, ])
  list_Bacteria_0.01[[i]] = rownames(table[table[,i] < 0.01, ])
  list_Bacteria_0.001[[i]] = rownames(table[table[,i] < 0.001, ])
}
#
list_Bacteria = list(P_all = list_Bacteria_all, P_0.05 = list_Bacteria_0.05, P_0.01 = list_Bacteria_0.01, P_0.001 = list_Bacteria_0.001)
#
n="all"
#Combine the significant bacteria and extract the correlation data between bacteria and phenotypes
Union_table_Bacteria = unique(Reduce(union, list_Bacteria[[paste0("P_", n)]]))
table_Bacteria_with_significant_correlation = na.omit(merge_table[,c(1:11)][Union_table_Bacteria,])
table_Bacteria_with_significant_correlation_P_value = na.omit(merge_table[,c(12:22)][Union_table_Bacteria,])
colnames(table_Bacteria_with_significant_correlation_P_value) = colnames(table_Bacteria_with_significant_correlation)
# BiocManager::install("pacman")
library(pacman)
pacman::p_load(tidyverse,psych,reshape,ggtree,aplot)
#
cor <- as.matrix(table_Bacteria_with_significant_correlation) 
pvalue <- as.matrix(table_Bacteria_with_significant_correlation_P_value) 
#
myfun <- function(pval) {
  stars = ""
  if(pval <= 0.001)
    stars = "***"
  if(pval > 0.001 & pval <= 0.01)
    stars = "**"
  if(pval > 0.01 & pval <= 0.05)
    stars = "*"
  if(pval > 0.05 & pval <= 0.1)
    stars = ""
  stars
}
#
heatmap <- melt(cor) %>% rename(replace=c("X1"="Sample",
                                          "X2"="Group", 
                                          "value"="Cor")) %>%
  mutate(pvalue = melt(pvalue)[,3]) %>%  
  mutate(signif = sapply(pvalue, function(x) myfun(x)))
save_name = paste0("genus_cor_heatmap_", n, ".csv")
write.csv(heatmap, file = save_name)

#Cluster tree for plotting correlation and P-value plots
P_cor_phr <- hclust(dist(cor)) %>% ggtree(layout="rectangular", branch.length="none") # row clustering tree
#
P_cor <- ggplot(heatmap, aes(Group, Sample, col = Cor))+ 
  geom_tile(color="grey70", fill="white", size=0.5)+ 
  geom_point(aes(size = as.numeric(abs(cor))), shape=16)+ 
  geom_text(aes(label=signif),size=6,color="white", 
            hjust=0.5,vjust=0.7)+
  labs(x = NULL,y = NULL,color=NULL)+ 
  scale_color_viridis_c()+
  scale_y_discrete(position="right")+xlab(NULL) + ylab(NULL)+
  scale_x_discrete(position="top")+
  #
  theme(axis.ticks.x = element_blank(),
        axis.ticks.y = element_blank(),
        panel.border = element_rect(fill=NA,color="grey70",
                                    size=2, linetype="solid"))+
  theme(
    axis.text.x = element_blank(), 
    axis.text.y = element_blank(), 
    legend.text = element_text(face="plain",family = "Times",
                               colour = "black",size = 12))+
  labs(fill= "")+ #
  scale_size(range=c(1,10), guide="legend")+
  scale_colour_gradient2(high="#fb8500", mid="white", low="#007f5f")+ 
  guides(color = guide_colorbar(direction = "horizontal", 
                                reverse = F,barwidth = unit(4.5, "cm"), 
                                barheight = unit(0.5, "cm")))
##
##
microbial_name = unique(heatmap$Sample) # Extract microbial names from abundance heatmaps
table_OTU_abundance = read.csv("genus_relative_abundance_filter.csv", header = TRUE, row.names = 1)表
#The scale() function standardizes data by columns
table_microbial_name_OTU_abundance = table_OTU_abundance[colnames(table_OTU_abundance) %in% microbial_name, ]
table_microbial_name_OTU_abundance = scale(table_microbial_name_OTU_abundance, center = TRUE, scale = TRUE) # Standardize bacterial abundance in different samples
table_microbial_name_OTU_abundance = as.data.frame(t(table_microbial_name_OTU_abundance)) 
table_microbial_name_OTU_abundance[table_microbial_name_OTU_abundance > 1] = 1 
table_microbial_name_OTU_abundance[table_microbial_name_OTU_abundance < -1] = -1 
#
#
OTU_abundance_heatmap <- melt(as.matrix(table_microbial_name_OTU_abundance)) %>% rename(replace=c("X1"="Sample",
                                                                                                  "X2"="Group", 
                                                                                                  "value"="value"))

#Clustering tree for plotting abundance heatmap
table_microbial_name_OTU_abundance_tran = read.csv("genus_relative_abundance_filter.csv", header = TRUE, row.names = 1) 

#
P_OTU_abundance_hclust <- hclust(dist(table_microbial_name_OTU_abundance_tran, method = "maximum"), method = "complete")
#
P_OTU_abundance_phc <- ggtree(P_OTU_abundance_hclust) + layout_dendrogram() 

#Draw sample classification diagram
Strain_classification_heatmap_data  = read.csv("CST_group.csv", header = T)
#
Strain_classification_heatmap_data$subCST = factor(Strain_classification_heatmap_data$subCST,
                                                   levels = c("I-A", "I-B", "II", "III-A", "III-B","IV-B","IV-C0","IV-C1","IV-C3","V"), ordered = TRUE)
#
Strain_classification_heatmap_data = arrange(Strain_classification_heatmap_data, subCST)
#unique(Strain_classification_heatmap_data$subCST)
P_Strain_classification_heatmap <- ggplot(Strain_classification_heatmap_data, aes(x=SampleID, y=Group, fill=subCST))+
  geom_tile()+
  theme_minimal()+
  scale_y_discrete(position="right")+ 
  xlab(NULL)+ylab(NULL)+ 
  theme(panel.grid = element_blank())+  
  #
  scale_fill_manual(values = c("I-A"="#e07a5f", "I-B"="#f2cc8f", "II"="#3d405b", "III-A"="#81b29a", "III-B"="#f4f1de","IV-B"="#d7d7d8","IV-C0"="#6d6e70","IV-C1"="#c4dff6","IV-C3"="#56a0d3","V"="#0091cd"))+
  #
  theme_minimal()+
  xlab(NULL)+ylab(NULL)+
  theme(panel.grid = element_blank())+  
  theme(axis.ticks.x = element_blank(), 
        axis.text.x = element_blank(), 
        axis.ticks.y = element_blank(), 
        axis.text.y = element_blank(), 
  )+
  labs(fill= "subCST")
  
## Bacterial abundance heatmap
OTU_abundance_heatmap$Group = factor(OTU_abundance_heatmap$Group, levels = c(Strain_classification_heatmap_data$SampleID), ordered = TRUE)
P_OTU_abundance <- ggplot(OTU_abundance_heatmap,aes(x=Group, y=Sample, fill=value))+
  geom_tile()+
  theme_minimal()+
  scale_fill_viridis_c()+
  scale_y_discrete(position="right")+
  scale_x_discrete(position="top")+
  xlab(NULL)+ylab(NULL)+
  theme(panel.grid = element_blank())+
  scale_fill_gradient2(high="#f94144", mid="white", low="#277da1", midpoint = 0)+
  #
  theme(axis.ticks.x = element_blank(),
        axis.ticks.y=element_blank(),
        panel.border = element_rect(fill=NA,color="grey70",
                                    size=2, linetype="solid"))+
  theme(axis.text.x = element_blank(), 
        axis.text.y = element_text(family= "Times",face = "plain",colour = "black",size=12), 
        legend.text=element_text(face="plain",family = "Times",colour = "black",size = 12))+
  labs(fill= "")+ #
  guides(fill = guide_colorbar(direction = "horizontal", 
                               reverse = F,barwidth = unit(4.5, "cm"),
                               barheight = unit(0.5, "cm")))

## Bacterial abundance stacked column chart data
table_abundance_barplot_last = as.data.frame(matrix(ncol = 3))
colnames(table_abundance_barplot_last) = c("value", "sample_name", "Bacteria")
table_abundance_barplot = read.csv("genus_relative_abundance_filter.csv", header = T, row.names = 1)
table_abundance_barplot_delete_others = t(table_abundance_barplot[colnames(table_abundance_barplot) != "Others", ])
#Loop to summarize the data of the top 10 bacteria in each sample and calculate the value of others
for (sample in colnames(table_abundance_barplot_delete_others)) {
  table_abundance_barplot_delete_others = table_abundance_barplot_delete_others[order(table_abundance_barplot_delete_others[ ,sample], decreasing = T), ] 
  table_each_sample = as.data.frame(table_abundance_barplot_delete_others[ ,sample][1:10])
  table_each_sample$sample_name = sample
  table_each_sample$Bacteria = rownames(table_abundance_barplot_delete_others)[1:10]
  colnames(table_each_sample) = c("value", "sample_name", "Bacteria")
  table_each_sample[11, ] = c(1-sum(table_each_sample$value), sample, "others") 
  #
  table_abundance_barplot_last = rbind(table_abundance_barplot_last, table_each_sample) 
  table_abundance_barplot_last = na.omit(table_abundance_barplot_last) 
}
#The top 10 bacteria with average abundance in all samples at the genus level
table_abundance_barplot_delete_others_Top10 = data.frame(table_abundance_barplot_delete_others) 
table_abundance_barplot_delete_others_Top10$Average = rowMeans(table_abundance_barplot_delete_others_Top10)
#
table_abundance_barplot_delete_others_Top10 = table_abundance_barplot_delete_others_Top10[order(table_abundance_barplot_delete_others_Top10$Average, decreasing = T), ]
name_abundance_barplot_delete_others_Top10_Bacteria = rownames(table_abundance_barplot_delete_others_Top10)[1:10]
table_abundance_barplot_last_2 = table_abundance_barplot_delete_others[rownames(table_abundance_barplot_delete_others) %in% name_abundance_barplot_delete_others_Top10_Bacteria, ]
table_abundance_barplot_last_2['others',] = 1 - colSums(table_abundance_barplot_last_2)
table_abundance_barplot_last_2 = as.matrix(table_abundance_barplot_last_2)
#
table_abundance_barplot_last_2 = melt(table_abundance_barplot_last_2)
colnames(table_abundance_barplot_last_2) = c("Bacteria", "sample_name", "value")
#
table_abundance_barplot_last_2$Bacteria = factor(table_abundance_barplot_last_2$Bacteria,
                                                 levels = c(name_abundance_barplot_delete_others_Top10_Bacteria, "others"), ordered = TRUE)

# Draw a stacked column chart of the Top 10 bacterial species in each sample
## Stacked column chart of abundance of the top 10 bacterial genera at the genus level saved individually
fig_P_abundance_barplot_name = paste0("Top 10 accumulation histogram of bacterial abundance at genus level", ".pdf")
pdf(fig_P_abundance_barplot_name,  height = 5, width = 25)
ggplot(data = table_abundance_barplot_last_2)+
  geom_bar(aes(x = sample_name, y = value, fill = Bacteria), stat = "identity", position = "stack")+
  theme_minimal()+ 
  xlab(NULL)+ylab(NULL)+ 
  theme(panel.grid = element_blank())+  
  theme(axis.ticks.x = element_blank(), 
        axis.text.x = element_blank(), 
        axis.ticks.y = element_blank(), 
        axis.text.y = element_blank(), 
  )+
  scale_fill_manual(values = c('#cbc9ad', "#2c4251", '#d16666', "#b6c649", '#d58936', "#ffc482", '#dcb8cb', "#66999b", '#023618', "#fde935", '#69140e'))
dev.off()
  
## Abundance stacked histogram of the top 10 bacterial genera at the genus level
P_abundance_barplot = ggplot(data = table_abundance_barplot_last_2)+
  geom_bar(aes(x = sample_name, y = value, fill = Bacteria), stat = "identity", position = "stack")+
  theme_minimal()+ #
  xlab(NULL)+ylab(NULL)+ 
  theme(panel.grid = element_blank())+  
  theme(axis.ticks.x = element_blank(), 
        axis.text.x = element_blank(), 
        axis.ticks.y = element_blank(), 
        axis.text.y = element_blank(), 
  )+
  scale_fill_manual(values = c('#cbc9ad', "#2c4251", '#d16666', "#b6c649", '#d58936', "#ffc482", '#dcb8cb', "#66999b", '#023618', "#fde935", '#69140e'))+
  guides(
    fill = guide_colorbar(direction = "horizontal", reverse = F, barwidth = unit(4.5, "cm"), barheight = unit(0.5, "cm")) 
  )
#
P = P_OTU_abundance %>%
  insert_left(P_cor, width=0.25) %>% 
  insert_left(P_cor_phr, width=0.05) %>% 
  insert_top(P_abundance_barplot, height=0.30) %>% 
  insert_top(P_Strain_classification_heatmap, height=0.15) 


## Draw a heat map of each phenotype below the abundance heat map.
Phenotypic_heatmap_table = read.csv("Traits_table.csv", header = T, row.names = 1) 
#
Phenotypic_heatmap_table = Phenotypic_heatmap_table[rownames(group_table), ]
#
Phenotypic_heatmap_table = as.matrix(t(Phenotypic_heatmap_table))
# 
Phenotypic_heatmap_table <- melt(Phenotypic_heatmap_table) %>% rename(replace=c("X1"="Sample",
                                                                                "X2"="Group", 
                                                                                "value"="value"))
#
# After each "phenotype=" is run, lines 275-313 are run repeatedly, a total of 11 times.
phenotype = "AMH"
phenotype = "JS_T"
phenotype = "JS_P"
phenotype = "JS_PRL"
phenotype = "JS_E2"
phenotype = "JS_LH"
phenotype = "JS_FSH"
phenotype = "BMI"
phenotype = "weight"
phenotype = "height"
phenotype = "age"
#
Phenotypic_heatmap_data = Phenotypic_heatmap_table[Phenotypic_heatmap_table$Sample == phenotype, ]
# 
D_value = max(Phenotypic_heatmap_data$value) - min(Phenotypic_heatmap_data$value) 
table_frequency_distribution = as.data.frame(matrix(ncol = 3, nrow = 25)) 
colnames(table_frequency_distribution) = c("number", "start", "end") 
table_frequency_distribution$number = table(cut(Phenotypic_heatmap_data$value, c(seq(min(Phenotypic_heatmap_data$value), max(Phenotypic_heatmap_data$value), D_value/25)))) 
# Loop to calculate the start point and end point of the interval
for (number in c(1:25)) {
  table_frequency_distribution$start[number] = min(Phenotypic_heatmap_data$value) + number*D_value/25 - D_value/25 # Calculate the starting position of each interval
  table_frequency_distribution$end[number] = min(Phenotypic_heatmap_data$value) + number*D_value/25 # Calculate the ending position of each interval
}
# Determine the midpoint of the interval with the most frequency as the separation point
midpoint = table_frequency_distribution[table_frequency_distribution$number == max(table_frequency_distribution$number), ]$start/2 + table_frequency_distribution[table_frequency_distribution$number == max(table_frequency_distribution$number), ]$end/2
#
# 
P_phenotypic_heatmap <- ggplot(Phenotypic_heatmap_data,aes(x=Group, y=Sample, fill=value))+
  geom_tile()+
  theme_minimal()+ 
  scale_fill_viridis_c()+
  scale_y_discrete(position="right")+ 
  xlab(NULL)+ylab(NULL)+
  theme(panel.grid = element_blank())+ 
  #
  scale_fill_gradient2(high="#00296b", mid="#fefcfb", low="#ffd500",
                       midpoint = midpoint[1]
  )+ 
  #
  theme(axis.ticks.x = element_blank(), 
        axis.text.x = element_blank(),
        axis.ticks.y = element_blank(), 
        axis.text.y = element_text(family= "Times",face = "plain",colour = "black",size=12), 
        panel.border = element_rect(fill=NA, color="grey70",size=1, linetype="solid"))+
  guides(
    fill = guide_colorbar(title = phenotype, direction = "horizontal", reverse = F, barwidth = unit(7.5, "cm"), barheight = unit(0.5, "cm")) 
  ) 
#
P %>% insert_bottom(P_phenotypic_heatmap, height=0.06)
#
P = P %>% insert_bottom(P_phenotypic_heatmap, height=0.06)
  
## Add shannon index histogram
table_shannon_index = read.csv("alpha_diversity_index_group.csv", header = T)[ ,c(1,3)]
table_shannon_index$color = "color"
colnames(table_shannon_index)[1] <- "Sample_name"
table_shannon_index = table_shannon_index[table_shannon_index[["Sample_name"]] %in% rownames(group_table), ]
library(ggplot2)
P_shannon_index = ggplot2::ggplot(data = table_shannon_index)+
  geom_bar(aes(x = Sample_name, y = shannon, fill = color), stat = "identity")+
  theme_minimal()+
  xlab(NULL)+ylab("Shannon index")+ 
  theme(panel.grid = element_blank())+ 
  theme(axis.ticks.x = element_blank(), 
        axis.text.x = element_blank(),
        axis.ticks.y = element_line(), 
        axis.text.y = element_text(family= "Times",face = "plain",colour = "black",size=12)
  )+
  scale_fill_manual(values = c('#86bbd8'))+
  guides(fill = "none", 
  ) 
P = P %>% insert_bottom(P_shannon_index, height=0.24)
#
pdf("genus_cor_heatmap_composite_graphics.pdf", height = 25.5, width = 42)
P
dev.off()
