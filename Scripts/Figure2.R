######## Figure2A PCoA ##############
#import data
table_otu = read.csv(file = "otu_relative_abundance.csv", header = TRUE, row.names = 1)
table_subCST = read.csv(file = "CST_group.csv", header = TRUE)
#
library(vegan)
library(ggplot2)
library(ggforce)
library(ade4)
library(gglayer)
library(ggrepel)
#
rownames(table_subCST) = table_subCST[["SampleID"]]
t_table_otu = as.data.frame(t(table_otu)) 
t_table_otu = t_table_otu[rownames(table_subCST), ]
PERMANOVA_test_res = adonis2(t_table_otu ~ subCST, data = table_subCST, permutations = 999, method="bray")
adonis_R2 = round(PERMANOVA_test_res$R2[1]*100, 3)
P_value = round(PERMANOVA_test_res$`Pr(>F)`[1], 3)
#
bray_dis = vegdist(t_table_otu, method="bray")
pcoa = dudi.pco(bray_dis, scan = FALSE, nf=2)
#
pcoa_eig = (pcoa$eig)[1:2] / sum(pcoa$eig)
#
sample_site = data.frame({pcoa$li})[1:2]
sample_site$names = rownames(sample_site)
names(sample_site)[1:2] = c('PCoA1', 'PCoA2') 
#
table_subCST = table_subCST[rownames(sample_site),]
merge_table = cbind(sample_site, table_subCST)
#
write.table(merge_table, file = "plot_data_PCoA.txt", row.names = FALSE, sep = "/t", quote = FALSE)
#
manual_order_1 =  c("V", "IV-C3", "IV-C1", "IV-B", "III-B", "III-A", "II", "I-B", "I-A")
merge_table[["subCST"]] = factor(merge_table[["subCST"]], levels = manual_order_1, ordered = T)
merge_table <- as.data.frame(merge_table)
#
pcoa_plot = ggplot(merge_table, aes(PCoA1, PCoA2,color=subCST))+
  theme_test()+
  geom_vline(xintercept = 0, color = 'gray', size = 0.25, linetype = "dashed")+ 
  geom_hline(yintercept = 0, color = 'gray', size = 0.25, linetype = "dashed")+
  geom_point(aes(color = subCST ), size = 2.0, alpha = 0.75)+
  scale_colour_manual(values = c("I-A"="#AD002A99", "I-B"="#f89588", "II"="#42B54099", "III-A"="#f6a34a", "III-B"="#f6ddb4", "IV-B"="#00468B99", "IV-C1"="#bbbdc0", "IV-C3"="#f8ac8c", "V"="#9AC9DB"))+
  stat_ellipse(aes(color = subCST), type = "norm", linetype = 2, level = 0.99, alpha = 1, show.legend = FALSE)+
  labs(x = paste('PCoA1: ', round(100 * pcoa_eig[1], 2), '%'), y = paste('PCoA2: ', round(100 * pcoa_eig[2], 2), '%'))+
  labs(title = paste0("adonis R2 = ", adonis_R2, "% ; ", "P value = ", P_value), col = "subCST")+
  theme(panel.grid = element_line(color = 'gray', linetype = 1, size = 0.1), 
        panel.background = element_rect(color = 'black', fill = 'transparent'))
pcoa_plot
ggsave(pcoa_plot, file = "PCoA_bray_curtis.pdf", width = 7.2, height =6.7) 
#
rm(list = ls())
dev.off()
#
#
######## Figure2B Stacked Bar ##############
library(tidyverse)
library(ggplot2)
library(gg.gap)
otu_table <- read.csv("species_top10_order.csv",header = T)
pivot_longer(otu_table,cols = 3:13)
long_otu_table <- pivot_longer(otu_table, cols = 3:13, names_to = "species", values_to = "relative_abundance") 
#
manual_order_1 =  c("others","Prevotella_timonensis","Sneathia_vaginalis","Streptococcus_anginosus","Prevotella_bivia","Fannyhessea_vaginae","Lactobacillus_hominis","Lactobacillus_jensenii","Gardnerella_vaginalis","Lactobacillus_crispatus","Lactobacillus_iners")
long_otu_table[["species"]] = factor(long_otu_table[["species"]], levels = manual_order_1, ordered = T)
long_otu_table <- as.data.frame(long_otu_table)
#
manual_order_2 <- as.vector(otu_table [,'SampleID'])
long_otu_table[[""]] = factor(long_otu_table[[""]], levels = manual_order_2, ordered = T)
long_otu_table <- as.data.frame(long_otu_table)
#
p <- ggplot(long_otu_table,aes(x=,y=relative_abundance,fill=species))+
  geom_bar(position = "fill", stat = "identity")+
  # 
  scale_fill_manual(values = c("Lactobacillus_iners"="#AD002A99","Lactobacillus_crispatus"="#42B54099","Gardnerella_vaginalis"="#f6a34a", "Lactobacillus_jensenii"="#f6ddb4","Lactobacillus_hominis"="#00468B99", "Fannyhessea_vaginae"="#bbbdc0","Prevotella_bivia"="#0099B499","Streptococcus_anginosus"="#f8ac8c","Sneathia_vaginalis"="#9AC9DB","Prevotella_timonensis"="#f89588","others"="#ADB6B699"))+
  # 
  theme_bw()+
  theme(panel.grid = element_blank(),
        panel.background = element_rect(fill = "#f3f4f4")+
        scale_x_discrete("")+
        theme(axis.text.x = element_blank(),
              axis.ticks.x = element_blank(),
              axis.line.x = element_blank()))
#
ggsave(p, file = "StackedBar_top10_order.pdf", width = 10, height =6.7)
#
rm(list = ls())
dev.off()
#
#
######## Figure2C Violin Plot ##############
#
table_Alpha_diversity = read.csv("alpha_diversity_index_group.csv", header = T)
rownames(table_Alpha_diversity) = table_Alpha_diversity$SampleID
#
alpha_index = colnames(table_Alpha_diversity)[2:8]
# 
table_SCT = read.csv("CST_group.csv", header = T)
rownames(table_SCT) = table_SCT$sampleID
#
table_combined = merge(x = table_SCT, y = table_Alpha_diversity, by = "sampleID")
#
library(tidyr)
library(survminer)
library(dplyr)
library(ggplot2)
library(ggpubr)
library(showtext)
library(sysfonts)
library(ggplot2)
library(tidyverse)
library(rstatix)
library(ggpubr)
library(ggtext)
#
for (index in alpha_index) { # index = "shannon"
  table_combined_index = table_combined[ , c("sampleID", "subCST", "score", "CST", index)]
  colnames(table_combined_index) =  c("sampleID", "subCST", "score", "CST", "value")
  table_combined_index[["value"]] = as.numeric(table_combined_index[["value"]])
  # 
  manual_order_1 =  c("V", "IV-C3", "IV-C1", "IV-B", "III-B", "III-A", "II", "I-B", "I-A")
  #
  table_combined_index$subCST = factor(table_combined_index$subCST, levels = manual_order_1, ordered = TRUE) 
  #
  data_2 = group_by(table_combined_index, subCST)
  data_3 = summarise(data_2, number = n(), median = round(median(value), 2), y_median = median(value)*1.5, y_number = median(value)*1.7+0.9, position_y_label =max(value), position_y =max(value)*0.995)
  data_3$subCST = factor(data_3$subCST, levels = manual_order_1, ordered = TRUE)
  #计算组间P值
  table_combined_index_p_val2 <- table_combined_index %>% 
    wilcox_test(formula = value ~ subCST) %>% 
    add_significance(p.col = 'p',cutpoints = c(0,0.001,0.01,0.05,1),symbols = c('***','**','*','ns')) %>% 
    add_xy_position()
  save_name_1 = paste0(index, "_", "subCST_", "diversity_p_value_1391.csv")
  #
  head(table_combined_index_p_val2)
  table_combined_index_p_val2_1 = table_combined_index_p_val2[ , -12]
  write.csv(table_combined_index_p_val2_1, file = save_name_1,row.names = F)
  #
  save_name_3 = paste0(index, "_", "subCST_", "diversity_violin_and_box.pdf")
  
  P3= ggplot(data = table_combined_index, aes(x = subCST, y = value, fill = subCST)) +
    geom_violin(data = table_combined_index, aes(x = subCST, y = value, fill = subCST), trim = FALSE, color="white") + 
    geom_boxplot(data = table_combined_index, aes(x = subCST, y = value), fill = "white",
                 outlier.shape = NA, notch = F, width = 0.05, position=position_dodge(0.9)) +
    geom_text(data = data_3, aes(x = subCST, y = y_median, label = median, group = subCST), position = position_dodge(0.9)) + 
    geom_text(data = data_3, aes(x = subCST, y = y_number, label = number, group = subCST), position = position_dodge(0.9)) +     
    theme_test()+ 
    scale_fill_manual(values = c("I-A"="#AD002A99", "I-B"="#f89588", "II"="#42B54099", "III-A"="#f6a34a", "III-B"="#f6ddb4", "IV-B"="#00468B99", "IV-C1"="#bbbdc0", "IV-C3"="#f8ac8c", "V"="#9AC9DB"))+
    theme(panel.grid.major = element_blank(),  
          panel.grid.minor = element_blank())+ 
    ylab(index)+xlab("") 
    coord_flip()
  
  ggsave(P3, file = save_name_3, width = 4, height = 6.7) 
}
# 
rm(list = ls())
dev.off()
#
#
######## Figure2D ggClusterNet ##############
#
##1、Calculate correlation coefficients between microbial abundances
# Detect and install dependent packages
local({r <- getOption("repos")  
r["CRAN"] <- "http://mirrors.tuna.tsinghua.edu.cn/CRAN/"   
r["BioC_mirror"] <- "https://mirrors.tuna.tsinghua.edu.cn/bioconductor"
options(repos=r)}) 

if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

a = rownames(installed.packages())

install_package <- c("psych", "reshape2", "igraph", "VennDiagram", "pheatmap")

for (i in install_package) {
  if (!i %in% a)
    BiocManager::install(i, update = F)
}

#
suppressPackageStartupMessages(library("psych"))
suppressPackageStartupMessages(library("reshape2"))
suppressPackageStartupMessages(library("igraph"))
suppressPackageStartupMessages(library("VennDiagram"))
suppressPackageStartupMessages(library("pheatmap"))
library(Hmisc)
#
for (i in c("CST1A", "CST1B", "CST3A", "CST3B", "CST4B")) { # i="CST1A"
  genus = read.csv(file = paste0(i,"_genus_counts_filter.csv"), header = TRUE, row.names = 1)
  #
  genus<-as.matrix(genus)
  dim(genus)
  head(genus)
  #
  genus_corr <- rcorr(t(genus), type = 'spearman')
  #
  r <- genus_corr$r
  #
  p <- genus_corr$P
  p <- p.adjust(p, method = 'BH')
  p[p>=0.05] <- -1
  p[p<0.05 & p>=0] <- 1
  p[p==-1] <- 0
  #
  z <- r * p
  diag(z) <- 0
  head(z)[1:6,1:6]
  #
  write.table(data.frame(z, check.names = FALSE), paste0(i,"_genus_corr.matrix_counts_filter.txt"), col.names = NA, sep = '\t', quote = FALSE)
  
  ##2、get network
  library(igraph)
  #
  g <- graph.adjacency(z, weighted = TRUE, mode = 'undirected')
  #
  g <- simplify(g)
  #
  E(g)$correlation <- E(g)$weight
  E(g)$weight <- abs(E(g)$weight)
  #
  tax = read.csv(file = "taxonomy_genus_filter_ggClusterNet.csv", header = TRUE, row.names = 1)
  tax <- tax[as.character(V(g)$name), ]
  
  V(g)$kingdom <- tax$Kingdom
  V(g)$phylum <- tax$Phylum
  V(g)$class <- tax$Class
  V(g)$order <- tax$Order
  V(g)$family <- tax$Family
  V(g)$genus <- tax$Genus
  #
  g
  plot(g)
  #
  adj_matrix <- as.matrix(get.adjacency(g, attr = 'correlation'))
  write.table(data.frame(adj_matrix, check.names = FALSE), paste0(i,"_network.adj_matrix_genus.txt"), col.names = NA, sep = '\t', quote = FALSE)
  #
  edge <- data.frame(as_edgelist(g))
  edge_list <- data.frame(
    source = edge[[1]],
    target = edge[[2]],
    weight = E(g)$weight,
    correlation = E(g)$correlation
  )
  #
  edge_list$cor <- 0 
  edge_list[which(edge_list[,4] < 0),"cor"] <- -1
  edge_list[which(edge_list[,4] > 0),"cor"] <- 1 
  edge_list$type <- "undirected" 
  head(edge_list)
  write.csv(edge_list, paste0(i,"_network.edge_list_filter_genus.csv"), row.names = FALSE, quote = F)
  
  # Node attribute list
  node_list <- data.frame(
    Id = names(V(g)),
    label = names(V(g)),
    kingdom = V(g)$kingdom,
    phylum = V(g)$phylum,
    class = V(g)$class,
    order = V(g)$order,
    family = V(g)$family,
    genus = V(g)$genus
  )
  head(node_list)
  # Calculate the weight of a node
  node.result <- data.frame(id = as.character(rownames(genus)), weight = rowMeans(genus))
  head(node.result)
  # Merge attributes and node weights by sample name (Id)
  node_list_combined = merge(x = node_list, y = node.result, by.x = "Id", by.y = "id")
  head(node_list_combined)
  write.csv(node_list_combined, paste0(i,"_network.node_list_filter_genus.csv"), row.names = FALSE)
  # The optional edge files and node files can be input to Gephi for drawing.
  #Output node connectivity file
  degree_CST4B = data.frame(genus = names(igraph::degree(g)), degree = igraph::degree(g))
  write.csv(degree_CST4B, paste0(i,"_network_node_degree.csv"), row.names = F, quote = F)
}
#
rm(list = ls())
dev.off()
}
#
rm(list = ls())
dev.off()
#
#
######## Figure2E Venn ##############
### Venn diagram shows shared edges
suppressPackageStartupMessages(library("VennDiagram"))
## Read 5 edge files and node files
edge_CST1A <- read.csv(file = "CST1A_network.edge_list_filter_genus.csv", header = T, row.names=NULL)
edge_CST1B <- read.csv(file = "CST1B_network.edge_list_filter_genus.csv", header = T, row.names=NULL)
edge_CST3A <- read.csv(file = "CST3A_network.edge_list_filter_genus.csv", header = T, row.names=NULL)
edge_CST3B <- read.csv(file = "CST3B_network.edge_list_filter_genus.csv", header = T, row.names=NULL)
edge_CST4B <- read.csv(file = "CST4B_network.edge_list_filter_genus.csv", header = T, row.names=NULL)
# Count the number of edges
num1 <- length(edge_CST1A [,1])
num2 <- length(edge_CST1B [,1])
num3 <- length(edge_CST3A [,1])
num4 <- length(edge_CST3B [,1])
num5 <- length(edge_CST4B [,1])
# Calculate intersection of edges
e1 <- data.frame(source = c(as.character(edge_CST1A[,1]),as.character(edge_CST1A[,2])), target = c(as.character(edge_CST1A[,2]),as.character(edge_CST1A[,1])), cor = c(as.character(edge_CST1A[,6]),as.character(edge_CST1A[,6])))
e2 <- data.frame(source = c(as.character(edge_CST1B[,1]),as.character(edge_CST1B[,2])), target = c(as.character(edge_CST1B[,2]),as.character(edge_CST1B[,1])), cor = c(as.character(edge_CST1B[,6]),as.character(edge_CST1B[,6])))
e3 <- data.frame(source = c(as.character(edge_CST3A[,1]),as.character(edge_CST3A[,2])), target = c(as.character(edge_CST3A[,2]),as.character(edge_CST3A[,1])), cor = c(as.character(edge_CST3A[,6]),as.character(edge_CST3A[,6])))
e4 <- data.frame(source = c(as.character(edge_CST3B[,1]),as.character(edge_CST3B[,2])), target = c(as.character(edge_CST3B[,2]),as.character(edge_CST3B[,1])), cor = c(as.character(edge_CST3B[,6]),as.character(edge_CST3B[,6])))
e5 <- data.frame(source = c(as.character(edge_CST4B[,1]),as.character(edge_CST4B[,2])), target = c(as.character(edge_CST4B[,2]),as.character(edge_CST4B[,1])), cor = c(as.character(edge_CST4B[,6]),as.character(edge_CST4B[,6])))
#
g1 <- paste(e1[,1],e1[,2],e1[,3])
g2 <- paste(e2[,1],e2[,2],e2[,3])
g3 <- paste(e3[,1],e3[,2],e3[,3])
g4 <- paste(e4[,1],e4[,2],e4[,3])
g5 <- paste(e5[,1],e5[,2],e5[,3])
#
x =list()
x =list(SCT1A = g1, SCT1B = g2, SCT3A = g3, SCT3B = g4 ,SCT4B = g5)
#
str(x)
FantasticFox1<-c("#d37a20","#dbcb09","#3a9cbc","#ADB6B699","#a30019")
color_list = FantasticFox1
#
fill_colors = color_list[1:length(x)]
#
venn <- venn.diagram(x,col="white",fill=fill_colors,
                      alpha= 0.4,   
                      lwd=.5,lty=1,      
                      cat.dist = c(0.2, 0.2, 0.2, 0.2,0.2),
                      filename=NULL,
                      cex=.5,
                      cat.cex=.5)
pdf("venn_CSTall_coedge.pdf",width = 5, height = 5)
grid.draw(venn)
dev.off()
#
#
######## Figure2F Heatmap##############
### Heatmap shows degree comparison of points in two network graphs
CST1A = "CST1A_network_node_degree.csv"
CST1B = "CST1B_network_node_degree.csv"
CST3A = "CST3A_network_node_degree.csv"
CST3B = "CST3B_network_node_degree.csv"
CST4B = "CST3B_network_node_degree.csv"

suppressPackageStartupMessages(library("pheatmap"))
#
degree.CST1A <- read.csv(CST1A, header = T)
degree.CST1B <- read.csv(CST1B, header = T)
degree.CST3A <- read.csv(CST3A, header = T)
degree.CST3B <- read.csv(CST3B, header = T)
degree.CST4B <- read.csv(CST4B, header = T)
#
CSTheat <- data.frame(CST1A=degree.CST1A[,2],CST1B=degree.CST1B[,2],CST3A=degree.CST3A[,2],CST3B=degree.CST3B[,2],CST4B=degree.CST4B[,2])
rownames(CSTheat) <- degree.CST1A[,1]
#绘图，自定义颜色
h = pheatmap(CSTheat, color = colorRampPalette(colors = c("#4292c6","#fec44f"))(100),
             cluster_rows = FALSE,
             cluster_cols = FALSE)
pdf("heatmap_CSTall_degree_compare.pdf",height = 5,width = 6)
h
dev.off()
rm(list = ls())
