library(DESeq2)
library(ggplot2)
library(kohonen)
library(tsne)
library(umap)
library(pheatmap)
library(EnhancedVolcano)
library(org.Hs.eg.db)
library(clusterProfiler)
library(AnnotationDbi)
library(ComplexHeatmap)
library(ggstance)
library(forcats)
library(esquisse)
library(dbscan)
library(factoextra)
library(colorRamp2)
library(gplots)
library(corrplot)
library(Rtsne)

raw_counts = read.csv('../data/gene_expression/prepared_data/raw_count_BRCA.csv', sep = ',')
clinical = read.csv('../data/gene_expression/prepared_data/clinical_BRCA.csv', sep = ',')

#clinical_t = as.data.frame(t(clinical))

genes = raw_counts[,1]

for (i in 1:length(genes)){
  if (genes[i] == "?"){
    genes[i] = paste('?',i,sep = '_')
  }
  if (genes[i] == "SLC35E2"){
    genes[i] = paste(genes[i],i,sep = '_')
  }
}

duplicates <- genes[duplicated(genes)]
if (length(duplicates) > 0) {
  print("Duplicates found:")
  print(duplicates)
} else {
  print("No duplicates found.")
}

raw_counts = raw_counts[,-1]
rownames(raw_counts) = genes

col_names = colnames(raw_counts)
rownames(clinical) = col_names
clinical = clinical[,-1]



hrdgenes<-read.gmt("../data/HRDgenes.gmt")
kegg_symbols = read.gmt('../data/c2.cp.kegg.v7.5.1.symbols.gmt')
hallmark = read.gmt('../data/h.all.v7.5.1.symbols.gmt')


write.csv(raw_counts,'../data/raw_counts_deseq.csv', row.names = TRUE)


#### Functions ####


# Creates Volcano plots 
volcano_plotter = function(res, hrd_genes, reference_type, compare_type, pathway){
  res_df = as.data.frame(res)
  res_df$symbol = rownames(res_df)
  
  hrdgenes_all <- subset(res_df, rownames(res_df) %in% hrdgenes$gene)
  hrdgenes_all <- subset(hrdgenes_all, hrdgenes_all$log2FoldChange<0.58 | hrdgenes_all$log2FoldChange> -0.58 )
  hrdgenes_significant <- subset(res_df, rownames(res_df) %in% hrdgenes$gene & res_df$padj<=0.05 )
  hrdgenes_significant <- subset(hrdgenes_significant, hrdgenes_significant$log2FoldChange>=0.58 | hrdgenes_significant$log2FoldChange<= -0.58 )
  sig_genes_label<-subset(hrdgenes_significant)
  
  volcano_plot = ggplot(res_df, aes(x = log2FoldChange, y = -log10(padj))) +
    geom_point(color="grey87") +
    ggtitle("Genes belonging to the HR pathway") +
    theme_bw() +
    geom_point(data= hrdgenes_all, col="dodgerblue2") +
    geom_point(data= hrdgenes_significant, col="red") +
    geom_text_repel(data = sig_genes_label,
                    aes(x=log2FoldChange,
                        y=-log10(padj),label= symbol),
                    max.overlaps = 40) +
    theme(legend.position = "none") +
    scale_x_continuous(name = paste("log2(fold change) ",compare_type," vs ", reference_type, sep = '')) +
    scale_y_continuous(name = "-log10 p-value") +
    geom_hline(yintercept = -log10(0.05), linetype="dashed") +
    geom_vline(xintercept = 0.58, linetype="dashed") +
    geom_vline(xintercept = -0.58, linetype="dashed")
  
  print(volcano_plot)
  
  ggsave(paste(pathway,'volcanoplot_',compare_type,'_',reference_type,'.png',sep = ''),volcano_plot, width = 8.63, height = 5.71)
}



## Kegg and Hallmark GSEA
kegg_hallmark_hrdgenes_gsea = function(res, hrdgenes, kegg_symbols, hallmark, reference_type, compare_type, pathway){
  
  gene_list = res$stat
  res$symbol = rownames(res)
  names(gene_list) = make.names(res$symbol, unique = T)
  gene_list = gene_list[order(gene_list, decreasing = T)]
  
  gene_names = make.names(res$symbol, unique = T)
  
  hrdgenes_gsea = GSEA(gene_list, TERM2GENE = hrdgenes,
                       pvalueCutoff = 0.1,
                       eps=1e-50,
                       seed=T)
  
  kegg_gsea = GSEA(gene_list, TERM2GENE = kegg_symbols,
                   pvalueCutoff = 0.1,
                   eps=1e-50,
                   seed=T)
  
  hallmark_gsea = GSEA(gene_list, TERM2GENE = hallmark,
                       pvalueCutoff = 0.1,
                       eps=1e-50,
                       seed=T)
  
  kegg_gsea_df = as.data.frame(kegg_gsea)
  hrdgenes_gsea_df = as.data.frame(hrdgenes_gsea)
  hallmark_gsea_df = as.data.frame(hallmark_gsea)
  
  y <- mutate(kegg_gsea, ordering = abs(NES)) %>%
    arrange(desc(ordering))
  n <- 15
  y_bar <- group_by(y, sign(NES)) %>%
    slice(1:n)
  kegg_pathways = ggplot(y_bar, aes(NES, fct_reorder(Description, NES), fill = qvalue), showCategory=(n*2)) +
    geom_barh(stat='identity') +
    scale_fill_continuous(low='red', high='blue', guide=guide_colorbar(reverse=TRUE)) +
    theme_minimal() + ylab(NULL) + ggtitle(paste("KEGG GSEA cell models, ", compare_type,' vs ', reference_type, sep = ''))
  
  y <- mutate(hallmark_gsea, ordering = abs(NES)) %>%
    arrange(desc(ordering))
  n <- 15
  y_bar <- group_by(y, sign(NES)) %>%
    slice(1:n)
  hallmark_pathways = ggplot(y_bar, aes(NES, fct_reorder(Description, NES), fill = qvalue), showCategory=(n*2)) +
    geom_barh(stat='identity') +
    scale_fill_continuous(low='red', high='blue', guide=guide_colorbar(reverse=TRUE)) +
    theme_minimal() + ylab(NULL) + ggtitle(paste("Hallmark GSEA cell models, ", compare_type,' vs ', reference_type, sep = ''))
  
  print(kegg_pathways)
  print(hallmark_pathways)
  
  ggsave(paste(pathway,'kegg_pathways_',compare_type,'_vs_',reference_type,'.png',sep = ''),kegg_pathways, width = 8.63, height = 5.71)
  ggsave(paste(pathway,'hallmark_pathways_',compare_type,'_vs_',reference_type,'.png',sep = ''),hallmark_pathways, width = 8.63, height = 5.71)
  
  return (list(hrdgenes_gsea = hrdgenes_gsea_df, kegg_gsea = kegg_gsea_df, hallmark_gsea = hallmark_gsea_df))
}







#### Differential expression analysis of each subtypes versus each other and PCA ####

dds <- DESeqDataSetFromMatrix(countData = raw_counts,
                              colData = clinical,
                              design= ~ Subtype)

keep <- rowSums(counts(dds)) >= 10
dds <- dds[keep,]

dds <- DESeq(dds)
resultsNames(dds)



vsd <- vst(dds, blind=FALSE)

#res_all = assay(vsd)
#write.csv(res_all,'../data/res_all.csv', row.names = TRUE)


plotPCA(vsd, intgroup=c("Subtype"))



pcaData <- plotPCA(vsd, intgroup=c("Subtype", 'HRD_sum'), returnData=TRUE)
percentVar <- round(100 * attr(pcaData, "percentVar"))

pca_sub_hrd = ggplot(pcaData, aes(PC1, PC2, color= as.numeric(HRD_sum), shape=Subtype)) +
  geom_point(size=3) +
  xlab(paste0("PC1: ",percentVar[1],"% variance")) +
  ylab(paste0("PC2: ",percentVar[2],"% variance")) +
  scale_color_gradient(low = alpha("white", 0.7), high = "darkred")+
  coord_fixed() +
  labs(color = "HRD sum")

pca_hrd = ggplot(pcaData, aes(PC1, PC2, color= as.numeric(HRD_sum))) +
  geom_point(size=3) +
  xlab(paste0("PC1: ",percentVar[1],"% variance")) +
  ylab(paste0("PC2: ",percentVar[2],"% variance")) +
  scale_color_gradient(low = alpha("gray", 0.7), high = "darkred")+
  coord_fixed() +
  labs(color = "HRD sum")

pca_subtype = ggplot(pcaData, aes(PC1, PC2, color= Subtype)) +
  geom_point(size=3) +
  xlab(paste0("PC1: ",percentVar[1],"% variance")) +
  ylab(paste0("PC2: ",percentVar[2],"% variance")) +
  coord_fixed()



ggsave('../data/figures/BRCA_cohort/PCA_subtypes.png',pca_subtype, width = 8.63, height = 5.71)
ggsave('../data/figures/BRCA_cohort/PCA_HRD.png',pca_hrd, width = 8.63, height = 5.71)
ggsave('../data/figures/BRCA_cohort/PCA_subtypes_HRD.png',pca_sub_hrd, width = 8.63, height = 5.71)




#### Differential Expression analysis, Basal versus other ####


dds_2 <- DESeqDataSetFromMatrix(countData = raw_counts,
                              colData = clinical,
                              design= ~ De)

keep <- rowSums(counts(dds_2)) >= 10
dds_2 <- dds_2[keep,]

dds_2$De = relevel(dds_2$De, ref = 'other')

dds_2 <- DESeq(dds_2)

resultsNames(dds_2)
res_2 = results(dds_2, alpha = 0.01)


resOrdered <- res_2[order(res_2$padj),]

reslfc <- lfcShrink(dds_2,coef='De_Basal_vs_other', type = 'normal')

plotMA(resOrdered)

EnhancedVolcano(resOrdered,
                lab = rownames(resOrdered),
                x = 'log2FoldChange',
                y = 'pvalue',
                FCcutoff = 2,
                pCutoff = 0.01,
                title = 'BRCA Basal versus others')






resOrdered_df = as.data.frame(resOrdered)
resOrdered_df = resOrdered_df[rownames(resOrdered_df) %in% hrdgenes$gene,]
resOrdered_df = resOrdered_df[order(resOrdered_df$padj),]

genes_sig = rownames(resOrdered_df)[1:20]

vst <- vst(dds_2, blind=FALSE)
vst <- assay(vst)
vst <- as.data.frame(vst)
vst_sig <- vst[rownames(vst) %in% hrdgenes$gene,]
tvst <- t(vst_sig)
rownames(tvst) <- NULL

heat <- t(scale(tvst))

#scaled_tvst <- 2 * (tvst - min(tvst)) / (max(tvst) - min(tvst))

# Rescale the data to a range of -1 to 1
#heat <- scaled_tvst - 1

#heat = t(heat)

colnames(heat) = paste0(dds_2@colData$De)

type = dds_2@colData$De
col_order = colnames(heat)


heat = heat[rownames(heat) %in% genes_sig,]

ha = HeatmapAnnotation(Subtype = colnames(heat), col = list(Subtype = c('Basal' = 'darkblue', 'other' = 'lightblue')))

hrds_values = dds_2@colData$HRD_sum
loh_values = dds_2@colData$LOH
lst_values = dds_2@colData$LST
tai_values = dds_2@colData$TAI

hrd_annotation = HeatmapAnnotation(HRDs = hrds_values, col = list(HRDs = colorRamp2(c(min(hrds_values),max(hrds_values)),c('beige','saddlebrown'))))
loh_annotation = HeatmapAnnotation(LOH = loh_values, col = list(LOH = colorRamp2(c(min(loh_values),max(loh_values)),c('lightgreen', 'darkgreen'))))
lst_annotation = HeatmapAnnotation(LST = lst_values, col = list(LST = colorRamp2(c(min(lst_values),max(lst_values)),c('lightblue','darkblue'))))
tai_annotation = HeatmapAnnotation(TAI = tai_values, col = list(TAI = colorRamp2(c(min(tai_values),max(tai_values)),c('lightcyan', 'darkcyan'))))

#png(file='../data/figures/BRCA_cohort/heatmap.png')
hm = Heatmap(heat,show_column_names = FALSE, row_names_gp = gpar(fontsize = 4), bottom_annotation = ha, column_km = 4, name = 'z-score') %v% 
  hrd_annotation %v% loh_annotation %v% lst_annotation %v% tai_annotation

draw(hm)
#dev.off()

ggsave('../data/figures/BRCA_cohort/heatmap.png',hm, width = 8.63, height = 5.71)

Heatmap(heat)

Heatmap(heat,
        row_names_gp = gpar(fontsize = 8))



# GSEA using hallmark and kegg terms (second try gmt files from Alicia)

gene_list = res_2$stat
res_2$symbol = rownames(res_2)
names(gene_list) = make.names(res_2$symbol, unique = T)
gene_list = gene_list[order(gene_list, decreasing = T)]

gene_names = make.names(res_2$symbol, unique = T)

hrdgenes_gsea = GSEA(gene_list, TERM2GENE = hrdgenes,
                     pvalueCutoff = 0.1,
                    eps=1e-50,
                    seed=T)

kegg_gsea = GSEA(gene_list, TERM2GENE = kegg_symbols,
                 pvalueCutoff = 0.1,
                 eps=1e-50,
                 seed=T)

hallmark_gsea = GSEA(gene_list, TERM2GENE = hallmark,
                     pvalueCutoff = 0.1,
                     eps=1e-50,
                     seed=T)

hrdgenes_enrich = enricher(gene_names, TERM2GENE = hrdgenes)
kegg_enrich = enricher(gene_names, TERM2GENE = kegg_symbols)
hallmark_enrich = enricher(gene_names, TERM2GENE = hallmark)

kegg_gsea_df = as.data.frame(kegg_gsea)
hrdgenes_gsea_df = as.data.frame(hrdgenes_gsea)
hallmark_gsea_df = as.data.frame(hallmark_gsea)

barplot(kegg_enrich)



y <- mutate(kegg_gsea, ordering = abs(NES)) %>%
  arrange(desc(ordering))
n <- 15
y_bar <- group_by(y, sign(NES)) %>%
  slice(1:n)
kegg_pathways = ggplot(y_bar, aes(NES, fct_reorder(Description, NES), fill = qvalue), showCategory=(n*2)) +
  geom_barh(stat='identity') +
  scale_fill_continuous(low='red', high='blue', guide=guide_colorbar(reverse=TRUE)) +
  theme_minimal() + ylab(NULL) + ggtitle("KEGG GSEA cell models")

y <- mutate(hallmark_gsea, ordering = abs(NES)) %>%
  arrange(desc(ordering))
n <- 15
y_bar <- group_by(y, sign(NES)) %>%
  slice(1:n)
hallmark_pathways = ggplot(y_bar, aes(NES, fct_reorder(Description, NES), fill = qvalue), showCategory=(n*2)) +
  geom_barh(stat='identity') +
  scale_fill_continuous(low='red', high='blue', guide=guide_colorbar(reverse=TRUE)) +
  theme_minimal() + ylab(NULL) + ggtitle("Hallmark GSEA cell models")

ggsave('../data/figures/BRCA_cohort/kegg_pathways.png',kegg_pathways, width = 8.63, height = 5.71)
ggsave('../data/figures/BRCA_cohort/hallmark_pathways.png',hallmark_pathways, width = 8.63, height = 5.71)


# GSEA KEGG (First try, did not really work)

resOrdered_lfc = res_2[order(-res_2$log2FoldChange),]

resOrdered_lfc$entrez <- mapIds(x = org.Hs.eg.db,
                         keys = row.names(resOrdered_lfc),
                         column = "ENTREZID",
                         keytype = "SYMBOL",
                         multiVals = "first")

gene_list = resOrdered_lfc$log2FoldChange

names(gene_list) = resOrdered_lfc$entrez


kk2 <- gseKEGG(geneList     = gene_list,
               organism     = 'hsa',
               keyType = 'ncbi-geneid',
               pvalueCutoff = 0.05,
               verbose      = FALSE)

df_kegg = as.data.frame(kk2)

pathview(gene.data = gene_list, 
         pathway.id = "03400", 
         species = "hsa")


# GSEA GO (First try also not really want we were looking for)

resOrdered_lfc = res_2[order(-res_2$log2FoldChange),]

gene_list = resOrdered_lfc$log2FoldChange

names(gene_list) = rownames(resOrdered_lfc)



gse = gseGO(gene_list,
            ont = 'BP',
            keyType = 'SYMBOL',
            OrgDb = 'org.Hs.eg.db',
            eps = 0)


df_gse = as.data.frame(gse)

# Volcano plot of result highlighting HRD genes (from Alicia)

resOrdered_df = as.data.frame(res_2)
resOrdered_df$symbol = rownames(resOrdered_df)

hrdgenes_all <- subset(resOrdered_df, rownames(resOrdered_df) %in% hrdgenes$gene)
hrdgenes_all <- subset(hrdgenes_all, hrdgenes_all$log2FoldChange<0.58 | hrdgenes_all$log2FoldChange> -0.58 )
hrdgenes_significant <- subset(resOrdered_df, rownames(resOrdered_df) %in% hrdgenes$gene & resOrdered_df$padj<=0.05 )
hrdgenes_significant <- subset(hrdgenes_significant, hrdgenes_significant$log2FoldChange>=0.58 | hrdgenes_significant$log2FoldChange<= -0.58 )
sig_genes_label<-subset(hrdgenes_significant)

volcano_basal_other = ggplot(resOrdered_df, aes(x = log2FoldChange, y = -log10(padj))) +
  geom_point(color="grey87") +
  ggtitle("Genes belonging to the HR pathway") +
  theme_bw() +
  geom_point(data= hrdgenes_all, col="dodgerblue2") +
  geom_point(data= hrdgenes_significant, col="red") +
  geom_text_repel(data = sig_genes_label,
                  aes(x=log2FoldChange,
                      y=-log10(padj),label= symbol),
                  max.overlaps = 40) +
  theme(legend.position = "none") +
  scale_x_continuous(name = "log2(fold change) Basal vs other") +
  scale_y_continuous(name = "-log10 p-value") +
  geom_hline(yintercept = -log10(0.05), linetype="dashed") +
  geom_vline(xintercept = 0.58, linetype="dashed") +
  geom_vline(xintercept = -0.58, linetype="dashed")

ggsave('../data/figures/BRCA_cohort/basal_other_volcano.png',volcano_basal_other, width = 8.63, height = 5.71)



#sub_df_clinical = clinical_t[,c('Subtype','HRD_sum')]
sub_df_clinical = clinical[,c('Subtype','HRD_sum')]








# Violin plot

subtype_violin = ggplot(sub_df_clinical, aes(x = Subtype, y = as.numeric(HRD_sum), fill = Subtype)) +
  geom_violin() + geom_boxplot(width = 0.2)+
  labs(x = "Subtype", y = "HRD Score")

ggsave('../data/figures/BRCA_cohort/subtypes_violin.png',subtype_violin, width = 8.63, height = 5.71)

# Correlation tests

kruskal.test(HRD_sum ~ Subtype, data = clinical)

pairwise.wilcox.test(as.numeric(clinical$HRD_sum), clinical$Subtype,
                     p.adjust.method = "BH")


my_anova <- aov(HRD_sum ~ Subtype, data = sub_df_clinical)

# Print ANOVA summary
summary(my_anova)


# Density plot per subtype

hrd_dist_subtypes = ggplot(sub_df_clinical, aes(x = as.numeric(HRD_sum), fill = Subtype)) +
  geom_density(alpha = 0.5) +
  facet_wrap(~Subtype, ncol = 1) +
  labs(x = "Score", y = "Density", title = "Score Distribution by Subtype")

ggsave('../data/figures/BRCA_cohort/hrd_dist_subtypes.png',hrd_dist_subtypes, width = 8.63, height = 5.71)


ggplot(sub_df_clinical, aes(x = as.numeric(HRD_sum), fill = Subtype)) +
  geom_density(alpha = 0.4) +
  labs(x = "Score", y = "Density", title = "Score Distribution by Subtype") +
  scale_fill_discrete(name = "Subtype")




#### Create Statisic dataframe for the subtypes ####

subtypes = unique(clinical$Subtype)

col_names = c('Subtype', 'n', 'HRDs_min', 'HRDs_q1', 'HRDs_med','HRDs_mean', 'HRDs_q3', 'HRDs_max',
                                 'LOH_min', 'LOH_q1', 'LOH_med','LOH_mean', 'LOH_q3', 'LOH_max', 
                                 'LST_min', 'LST_q1', 'LST_med','LST_mean', 'LST_q3', 'LST_max',
                                 'TAI_min', 'TAI_q1', 'TAI_med','TAI_mean', 'TAI_q3', 'TAI_max')

subtype_statistics <- data.frame(matrix(ncol = length(col_names), nrow = 0))
colnames(subtype_statistics) <- col_names


for (subtype in subtypes){
  sub_df = clinical[clinical$Subtype == subtype,]
  stat_HRDs = as.array(summary(sub_df$HRD_sum))
  stat_LOH = as.array(summary(sub_df$LOH))
  stat_LST = as.array(summary(sub_df$LST))
  stat_TAI = as.array(summary(sub_df$TAI))
  number_cases = nrow(sub_df)
  
  
  subtype_stats = data.frame(
    Subtype = subtype,
    n = number_cases,
    HRDs_min = stat_HRDs[1],
    HRDs_q1 = stat_HRDs[2],
    HRDs_med = stat_HRDs[3],
    HRDs_mean = stat_HRDs[4],
    HRDs_q3 = stat_HRDs[5],
    HRDs_max = stat_HRDs[6],
    LOH_min = stat_LOH[1],
    LOH_q1 = stat_LOH[2],
    LOH_med = stat_LOH[3],
    LOH_mean = stat_LOH[4],
    LOH_q3 = stat_LOH[5],
    LOH_max = stat_LOH[6],
    LST_min = stat_LST[1],
    LST_q1 = stat_LST[2],
    LST_med = stat_LST[3],
    LST_mean = stat_LST[4],
    LST_q3 = stat_LST[5],
    LST_max = stat_LST[6],
    TAI_min = stat_TAI[1],
    TAI_q1 = stat_TAI[2],
    TAI_med = stat_TAI[3],
    TAI_mean = stat_TAI[4],
    TAI_q3 = stat_TAI[5],
    TAI_max = stat_TAI[6]
  )
  
  subtype_statistics = rbind(subtype_statistics, subtype_stats)
  #sub_df <- reshape2::melt(sub_df)
  #name = paste('sub_df_',subtype,sep = '')
  #assign(name,sub_df)
}
write.csv(subtype_statistics, '../data/statistics/BRCA_subtype_statistics.csv', row.names = FALSE)



#### Tried out to plot statisitc (not done yet) ####
ggplot(sub_df_BRCA_Basal) +
  aes(x = value, fill = variable) +
  geom_density(adjust = 1L, alpha = 0.5) +
  scale_fill_hue(direction = 1) +
  theme_minimal() +
  facet_wrap(vars(variable))


#### Subselect only hrd genes

raw_count_hrdgenes = raw_counts[rownames(raw_counts) %in% hrdgenes$gene,]

write.csv(raw_count_hrdgenes,'../data/raw_counts_hrdgenes.csv', row.names = TRUE)


dds_sub <- DESeqDataSetFromMatrix(countData = raw_count_hrdgenes,
                              colData = clinical,
                              design= ~ Subtype)

keep <- rowSums(counts(dds_sub)) >= 10
dds_sub <- dds_sub[keep,]

dds_sub <- DESeq(dds_sub)
resultsNames(dds_sub)




vsd_sub <- varianceStabilizingTransformation(dds_sub, blind=FALSE)

#ress = assay(vsd_sub)
#write.csv(ress,'../data/de_hrdgenes_vsd.csv', row.names = TRUE)

plotPCA(vsd_sub, intgroup=c("Subtype"))

pcaData <- plotPCA(vsd_sub, intgroup=c("Subtype", 'HRD_sum'), returnData=TRUE)
percentVar <- round(100 * attr(pcaData, "percentVar"))

pca_sub_hrd = ggplot(pcaData, aes(PC1, PC2, color= as.numeric(HRD_sum), shape=Subtype)) +
  geom_point(size=3) +
  xlab(paste0("PC1: ",percentVar[1],"% variance")) +
  ylab(paste0("PC2: ",percentVar[2],"% variance")) +
  scale_color_gradient(low = alpha("white", 0.7), high = "darkred")+
  coord_fixed()

pca_hrd = ggplot(pcaData, aes(PC1, PC2, color= as.numeric(HRD_sum))) +
  geom_point(size=3) +
  xlab(paste0("PC1: ",percentVar[1],"% variance")) +
  ylab(paste0("PC2: ",percentVar[2],"% variance")) +
  scale_color_gradient(low = alpha("gray", 0.7), high = "darkred")+
  coord_fixed()

pca_subtype = ggplot(pcaData, aes(PC1, PC2, color= Subtype)) +
  geom_point(size=3) +
  xlab(paste0("PC1: ",percentVar[1],"% variance")) +
  ylab(paste0("PC2: ",percentVar[2],"% variance")) +
  coord_fixed()

print(pca_hrd)

res_3 = assay(vsd_sub)

##### DBSCAN test #####
res_3 = assay(vsd_sub)

db = dbscan(t(res_3), eps = 0.15)

fviz_cluster(db, t(res_3), geom = 'point')

#### Hierachical clustering ####

dist_mat = dist(t(res_3))
hclust_avg <- hclust(dist_mat, method = 'average')
plot(hclust_avg)

#### SVD ####
svd_result = svd(t(res_3))

k = 2

S_k = diag(svd_result$d[1:k])
U_k = svd_result$u[, 1:k]
V_k = svd_result$v[, 1:k]

A_k = U_k %*% S_k %*% t(V_k)

df <- data.frame(x = A_k[, 1], y = A_k[, 2])
ggplot(df, aes(x = x, y = y)) + geom_point()

### SAVE results
write.csv(res_3,'../data/res_hrdgenes.csv', row.names = TRUE)



# Compare Basal to each other type individually
bahe = results(dds, contrast = c('Subtype', 'BRCA_Basal', 'BRCA_Her2'), alpha = 0.05)
bano = results(dds, contrast = c('Subtype', 'BRCA_Basal', 'BRCA_Normal'), alpha = 0.05)
bala = results(dds, contrast = c('Subtype', 'BRCA_Basal', 'BRCA_LumA'), alpha = 0.05)
balb = results(dds, contrast = c('Subtype', 'BRCA_Basal', 'BRCA_LumB'), alpha = 0.05)
baun = results(dds, contrast = c('Subtype', 'BRCA_Basal', 'undefined'), alpha = 0.05)


volcano_plotter(bahe, hrd_genes = hrd_genes, 'Her2', 'Basal', '../data/figures/BRCA_cohort/')
volcano_plotter(bano, hrd_genes = hrd_genes, 'Normal', 'Basal', '../data/figures/BRCA_cohort/')
volcano_plotter(bala, hrd_genes = hrd_genes, 'LumA', 'Basal', '../data/figures/BRCA_cohort/')
volcano_plotter(balb, hrd_genes = hrd_genes, 'LumB', 'Basal', '../data/figures/BRCA_cohort/')
volcano_plotter(baun, hrd_genes = hrd_genes, 'Undefined', 'Basal', '../data/figures/BRCA_cohort/')


res_gsea_bahe = kegg_hallmark_hrdgenes_gsea(bahe, hrdgenes, kegg_symbols, hallmark, 'Her2', 'Basal', '../data/figures/BRCA_cohort/')
res_gsea_bala = kegg_hallmark_hrdgenes_gsea(bala, hrdgenes, kegg_symbols, hallmark, 'LumA', 'Basal', '../data/figures/BRCA_cohort/')
res_gsea_balb = kegg_hallmark_hrdgenes_gsea(balb, hrdgenes, kegg_symbols, hallmark, 'LumB', 'Basal', '../data/figures/BRCA_cohort/')
res_gsea_bano = kegg_hallmark_hrdgenes_gsea(bano, hrdgenes, kegg_symbols, hallmark, 'Normal', 'Basal', '../data/figures/BRCA_cohort/')
res_gsea_baun = kegg_hallmark_hrdgenes_gsea(baun, hrdgenes, kegg_symbols, hallmark, 'Undefined', 'Basal', '../data/figures/BRCA_cohort/')


#### Quantile Analysis all subtypes ####

clinical_q1 = clinical[which(clinical$HRD_sum <= quantile(clinical$HRD_sum, 0.25)), ]
clinical_q3 = clinical[which(clinical$HRD_sum >= quantile(clinical$HRD_sum, 0.75)), ]

clinical_q1$Quantile = 'Q1'
clinical_q3$Quantile = 'Q3'

clinical_quantile = rbind(clinical_q1,clinical_q3)
cases_quantile = rownames(clinical_quantile)

raw_counts_quantile = raw_counts[,colnames(raw_counts) %in% cases_quantile]

cols = colnames(raw_counts_quantile)
sorted_rownames = intersect(cols,cases_quantile)
clinical_quantile = clinical_quantile[sorted_rownames,]


dds_quantile <- DESeqDataSetFromMatrix(countData = raw_counts_quantile, colData = clinical_quantile, design= ~ Quantile)


keep <- rowSums(counts(dds_quantile)) >= 10
dds_quantile <- dds_quantile[keep,]

dds_quantile <- DESeq(dds_quantile)
resultsNames(dds_quantile)

vsd <- varianceStabilizingTransformation(dds_quantile, blind=FALSE)

plotPCA(vsd, intgroup=c("Quantile"))

pcaData <- plotPCA(vsd, intgroup=c("Quantile", 'HRD_sum'), returnData=TRUE)
percentVar <- round(100 * attr(pcaData, "percentVar"))

pca_sub_hrd = ggplot(pcaData, aes(PC1, PC2, color= as.numeric(HRD_sum), shape=Quantile)) +
  geom_point(size=3) +
  xlab(paste0("PC1: ",percentVar[1],"% variance")) +
  ylab(paste0("PC2: ",percentVar[2],"% variance")) +
  scale_color_gradient(low = alpha("white", 0.7), high = "darkred")+
  coord_fixed()

pca_hrd = ggplot(pcaData, aes(PC1, PC2, color= as.numeric(HRD_sum))) +
  geom_point(size=3) +
  xlab(paste0("PC1: ",percentVar[1],"% variance")) +
  ylab(paste0("PC2: ",percentVar[2],"% variance")) +
  scale_color_gradient(low = alpha("gray", 0.7), high = "darkred")+
  coord_fixed()

pca_subtype = ggplot(pcaData, aes(PC1, PC2, color= Quantile)) +
  geom_point(size=3) +
  xlab(paste0("PC1: ",percentVar[1],"% variance")) +
  ylab(paste0("PC2: ",percentVar[2],"% variance")) +
  coord_fixed()

q1q3 = results(dds_quantile, contrast = c('Quantile', 'Q3', 'Q1'), alpha = 0.05)

volcano_plotter(q1q3, hrd_genes = hrd_genes, 'Q1', 'Q3', '../data/figures/BRCA_cohort/')
res_gsea_q1q3 = kegg_hallmark_hrdgenes_gsea(q1q3, hrdgenes, kegg_symbols, hallmark, 'Q1', 'Q3', '../data/figures/BRCA_cohort/')

res_q1q3 = assay(vsd)
write.csv(res_q1q3,'../data/res_q1q3.csv', row.names = TRUE)

pca_hrd


q1q3_df = as.data.frame(q1q3)
q1q3_df = q1q3_df[rownames(q1q3_df) %in% hrdgenes$gene,]
q1q3_df = q1q3_df[order(q1q3_df$padj),]

genes_sig = rownames(q1q3_df)[0:40]

vst <- vst(dds_quantile, blind=FALSE)
vst <- assay(vst)
vst <- as.data.frame(vst)
vst_sig <- vst[rownames(vst) %in% hrdgenes$gene,]
tvst <- t(vst_sig)
rownames(tvst) <- NULL

heat <- t(scale(tvst))

colnames(heat) = paste0(dds_quantile@colData$Quantile)

type = dds_quantile@colData$Quantile
col_order = colnames(heat)


heat = heat[rownames(heat) %in% genes_sig,]

ha = HeatmapAnnotation(Subtype = colnames(heat), col = list(Subtype = c('Q3' = 'darkblue', 'Q1' = 'lightblue')))

hrds_values = dds_quantile@colData$HRD_sum
loh_values = dds_quantile@colData$LOH
lst_values = dds_quantile@colData$LST
tai_values = dds_quantile@colData$TAI

hrd_annotation = HeatmapAnnotation(HRDs = hrds_values, col = list(HRDs = colorRamp2(c(min(hrds_values),max(hrds_values)),c('beige','saddlebrown'))))
loh_annotation = HeatmapAnnotation(LOH = loh_values, col = list(LOH = colorRamp2(c(min(loh_values),max(loh_values)),c('lightgreen', 'darkgreen'))))
lst_annotation = HeatmapAnnotation(LST = lst_values, col = list(LST = colorRamp2(c(min(lst_values),max(lst_values)),c('lightblue','darkblue'))))
tai_annotation = HeatmapAnnotation(TAI = tai_values, col = list(TAI = colorRamp2(c(min(tai_values),max(tai_values)),c('lightcyan', 'darkcyan'))))

#png(file='../data/figures/BRCA_cohort/heatmap.png')
hm = Heatmap(heat,show_column_names = FALSE, row_names_gp = gpar(fontsize = 4), bottom_annotation = ha, column_km = 4, name = 'z-score') %v% 
  hrd_annotation %v% loh_annotation %v% lst_annotation %v% tai_annotation

draw(hm)

# Gene list for mutation data
genes = paste(hrdgenes$gene, collapse = ',')



#### Trying to construct score ####


dds <- DESeqDataSetFromMatrix(countData = raw_counts,
                              colData = clinical,
                              design= ~ Subtype)

keep <- rowSums(counts(dds)) >= 10
dds <- dds[keep,]

dds <- DESeq(dds)
resultsNames(dds)

res = results(dds, alpha = 0.01)

res_df = as.data.frame(res)

res_df_hrd = subset(res_df, rownames(res_df) %in% hrdgenes$gene)

norm_counts = counts(dds, normalized = TRUE)

norm_counts_hrd = subset(norm_counts, rownames(norm_counts) %in% hrdgenes$gene)

# Using quantile analysis (Run quantile dds first)

res_quant = results(dds_quantile, alpha = 0.01)

res_quant_df = as.data.frame(res_quant)

res_quant_df_hrd = subset(res_quant_df, rownames(res_quant_df) %in% hrdgenes$gene)

norm_counts_quant = counts(dds_quantile, normalized = TRUE)

norm_counts_quant_hrd = subset(norm_counts_quant, rownames(norm_counts_quant) %in% hrdgenes$gene)



hrd_genes_scores = as.data.frame(hrdgenes)

hrd_genes_scores$score = ''

for (gene in rownames(res_quant_df_hrd)){
  score = abs(res_quant_df_hrd[gene,'padj']) * abs(res_quant_df_hrd[gene,'log2FoldChange'])
  #print(paste(gene,score,sep = ' '))
  hrd_genes_scores[hrd_genes_scores['gene']== gene,'score'] = score
}

hrd_genes_scores = subset(hrd_genes_scores, gene != 'GIYD1')


clinical$hhrscore = apply(norm_counts_hrd, 2, function(x) sum(as.numeric(x) * as.numeric(hrd_genes_scores$score)))

clinical_quantile$hhrscore = apply(norm_counts_quant_hrd, 2, function(x) sum(as.numeric(x) * as.numeric(hrd_genes_scores$score)))

clinical$sumgene = apply(norm_counts_hrd, 2, function(x) sum(as.numeric(x)))

ggplot(clinical_quantile, aes(x = HRD_sum, y = hhrscore)) + 
  geom_point() +
  labs(x = "HRD_sum", y = "HRR_score") +
  theme_bw()

ggplot(clinical, aes(x = HRD_sum, y = hhrscore)) + 
  geom_point() +
  labs(x = "HRD_sum", y = "HRR_score") +
  theme_bw()

ggplot(clinical, aes(x = HRD_sum, y = sumgene)) + 
  geom_point() +
  labs(x = "HRD_sum", y = "SumGeneCount") +
  theme_bw()


#### PCA with HRD scores and gene counts
norm_counts_hrd_t = t(norm_counts_hrd)
scores = clinical[,c('HRD_sum','LOH','LST','TAI')] 
scores_ncount = cbind(scores,norm_counts_hrd_t)


pca = prcomp(scores_ncount, scale. = TRUE)

plot(pca, type = "l", main = "Scree Plot")
biplot(pca, scale = 0)

plot(pca$x[,1], pca$x[,2], main = "Scores Plot")

pca_df <- data.frame(Case = rownames(pca$x), PC1 = pca$x[,1], PC2 = pca$x[,2], HRD_sum = clinical$HRD_sum)

ggplot(pca_df, aes(x = PC1, y = PC2, color = HRD_sum)) +
  geom_point() +
  labs(title = "PCA normalized counts HRD genes and HRD score", x = "PC1", y = "PC2")

hclust.out <- hclust(dist(norm_counts_hrd_t))
plot(hclust.out)
abline(h = 7, col = "red")


#### Factor Analysis ####
count.fa = factanal(norm_counts_hrd_t, factors = 2)
count.fa

Lambda <- count.fa$loadings
Psi <- diag(count.fa$uniquenesses)
S <- count.fa$correlation
Sigma <- Lambda %*% t(Lambda) + Psi
round(S - Sigma, 6)


plot(count.fa$loadings[,1], 
     count.fa$loadings[,2],
     xlab = "Factor 1", 
     ylab = "Factor 2", 
     ylim = c(-1,1),
     xlim = c(-1,1),
     main = "Varimax rotation")
text(count.fa$loadings[,1]-0.04, 
     count.fa$loadings[,2]+0.04,
     colnames(norm_counts_hrd_t),
     col="blue",
     cex=0.55)
abline(h = 0, v = 0)

#### Corr Matrix Counts ####

cov_matrix = cor(norm_counts_hrd_t)

my_palette <- colorRampPalette(c("#FDFD96", "#FF0000"))

corrplot(cov_matrix, method="number")

heatmap.2(cov_matrix, 
          main = "Covariance Matrix",
          xlab = "Variables", ylab = "Variables",
          key = TRUE, keysize = 1.2, key.title = NA,
          trace = "none", dendrogram = "none",
          col = my_palette(256))

#### T-sne ####

tsne_out <- Rtsne(norm_counts_hrd_t)

tsne_plot <- data.frame(x = tsne_out$Y[,1],
                        y = tsne_out$Y[,2])

ggplot(tsne_plot)+ geom_point(aes(x=x,y=y))


#### Z-Score as weight ####

vst <- vst(dds_quantile, blind=FALSE)
vst <- assay(vst)
vst <- as.data.frame(vst)
vst_sig <- vst[rownames(vst) %in% hrdgenes$gene,]

#tvst <- t(vst_sig)
#rownames(tvst) <- NULL

z_score <- scale(vst_sig)

#heat = heat[rownames(heat) %in% genes_sig,]


res_quant = results(dds_quantile, alpha = 0.01)
res_quant_df = as.data.frame(res_quant)
res_quant_df_hrd = subset(res_quant_df, rownames(res_quant_df) %in% rownames(z_score))

norm_counts_quant = counts(dds_quantile, normalized = TRUE)
norm_counts_quant_hrd = subset(norm_counts_quant, rownames(norm_counts_quant) %in% rownames(z_score))


hrr_score_df = norm_counts_quant_hrd * z_score

clinical_quantile$hrr_score = apply(hrr_score_df, 2, function(x) sum(as.numeric(x)))

ggplot(clinical_quantile, aes(x = HRD_sum, y = hrr_score)) + 
  geom_point() +
  labs(x = "HRD_sum", y = "HRR_score") +
  theme_bw()

pca = prcomp(t(hrr_score_df), scale. = TRUE)

plot(pca, type = "l", main = "Scree Plot")

plot(pca$x[,1], pca$x[,2], main = "Scores Plot")

pca_df <- data.frame(Case = rownames(pca$x), PC1 = pca$x[,1], PC2 = pca$x[,2], HRD_sum = clinical_quantile$HRD_sum)

ggplot(pca_df, aes(x = PC1, y = PC2, color = HRD_sum)) +
  geom_point() +
  labs(title = "PCA normalized counts HRD genes and HRD score", x = "PC1", y = "PC2")

#### For YAN ####

clinical$HRD = ifelse(clinical$HRD_sum >= 42, 'High', 'Low')

dds <- DESeqDataSetFromMatrix(countData = raw_counts,
                              colData = clinical,
                              design= ~ HRD)


keep <- rowSums(counts(dds)) >= 10
dds <- dds[keep,]

dds$HRD = relevel(dds$HRD, ref = 'Low')

dds <- DESeq(dds)
res = results(dds, alpha = 0.01)

gene_list = res$stat
res$symbol = rownames(res)
names(gene_list) = make.names(res$symbol, unique = T)
gene_list = gene_list[order(gene_list, decreasing = T)]

gene_names = make.names(res$symbol, unique = T)

kegg_gsea = GSEA(gene_list, TERM2GENE = kegg_symbols,
                 pvalueCutoff = 0.1,
                 eps=1e-50,
                 seed=T)

kegg_gsea_df = as.data.frame(kegg_gsea)

hallmark_gsea = GSEA(gene_list, TERM2GENE = hallmark,
                     pvalueCutoff = 0.1,
                     eps=1e-50,
                     seed=T)

hallmark_gsea_df = as.data.frame(hallmark_gsea)

write.csv(hallmark_gsea_df,'../data/Hallmark_Yan_v2.csv', row.names = TRUE)


#### For Yan OV ####

raw_counts_ov = read.csv('../data/gene_expression/prepared_data/raw_count_OV.csv', sep = ',')
clinical_ov = read.csv('../data/gene_expression/prepared_data/clinical_OV.csv', sep = ',')

genes = raw_counts_ov[,1]

counters = rep(1, length(genes))

duplicates <- genes[duplicated(genes)]
if (length(duplicates) > 0) {
  print("Duplicates found:")
  print(duplicates)
  for (i in 1:length(genes)){
    if (genes[i] %in% duplicates){
      genes[i] = paste(genes[i],i,sep = '_')
      counters[i] = counters[i] + 1
    }
  }
} else {
  print("No duplicates found.")
}

duplicates <- genes[duplicated(genes)]
if (length(duplicates) > 0) {
  print("Duplicates found:")
  print(duplicates)
} else {
  print("No duplicates found.")
}

raw_counts_ov = raw_counts_ov[,-1]
rownames(raw_counts_ov) = genes

col_names = colnames(raw_counts_ov)
rownames(clinical_ov) = col_names
clinical_ov = clinical_ov[,-1]

clinical_ov$HRD = ifelse(clinical_ov$HRD_sum >= 42, 'High', 'Low')


high = clinical_ov[clinical_ov$HRD == 'High',]
low = clinical_ov[clinical_ov$HRD == 'Low',]

dds <- DESeqDataSetFromMatrix(countData = raw_counts_ov,
                              colData = clinical_ov,
                              design= ~ HRD)


keep <- rowSums(counts(dds)) >= 10
dds <- dds[keep,]

dds$HRD = relevel(dds$HRD, ref = 'Low')

dds <- DESeq(dds)
res = results(dds, alpha = 0.01)

gene_list = res$stat
res$symbol = rownames(res)
names(gene_list) = make.names(res$symbol, unique = T)
gene_list = gene_list[order(gene_list, decreasing = T)]

gene_names = make.names(res$symbol, unique = T)

kegg_gsea = GSEA(gene_list, TERM2GENE = kegg_symbols,
                 pvalueCutoff = 0.1,
                 eps=1e-50,
                 seed=T)

kegg_gsea_df = as.data.frame(kegg_gsea)

hallmark_gsea = GSEA(gene_list, TERM2GENE = hallmark,
                     pvalueCutoff = 0.1,
                     eps=1e-50,
                     seed=T)

hallmark_gsea_df = as.data.frame(hallmark_gsea)

write.csv(hallmark_gsea_df,'../data/Hallmark_Yan_OV.csv', row.names = TRUE)

