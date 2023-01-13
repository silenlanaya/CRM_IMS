library(SingleCellSignalR)
library(Seurat)
library(tidyverse)
library(glmnet)
library(survminer)
library(data.table)
library(reshape2)
library(plyr)
library(dplyr)
library(edgeR)
library(org.Hs.eg.db)
library(clusterProfiler)
library(GSVA)
library(GSEABase)
library(survival)
library(ConsensusClusterPlus)
library(dendsort)
library(survminer)
library(CMScaller)
library(sva)
library(ArchR)
library(ggExtra)
library(ggtern)
library(ggplot2)
library(Seurat)
library(ggrepel)
library(ggbreak)
library(ggbiplot)
library(ggpubr)
library(forestplot)
library(pheatmap)
library(ggthemes)
library(ggforce)
library(circlize)
options(scipen = 50)
Sys.setenv(LANGUAGE = "en")
options(stringsAsFactors = FALSE) 
paletteLength <- 100
bk = unique(c(seq(-1, 1, length = 100)))
color_subtype <- c('#6699FF', '#FFCC00', '#FF6633', '#99CC66')

##GEO training cohort
test_infor <- read.csv('Test_GEO_infor.csv', row.names = 1)
test_expr <- read.csv('Test_GEO_expr.csv', row.names = 1)

mpi <- read.csv('MPI_network.csv')
mpi_gene <- unique(mpi$b)
mpi_gene %<>% str_split(' ') %>% map(., function(x){x[[1]]}) %>% unlist

meta_gene <- test_expr[mpi_gene, ]  %>%.[complete.cases(.), ] %>% t %>% data.frame
surv_data <- merge(test_infor, meta_gene, by.x = 1, by.y = 0)
surv_data$rfs.event <- as.numeric(surv_data$rfs.event)
res.cut <- surv_cutpoint(surv_data, time = "rfs", 
                         event = "rfs.event", 
                         variables = names(surv_data)[4:ncol(surv_data)], 
                         minprop = 0.25)
res.cat <- surv_categorize(res.cut)
for (i in names(res.cat)[3:ncol(res.cat)]){
  res.cat[[i]] <- factor(res.cat[[i]], levels = c('low', 'high'))
}
res.cat[1:3, 1:4]
res <- data.frame()
my.surv <- Surv(res.cat$rfs, res.cat$rfs.event)
uni <- function(x){
  group <- res.cat[[x]]
  survival_dat <- data.frame(group=group)
  fit <- survfit(my.surv ~ group)
  m=coxph(my.surv ~ group, data = survival_dat)
  tem <- summary(m)
  coef <- tem$coefficients[[1]]
  aoe <- cbind(as.data.frame(summary(m)[["coefficients"]]),as.data.frame(summary(m)[["conf.int"]]))
  abc <- t(aoe)[,1]
  HR <- round(abc[6],2)
  CI <- paste(paste(round(abc[8],2), round(abc[9],2), sep = " - "), sep = "")
  res_now <- data.frame(symbol = names(res.cat[x]), pvalue = abc[5], HR = HR, CI = CI, coef = coef)
  res_now
}
dim(res.cat)
tem <- colnames(res.cat)[3:(ncol(res.cat))]
res <- map(tem, uni)
tem_res <- res

res <- c()
for (i in c(1:length(tem_res))){
  res <- rbind(res, tem_res[[i]])
}
head(res)
res <- res[order(res$pvalue), ]
res_sig <- subset(res, res$pvalue < 0.01)
metabolism_sig <- meta_protein[,res_sig$symbol] 
result_meta_threshold <- mutate(res, threshold=ifelse(coef > 0 & pvalue < 0.05,"A", ifelse(coef < 0 & pvalue <0.05, "B", "C")))
result_meta_threshold$logp <- -log10(result_meta_threshold$pvalue)


#immune-related genes with prognosis significance
immune_gene <- read.csv('InnateDB_genes.csv')
immune_gene <- immune_gene$name
immu_gene <- test_expr[immune_gene, ]  %>%.[complete.cases(.), ] %>% t %>% data.frame
surv_data <- merge(test_infor, immu_gene, by.x = 1, by.y = 0)
surv_data$rfs.event <- as.numeric(surv_data$rfs.event)
res.cut <- surv_cutpoint(surv_data, time = "rfs", 
                         event = "rfs.event", 
                         variables = names(surv_data)[4:ncol(surv_data)], 
                         minprop = 0.25)
res.cat <- surv_categorize(res.cut)
for (i in names(res.cat)[3:ncol(res.cat)]){
  res.cat[[i]] <- factor(res.cat[[i]], levels = c('low', 'high'))
}
res.cat[1:3, 1:4]
res_immu <- data.frame()
my.surv <- Surv(res.cat$rfs, res.cat$rfs.event)

tem <- colnames(res.cat)[3:(ncol(res.cat))]
res_immu <- map(tem, uni)
tem_res <- res_immu

res_immu <- c()
for (i in c(1:length(tem_res))){
  res_immu <- rbind(res_immu, tem_res[[i]])
}
nrow(res_immu)
res_immu <- res_immu[order(res_immu$pvalue), ]
res_immu_sig <- subset(res_immu, res_immu$pvalue < 0.01)

###CMS subtyping
cms_geo <- CMScaller(test_expr,
                     rowNames = "symbol" ,   ##"entrez"; "symbol"; "ensg" 
                     RNAseq=F,
                     doPlot=TRUE)
sum(is.na(cms_geo$prediction))
cms_geo <- cms_geo[complete.cases(cms_geo$prediction), ]
cms_geo <- merge(cms_geo, test_infor, by.x = 0, by.y = 1)
cms_geo$rfs.event <- as.numeric(cms_geo$rfs.event)
names(cms_geo)[names(cms_geo) == 'prediction'] <- 'subtype_CMS'
cms_geo$subtype_CMS <- str_sub(cms_geo$subtype_CMS, 4)
cms_geo$subtype_CMS <- as.numeric(cms_geo$subtype_CMS)

###IMS subtyping
metabolism_pathway <- read.csv('metabolism_pw.csv')
stromal_pathway <- read.csv('stromal_oncogenic.csv')
names(stromal_pathway) <- c('pathway', 'symbol')
tme_pathway <- read.csv('TME_related.csv')
names(tme_pathway) <- c('pathway', 'symbol')
all_pathway <- rbind(metabolism_pathway, c7) %>%
  rbind(., stromal_pathway) %>%
  rbind(., tme_pathway)
metabolism_pathway <- read.csv('metabolism_pw.csv')
metabolism_pathway <- metabolism_pathway %>% 
  split(., .$pathway) %>% 
  lapply(., function(x)(x$symbol))
metabolism_pathway <- lapply(metabolism_pathway, unique)
metabolism_pathway <- lapply(metabolism_pathway, as.character)
metabolism_GSEA_geo <- as.data.frame(t(gsva(as.matrix(combat_edata), metabolism_pathway, method = "ssgsea")))
consensus_meta_geo <- data.frame(t(metabolism_GSEA_geo))
c7 <- read.csv('immune_c7.csv')
c7 <- c7 %>% 
  split(., .$signature ) %>% 
  lapply(., function(x)(x$gene))
c7 <- lapply(c7, unique)
c7 <- lapply(c7, as.character)
c7_score <- as.data.frame(t(gsva(as.matrix(test_expr), c7, method = "ssgsea")))
stromal_pathway <- read.csv('stromal_oncogenic.csv')
stromal_pathway <- stromal_pathway %>% 
  split(., .$pathway) %>% 
  lapply(., function(x)(x$gene))
stromal_pathway <- lapply(stromal_pathway, unique)
stromal_pathway <- lapply(stromal_pathway, as.character)
stromal_score <- as.data.frame(t(gsva(as.matrix(test_expr), stromal_pathway, method = "ssgsea")))
tme_pathway <- read.csv('TME_related.csv')
tme_pathway <- tme_pathway %>% 
  split(., .$pathway) %>% 
  lapply(., function(x)(x$gene))
tme_pathway <- lapply(tme_pathway, unique)
tme_pathway <- lapply(tme_pathway, as.character)
tme_score <- as.data.frame(t(gsva(as.matrix(test_expr), tme_pathway, method = "ssgsea")))
meta_expr_geo <- rbind(consensus_meta_geo, data.frame(t(c7_score)))
meta_expr_geo <- rbind(meta_expr_geo, data.frame(t(tme_score)))
meta_expr_geo <- rbind(meta_expr_geo, data.frame(t(stromal_score)))
meta_expr_geo <- rbind(meta_expr_geo, all_keyexpr)
meta_expr_geo <- meta_expr_geo[complete.cases(meta_expr_geo), ]
meta_expr_geo <- data.frame(scale(t(meta_expr_geo))) %>% t()%>% data.frame()
results_meta = ConsensusClusterPlus(as.matrix(meta_expr_geo), maxK = 6, reps = 50, pItem = 0.8, pFeature = 1,
                                    title = title, clusterAlg = "hc", distance = "pearson", 
                                    seed = 1262118388.71279, plot = "png")
cluster_meta <- results_meta[[3]][["consensusClass"]]
cluster_meta <- data.frame(sampleID = names(cluster_meta), cluster_meta = cluster_meta)
cluster_meta %>%
  group_by(cluster_meta) %>%
  dplyr::summarise(n = n())
cluster_meta <- merge(cluster_meta, test_infor, by.x = 1, by.y = 1)

##immune-related cell abundance
bulk_immunity <- as.data.frame(t(gsva(as.matrix(test_expr), immunity, method = "ssgsea")))
bulk_immunity <- data.frame(t(bulk_immunity))


tri_immunity <- data.frame(t(bulk_immunity))
tri_immunity <- merge(tri_immunity, cluster_col, by.x = 0, by.y = 0)
rownames(tri_immunity) <- tri_immunity$Row.names
tri_immunity <- tri_immunity[,-1]
tri_immunity <- aggregate(tri_immunity[,c(1:(ncol(tri_immunity)-1))], by = list(tri_immunity$Cluster), FUN = mean)
rownames(tri_immunity) <- tri_immunity$Group.1
tri_immunity <- tri_immunity[,-1]
tri_immunity <- data.frame(t(tri_immunity))
result1 <- tri_immunity

for(i in 1:ncol(tri_immunity)){ 
  
  df<-((tri_immunity[,i]-min(tri_immunity[,i]))/(max(tri_immunity[,i])-min(tri_immunity[,i])))
  
  result1[,i]<-as.data.frame(df)
  
}

tri_immunity <- result1


tri_immunity$total <- rowSums(tri_immunity)
ggtern(data=tri_immunity, aes(x=tri_immunity$C1,
                              y=tri_immunity$C2,
                              z=tri_immunity$C3)) +  
  geom_point(size = tri_immunity$total)  

hist(tri_immunity$total)

tri_immunity$totalfrequency <- ifelse(tri_immunity$total >= 1.5,
                                      1.5^(tri_immunity$total) * 3, 
                                      3 * (tri_immunity$total))

ggtern(data=tri_immunity, aes(x=C1,
                              y=C2,
                              z=C3)) +  
  geom_point(size = tri_immunity$totalfrequency)  

tri_immunity$cell <- rownames(tri_immunity)

p <- ggtern(data=tri_immunity, aes(x=C1,
                                   y=C2,
                                   z=C3)) +  
  geom_point(size = tri_immunity$totalfrequency) +
  theme_custom(20, '') +
  theme(tern.panel.background = element_rect(fill = "white"),
        tern.panel.grid.minor = element_line(color = "gray90"), 
        tern.axis.arrow.show = TRUE,
        
        tern.axis.arrow.T = element_line(color ='#0000E3', size = 2.5), 
        tern.axis.arrow.L = element_line(color = '#FF44FF', size = 2.5),
        tern.axis.arrow.R = element_line(color = 'red', size = 2.5),
        
        tern.axis.arrow.text.L = element_text(color = 'black'),  
        tern.axis.arrow.text.T = element_text(color = 'black'),
        tern.axis.arrow.text.R = element_text(color = 'black'),
        
        tern.axis.arrow.sep = 0.1, 
        
        tern.panel.grid.major.T = element_line(color = 'gray92', linetype = 1, size = 0.8), 
        tern.panel.grid.major.L = element_line(color = 'gray92', linetype = 1, size = 0.8),
        tern.panel.grid.major.R = element_line(color = 'gray92', linetype = 1, size = 0.8),
        tern.panel.grid.minor.T = element_line(color = 'gray94', linetype = 1, size = 0.8), 
        tern.panel.grid.minor.L = element_line(color = 'gray94', linetype = 1, size = 0.8),
        tern.panel.grid.minor.R = element_line(color = 'gray94', linetype = 1, size = 0.8),
        
        tern.axis.title.L = element_text(color = '#FF44FF', size = 15),
        tern.axis.title.T = element_text(color = '#0000E3', size = 15),
        tern.axis.title.R = element_text(color = 'red', size = 15),
        
        tern.axis.text.L = element_text(size = 17,face = 'bold'),
        tern.axis.text.R = element_text(size = 17,face = 'bold'),
        tern.axis.text.T = element_text(size = 17,face = 'bold'),
        
        tern.axis.vshift = 0.04,
        
        tern.axis.line.T = element_line(size = 0.8),
        tern.axis.line.R = element_line(size = 0.8),
        tern.axis.line.L = element_line(size = 0.8)) + 
  
  geom_Lisoprop(color='darkgrey', value=.5, linetype=4, size=1) 
p
p1 <- p + 
  geom_point(aes(color=tri_immunity$C1), size=tri_immunity$totalfrequency, alpha=0.7) + #半透明
  geom_point(aes(color=tri_immunity$C2), size=tri_immunity$totalfrequency, alpha=0.7) +
  geom_point(aes(color=tri_immunity$C3), size=tri_immunity$totalfrequency, alpha=0.7) +
  scale_color_gradient2(low='red', mid = '#0000E3', high ='purple', midpoint = 0.33, #三种颜色
                        guide = FALSE) + 
  geom_point(size=tri_immunity$totalfrequency, shape = 1, alpha = 0.8,
             stroke = 0.7, 
             color = "black")
p1
p1 + geom_text(aes(label=tri_immunity$cell)) 




##MSI validation
msi_infor <- read.csv('MSI_infor.csv')
msi_expr <- read.table('MSI_protein_matrix.txt', head = TRUE)

msi_expr$Description %<>%
  str_split('=') %>%
  map(., function(x){x[4]}) %>%
  unlist

msi_expr$Description %<>%
  str_split(' ') %>%
  map(., function(x){x[1]}) %>%
  unlist

metabolism_score_msi <- as.data.frame(t(gsva(as.matrix(msi_expr), metabolism_pathway, method = "ssgsea")))
c7_score_msi <- as.data.frame(t(gsva(as.matrix(msi_expr), c7, method = "ssgsea")))
stromal_score_msi <- as.data.frame(t(gsva(as.matrix(msi_expr), 
                                          stromal_pathway, method = "ssgsea")))
tme_score_msi <- as.data.frame(t(gsva(as.matrix(msi_expr), tme_pathway, method = "ssgsea")))
all_keyexpr_msi <- msi_expr[all_keygene,]
all_keyexpr_msi <- all_keyexpr_msi[complete.cases(all_keyexpr_msi), ]
consensus_meta_msi <- data.frame(t(metabolism_score_msi))
meta_expr_msi <- list(meta_expr_msi, data.frame(t(c7_score_msi)),
                      data.frame(t(tme_score_msi)), data.frame(t(stromal_score_msi)),
                      all_keyexpr_msi)
meta_expr_msi <- do.call(meta_expr_msi, rbind)
names(meta_expr_msi) <- str_replace_all(names(meta_expr_msi), '[.]', '-')
meta_expr_msi <- meta_expr_msi[complete.cases(meta_expr_msi), ]
results_meta_msi = ConsensusClusterPlus(as.matrix(meta_expr_msi), maxK = 6, reps = 50, pItem = 0.8, pFeature = 1,
                                        title = title, clusterAlg = "hc", distance = "pearson", 
                                        seed = 1262118388.71279, plot = "png")

cluster_meta_msi <- results_meta_msi[[3]][["consensusClass"]]
cluster_meta_msi <- data.frame(sampleID = names(cluster_meta_msi), cluster_meta = cluster_meta_msi)
cluster_meta_msi %>%
  group_by(cluster_meta) %>%
  dplyr::summarise(n = n())

cluster_meta_msi <- merge(cluster_meta_msi, msi_infor, by.x = 'sampleID', by.y = 'sampleID')
rownames(cluster_meta_msi) <- cluster_meta_msi$sampleID
names(cluster_meta_msi)[names(cluster_meta_msi) == 'cluster_meta'] <- 'Cluster'
names(cluster_meta_msi)[names(cluster_meta_msi) == 'age'] <- 'Age'
names(cluster_meta_msi)[names(cluster_meta_msi) == 'gender'] <- 'Gender'
names(cluster_meta_msi)[names(cluster_meta_msi) == 'tnm'] <- 'TNM'
names(cluster_meta_msi)
anno_msi <- cluster_meta_msi[,c('Cluster', "Age", "Gender", "pT", "pN","M","TNM")]
anno_msi <- anno_msi[order(anno_msi$Cluster), ]

msi_immunity <- as.data.frame(t(gsva(as.matrix(msi_expr), immunity, method = "ssgsea")))
msi_immunity <- data.frame(t(msi_immunity))
msi_immunity <- msi_immunity[,match(rownames(anno_msi), names(msi_immunity))]
anno_msi <- anno_msi[,c(1, 2, 3, 4, 5, 7)]
anno_msi$Gender <- ifelse(anno_msi$Gender == 0, 'MALE', 'FEMALE')
anno_msi$TNM[anno_msi$TNM == 1] <- 'I'
anno_msi$TNM[anno_msi$TNM == 2] <- 'II'
anno_msi$TNM[anno_msi$TNM == 3] <- 'III'
anno_msi$TNM[anno_msi$TNM == 4] <- 'IV'
anno_msi <- anno_msi[,c('Cluster', 'Age', 'Gender')]
all(rownames(anno_msi) == names(msi_immunity))
anno_msi$Cluster <- str_c('C', anno_msi$Cluster)
anno_msi <- anno_msi[order(anno_msi$Cluster), ]
msi_immunity <- msi_immunity[,match(rownames(anno_msi), names(msi_immunity))]
anno_msi %>%
  group_by(Cluster) %>%
  dplyr::summarise(n = n())
pheatmap(msi_immunity, 
         annotation_colors =  anno_colors,
         gaps_col =  c(12, 20),
         breaks = unique(c(seq(-1, 1, length = 100))), 
         annotation_col = anno_msi, 
         colorRampPalette(ArchRPalettes$greenBlue)(paletteLength),
         cluster_col = F, 
         cluster_row = T,
         show_rownames = T,  
         show_colnames = T,
         fontsize = 31,
         scale = 'row',
         cellwidth = 30, 
         cellheight = 30,
         border_color = "white")

##single cell data
sc_data <- readRDS('sc156_seurat_obj.rds')
sc_data@meta.data %>%
  group_by(cellType) %>%
  summarize(n = n())
sc_color <- c( "#C7EAB2", "#5FC1C2", "#1B90BE",  '#FFCC00', '#FF6633','#CC9900', '#009900',
               "#6E568C","#E0367A","#D8D155","#64495D","#7CC767",
               "#D20A13","#FFD121","#088247","#11AA4D")
tem_num <- sc_data@meta.data %>%
  group_by(cellType) %>%
  dplyr::summarise(n = n())
tem_num <- tem_num[order(tem_num$n), ]
tem_num$cell <- 'Cell'
tem_num$cellType <- factor(tem_num$cellType, levels = tem_num$cellType)
names(tem_num)[names(tem_num) == 'cellType'] <- 'CellType'
ggplot(tem_num, aes(cell, n, fill = CellType, group = CellType))+ theme_classic() + 
  geom_bar(stat="identity",position="fill")+
  ggtitle("") + big +
  xlab('') + ylab('Fraction') + 
  theme( legend.title = element_text(size = 20),
         legend.text = element_text(size = 20),
         title=element_text(size=18))  +
  theme(axis.text.x = element_text( angle = 90)) +
  theme(strip.text.x = element_text(size = 17)) + 
  scale_y_continuous(labels = scales::percent) +
  scale_fill_manual(values = sc_color)

sc_data$cellType <- factor(sc_data$cellType, levels = tem_num$CellType)
DimPlot(sc_data, reduction = "tsne", 
        group.by = "cellType", label.size = 7,
        pt.size = 1,
        cols = sc_color)+ 
  ggtitle('Cell Cluster')+
  guides(colour = guide_legend(override.aes = list(size = 3)))  +  
  theme(legend.title = element_blank(),
        legend.text = element_text(size = 17), 
        legend.position = 'right',
        plot.title = element_text(size = 27, face = "bold"))+
  guides(colour = guide_legend(override.aes = list(size = 8)),
         size = guide_legend(order=3))+
  theme(axis.text.x = element_text(angle = 0, hjust = 1, size = 20),
        axis.text.y = element_text(angle = 0, hjust = 1, size = 20),
        axis.title.x=element_text(size = 20),
        axis.title.y=element_text(size = 20),)

detailed_p <- DimPlot(sc_data, reduction = "tsne", 
                      group.by = "cellType_anno", label.size = 7,
                      pt.size = 1,
                      cols = sc_color)+ 
  ggtitle('NR to Immunotherapy')+
  guides(colour = guide_legend(override.aes = list(size = 3)))  +  
  theme(legend.title = element_blank(),
        legend.text = element_text(size = 17), 
        legend.position = 'right',
        plot.title = element_text(size = 27, face = "bold"))+
  guides(colour = guide_legend(override.aes = list(size = 8)),
         size = guide_legend(order=3))+
  theme(axis.text.x = element_text(angle = 0, hjust = 1, size = 20),
        axis.text.y = element_text(angle = 0, hjust = 1, size = 20),
        axis.title.x=element_text(size = 20),
        axis.title.y=element_text(size = 20),)

Idents(sc_data) <- sc_data$cellType
all.markers <- FindAllMarkers(sc_data, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
top5 <- all.markers %>% group_by(cluster) %>% top_n(n = 5, wt = avg_log2FC)
tem <- DoHeatmap(object = sc_data, features = top5$gene, label = TRUE)+ 
  theme(text = element_text(size = 20))

Idents(sc_data) <- sc_data$cellType
all.markers <- FindAllMarkers(sc_data, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
top5 <- all.markers %>% group_by(cluster) %>% top_n(n = 5, wt = avg_log2FC)
DoHeatmap(object = sc_data, features = top5$gene, label = TRUE)+ 
  theme(text = element_text(size = 20))

FeaturePlot(sc_data, features = 'S100A9', pt.size = 1,  
            min.cutoff = "q10", max.cutoff = "q90", reduction = "tsne")+ 
  ggtitle('') +
  guides(colour = guide_legend(override.aes = list(size = 10)))+
  theme(legend.text = element_text(size = 25),
        axis.title.x=element_text(size = 25),
        axis.title.y=element_text(size = 25),
        axis.text.x=element_text(size = 25,color = "black"),
        axis.text.y=element_text(size = 25,color = "black"),
        plot.title = element_text(size = 25, face = "bold")) +  
  theme(legend.title = element_blank(),
        legend.text = element_text(size = 17), 
        legend.position = 'right',
        plot.title = element_text(size = 27, face = "bold"))+ 
  ggtitle('S100A9 Expression')


metabolism_score_sc <- as.data.frame(t(gsva(as.matrix(sc_data@assays$SCT@scale.data),
                                            metabolism_pathway, method = "ssgsea")))
c7_score_sc <- as.data.frame(t(gsva(as.matrix(sc_data@assays$SCT@scale.data), 
                                    c7, method = "ssgsea")))
stromal_score_sc <- as.data.frame(t(gsva(as.matrix(sc_data@assays$SCT@scale.data), 
                                         stromal_pathway, method = "ssgsea")))
tme_score_sc <- as.data.frame(t(gsva(as.matrix(sc_data@assays$SCT@scale.data), 
                                     tme_pathway, method = "ssgsea")))
tem_gene <- intersect(rownames(sc_data@assays$SCT@scale.data), all_keygene)
all_keyexpr_sc <- sc_data@assays$SCT@scale.data[tem_gene,]
all_keyexpr_sc <- all_keyexpr_sc[complete.cases(all_keyexpr_sc), ]

sc_list <- list(data.frame(t(metabolism_score_sc)), data.frame(t(c7_score_sc)),
                data.frame(t(tme_score_sc)), data.frame(t(stromal_score_sc)),
                all_keyexpr_sc)
meta_expr_sc <- do.call(sc_list, rbind)
results_meta_sc = ConsensusClusterPlus(as.matrix(meta_expr_sc), maxK = 6, reps = 50, pItem = 0.8, pFeature = 1,
                                       title = title, clusterAlg = "hc", distance = "pearson", 
                                       seed = 1262118388.71279, plot = "png")

cluster_meta_sc <- results_meta_sc[[3]][["consensusClass"]]
cluster_meta_sc <- data.frame(sampleID = names(cluster_meta_sc), Cluster = cluster_meta_sc)
cluster_meta_sc %>%
  group_by(Cluster) %>%
  dplyr::summarise(n = n())
cluster_meta_sc$Cluster <- str_c('C' ,cluster_meta_sc$Cluster)
sc_data@meta.data <- merge(sc_data@meta.data, cluster_meta_sc, by.x = 0, by.y = 0)
rownames(sc_data@meta.data) <- sc_data@meta.data$Row.names
sc_data@meta.data <- sc_data@meta.data[,-1]
DimPlot(sc_data, reduction = "tsne", 
        group.by = "Cluster", label.size = 7,
        pt.size = 1,
        cols = color_subtype)+ 
  ggtitle('Cluster')+
  guides(colour = guide_legend(override.aes = list(size = 3)))  +  
  theme(legend.title = element_blank(),
        legend.text = element_text(size = 17), 
        legend.position = 'right',
        plot.title = element_text(size = 27, face = "bold"))+
  guides(colour = guide_legend(override.aes = list(size = 8)),
         size = guide_legend(order=3))+
  theme(axis.text.x = element_text(angle = 0, hjust = 1, size = 20),
        axis.text.y = element_text(angle = 0, hjust = 1, size = 20),
        axis.title.x=element_text(size = 20),
        axis.title.y=element_text(size = 20))

tem_num <- sc_data@meta.data %>%
  group_by(Cluster) %>%
  dplyr::summarise(n = n())
tem_num <- tem_num[order(tem_num$n), ]
tem_num$cell <- 'Cell'
ggplot(tem_num, aes(cell, n, fill = Cluster, group = Cluster))+ theme_classic() + 
  geom_bar(stat="identity",position="fill")+
  ggtitle("") + big +
  xlab('') + ylab('Fraction') + 
  theme( legend.title = element_text(size = 20),
         legend.text = element_text(size = 20),
         title=element_text(size=18))  +
  theme(axis.text.x = element_text( angle = 90)) +
  theme(strip.text.x = element_text(size = 17)) + 
  scale_y_continuous(labels = scales::percent) +
  scale_fill_manual(values = color_subtype)

sc_data@meta.data$Cluster_new <- ifelse(sc_data@meta.data$Cluster == 'C3', 'C3', 'Non-C3')
tem <- sc_data@meta.data %>%
  group_by(Cluster_new,cellType) %>%
  dplyr::summarise(n = n())
ggplot(tem, aes(Cluster_new, n, fill = cellType, group = cellType))+ theme_classic() + 
  geom_bar(stat="identity",position="fill")+
  ggtitle("") + big +
  xlab('') + ylab('Fraction') + 
  theme( legend.title = element_text(size = 20),
         legend.text = element_text(size = 20),
         title=element_text(size=18))  +
  theme(axis.text.x = element_text( angle = 90)) +
  theme(strip.text.x = element_text(size = 17)) + 
  scale_y_continuous(labels = scales::percent) +
  scale_fill_manual(values = sc_color)

sc_s100a9 <- sc_data@assays$SCT@scale.data['S100A9',] %>% data.frame
head(sc_s100a9)
names(sc_s100a9) <- 'S100A9'
sc_s100a9 <- merge(sc_s100a9, sc_data@meta.data, by.x= 0, by.y = 0)
sc_s100a9 %>%
  group_by(Cluster) %>%
  summarise(mean = mean(S100A9))

ggplot(sc_s100a9, map = aes(x=Cluster_new, y = S100A9)) + 
  geom_boxplot(position=position_dodge(0.9), size = 2, width = 0.5, 
               aes(color = Cluster_new, x = Cluster_new, y = S100A9))+ 
  scale_color_manual(values = color_subtype[3:4])+ 
  geom_jitter(width=0.25, aes(color = Cluster_new)) +
  theme_bw()+
  theme(axis.text.x=element_text(angle=90,hjust = 1,colour="black",size=30), 
        axis.text.y=element_text(size=30,face="plain"), 
        axis.title.y=element_text(size = 30,face="bold"),
        legend.text=element_text(face="italic",  colour="black",  
                                 size=30),
        legend.title=element_text(face="bold",  colour="black",
                                  size=30),
        panel.grid.major = element_blank(),   
        panel.grid.minor = element_blank(),
        panel.border = element_rect(size=2, linetype="solid"),
        title=element_text(size=25))+
  ylab("S100A9 expression")+xlab("") + 
  stat_compare_means(aes(group = Cluster_new), label = "p.signif", size = 15,
                     label.y = 7,  label.x = 1.83) +
  theme(legend.position ="None",  title=element_text(size=15)) 

sc_s100a9$roc_value <- ifelse(sc_s100a9$Cluster_new == 'C3', 1, 0)
rocobj <- roc(sc_s100a9$roc_value, sc_s100a9$S100A9)

ggroc(rocobj, colour = mycol[12], size = 1) +
  ggtitle(paste0('ROC Curve ', '(AUC = ', '0.87', ')')) +  theme_classic() + sbig +
  theme(axis.ticks = element_blank()) +
  theme(plot.title=element_text(hjust = 0.5, size= 30, face = 'bold'),
        axis.title.x=element_text(size = 30, face = 'plain'),
        axis.title.y=element_text(size = 30, face = 'plain'),
        axis.text.x=element_text(size = 30,color = "black"),
        axis.text.y=element_text(size = 30,color = "black"))

###feature gene in scRNA-seq
Idents(sc_data) <- sc_data@meta.data$cellType
all.markers <- FindAllMarkers(sc_data, only.pos = FALSE, min.pct = 0.25, logfc.threshold = 0.25)
head(all.markers)
df <- all.markers
df$label <- ifelse(df$p_val_adj<0.01,"adjust P-val<0.01","adjust P-val>=0.01")
names(df)[names(df) == "gene"] <- "geneID"
top10sig <- df %>% group_by(cluster) %>% top_n(n = 5, wt = avg_log2FC)

df$size <- case_when(!(df$geneID %in% top10sig$geneID)~ 1,
                     df$geneID %in% top10sig$geneID ~ 2)
dt <- filter(df,size==1)

bar_hight <- df %>%
  group_by(cluster) %>%
  summarise(max = max(avg_log2FC),
            min = min(avg_log2FC))
dfbar<-data.frame(x= bar_hight$cluster,
                  y= bar_hight$max)
dfbar1<-data.frame(x=bar_hight$cluster,
                   y= bar_hight$min)
p1 <- ggplot()+
  geom_col(data = dfbar,
           mapping = aes(x = x,y = y),
           fill = "#dcdcdc",alpha = 0.6)+
  geom_col(data = dfbar1,
           mapping = aes(x = x,y = y),
           fill = "#dcdcdc",alpha = 0.6)
p1
p2 <- ggplot()+
  geom_col(data = dfbar,
           mapping = aes(x = x,y = y),
           fill = "#dcdcdc",alpha = 0.6)+
  geom_col(data = dfbar1,
           mapping = aes(x = x,y = y),
           fill = "#dcdcdc",alpha = 0.6)+
  geom_jitter(data = dt,
              aes(x = cluster, y = avg_log2FC, color = label),
              size = 0.85,
              width =0.4)+
  geom_jitter(data = top10sig,
              aes(x = cluster, y = avg_log2FC, color = label),
              size = 1,
              width =0.4)
p2

dfcol<-data.frame(x= 1:8,
                  y= 0,
                  label= c('Endothelial', 'CD4+ T', 'Fibroblast', 'Epithelial', 'Stromal', 'CD8+ T', 'Neutrophil', 'Mono/Mac'))
mycol <- c("#E64B357F","#4DBBD57F","#00A0877F","#3C54887F","#F39B7F7F","#8491B47F","#91D1C27F","#DC00007F")
p3 <- p2 + geom_tile(data = dfcol,
                     aes(x=x,y=y),
                     height=0.4,
                     color = "black",
                     fill = mycol,
                     alpha = 0.5,
                     show.legend = F)
p3

p4 <- p3+
  geom_text_repel(
    data=top10sig,
    aes(x=cluster,y=avg_log2FC,label=geneID),
    force = 1.2,
    arrow = arrow(length = unit(0.008, "npc"),
                  type = "open", ends = "last")
  )
p4

p5 <- p4 +
  scale_color_manual(name=NULL,
                     values = c("#FF6633","#996633"))
p5

p6 <- p5+
  labs(x="Cluster",y="average logFC")+
  geom_text(data=dfcol,
            aes(x=x,y=y,label=label),
            size =4,
            color ="white")
p6

p7 <- p6+
  theme_minimal()+
  theme(
    axis.title = element_text(size = 16,
                              color = "black",
                              face = "bold"),
    axis.line.y = element_line(color = "black",
                               size = 1.2),
    axis.line.x = element_blank(),
    axis.text.x = element_blank(),
    panel.grid = element_blank(),
    legend.position = "top",
    legend.direction = "vertical",
    legend.justification = c(1,0),
    legend.text = element_text(size = 15),
    axis.text.y=element_text(size = 15)
  )+
  guides(colour = guide_legend(override.aes = list(size = 5)))  
p7
