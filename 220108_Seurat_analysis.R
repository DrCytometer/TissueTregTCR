library(SingleCellExperiment)
library(ggplot2)
library(pheatmap)
library(ggsci)
library(scales)
library(scater)
library(Matrix)
library(viridis)
library(Seurat)
library(dplyr)

setwd("//central-cluster/algroup/Orian Bricard/201130_ORB596_Oliver_2019_Tissue_treg_10x_GE_VDJ_FB_reanalysis/R/220108_Seurat_analysis")

## Define Colors ####

#show_col(pal_npg("nrc")(10))
mouseId_color = c("ALFX24_1b" = pal_npg("nrc")(10)[1], "ALFX24_1c" = pal_npg("nrc")(10)[2], "ALFX24_1d" = pal_npg("nrc")(10)[3], "ALFX24_1e" = pal_npg("nrc")(10)[4])
tissue_origin_color = c("Blood" = pal_npg("nrc")(10)[8],"Kidney"= pal_npg("nrc")(10)[6], "LPL" = pal_npg("nrc")(10)[5], "Liver" = "#A52A2AFF", "Pancreas"  = pal_npg("nrc")(10)[10] )



## Load data ####
sce = readRDS("//central-cluster/algroup/Orian Bricard/201130_ORB596_Oliver_2019_Tissue_treg_10x_GE_VDJ_FB_reanalysis/R/201202_data_preprocessing/tgtfsce.rds")

so = as.Seurat(sce, data = "counts")

## Seurat QC and cleaning ####

rownames(so)
so[["percent.mt"]] <- PercentageFeatureSet(so, pattern = "^mt-")
head(so@meta.data)

VlnPlot(so, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
FeatureScatter(so, feature1 = "nFeature_RNA", feature2 = "percent.mt")
so <- subset(so, subset = nFeature_RNA > 200  & percent.mt < 5)

## Normalisation ####

so <- NormalizeData(so, normalization.method = "LogNormalize", scale.factor = 10000)
so[["originalexp"]]@data

## Identification of highly variable features ####

so <- FindVariableFeatures(so, selection.method = "vst", nfeatures = 2000)

# Identify the 10 most highly variable genes
top10 <- head(VariableFeatures(so), 10)

# plot variable features with and without labels
plot1 <- VariableFeaturePlot(so)
plot2 <- LabelPoints(plot = plot1, points = top10, repel = TRUE)
plot2

## Scaling data ####
all.genes <- rownames(so)
so <- ScaleData(so, features = all.genes)
so[["originalexp"]]@scale.data

## Perform linear dimension reduction ####
so <- RunPCA(so, features = VariableFeatures(object = so))
print(so[["pca"]], dims = 1:5, nfeatures = 5)
VizDimLoadings(so, dims = 1:2, reduction = "pca")
DimPlot(so, reduction = "pca")
DimHeatmap(so, dims = 1, cells = 500, balanced = TRUE)

## Determine dimensionality of the dataset #####
so <- JackStraw(so, num.replicate = 100)
so <- ScoreJackStraw(so, dims = 1:20)
JackStrawPlot(so, dims = 1:20)

ElbowPlot(so)
### => choose 20, not 100% confident ###

## Cluster the cells ####

so <- FindNeighbors(so, dims = 1:20)
so <- FindClusters(so, resolution = 0.5)

head(Idents(so), 5)

## UMAP ####

so <- RunUMAP(so, dims = 1:20, seed.use = 1)
DimPlot(so, reduction = "umap",  label = TRUE) + NoLegend()
ggsave("10_clusters_UMAP.pdf", width = 100, height = 100, units = "mm")

## tSNE ####

so <- RunTSNE(so, dims = 1:20, seed.use = 1)
DimPlot(so, reduction = "tsne",  label = TRUE) + NoLegend()
ggsave("10_clusters_TSNE.pdf", width = 100, height = 100, units = "mm")

DimPlot(so, reduction = "tsne", group.by = "tissue_origin")
ggsave("tissue_origin_TSNE.pdf", width = 120, height = 100, units = "mm")

# alternative
#tsne_df = so@reductions$tsne@cell.embeddings %>% 
#  as.data.frame() %>% 
#  cbind(so@meta.data[,c("seurat_clusters", "tissue_origin", "mouseId")])
#  cbind(tissue_origin = so@meta.data$tissue_origin)
#ggplot(data = tsne_df)+
#  geom_point(aes(x=tSNE_1, y = tSNE_2, color = tissue_origin))+
#  theme_bw()


df = so@meta.data[,c("seurat_clusters", "tissue_origin")] %>%
  group_by(seurat_clusters, tissue_origin) %>%
  summarise(count =n()) %>%
  as.data.frame()

ggplot(data = df)+
  geom_bar(stat = "identity", aes(x = seurat_clusters, y = count, fill = tissue_origin), position = position_dodge())

ggplot(data = df)+
  geom_bar(stat = "identity", aes(x = tissue_origin, y = count, fill = seurat_clusters), position = position_dodge())

df = df  %>%
  group_by(tissue_origin) %>%
  mutate(tissue_origin_freq = count / sum(count)) %>%
  as.data.frame()

df = df  %>%
  group_by(seurat_clusters) %>%
  mutate(seurat_clusters_freq = count / sum(count)) %>%
  as.data.frame()

ggplot(data = df)+
  scale_x_discrete(expand=c(0,0))+
  scale_y_discrete(expand=c(0,0))+
  # scale_fill_viridis(trans = "log10")+
  scale_fill_gradient(trans = "log10", low = "white", high = "red") +
  geom_tile( aes(x = tissue_origin, y = seurat_clusters, fill = count)) +
  theme_bw()+
  theme(panel.grid = element_blank())
ggsave("cluster_vs_tissue_origin_count.pdf", width = 110, height = 100, units = "mm")

ggplot(data = df)+
  ggtitle("Tissue origin proportion heatmap with absolute count number")+
  scale_x_discrete(expand=c(0,0))+
  scale_y_discrete(expand=c(0,0))+
  # scale_fill_viridis(trans = "log10")+
  scale_fill_gradient(low = "white", high = "red", limits = c(0,100)) +
  geom_tile( aes(x = tissue_origin, y = seurat_clusters, fill = tissue_origin_freq*100)) +
  geom_text(aes(x = tissue_origin, y = seurat_clusters, label = count))+
  theme_bw()+
  theme(panel.grid = element_blank(),
        legend.title=element_blank())
ggsave("cluster_vs_tissue_origin_tissue_origin_freq.pdf", width = 110, height = 100, units = "mm")

ggplot(data = df)+
  ggtitle("Cluster proportion heatmap with absolute count number")+
  scale_x_discrete(expand=c(0,0))+
  scale_y_discrete(expand=c(0,0))+
  # scale_fill_viridis(trans = "log10")+
  scale_fill_gradient(low = "white", high = "red", limits = c(0,100)) +
  geom_tile( aes(x = tissue_origin, y = seurat_clusters, fill = seurat_clusters_freq*100)) +
  geom_text(aes(x = tissue_origin, y = seurat_clusters, label = count))+
  theme_bw()+
  theme(panel.grid = element_blank(),
        legend.title=element_blank())
ggsave("cluster_vs_tissue_origin_seurat_clusters_freq.pdf", width = 110, height = 100, units = "mm")




## Find cluster Markers ####

so.markers <- FindAllMarkers(so, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)

so.markers %>%
  group_by(cluster) %>%
  top_n(n = 5, wt = avg_log2FC) -> top10
DoHeatmap(so, features = top10$gene) + NoLegend()
ggsave("top5marker_for_clusters.pdf", width = 297, height = 210, units = "mm")



sessionInfo()




 