#---- General ------------------------------------------------------------------
# Análisis scRNAseq Spatial: el objetivo es establecer un flujo de trabajo
# para analizar datos scRNAseq y posteiormente integrarlos con datos espaciales
# Programador: J.Raul.Mtz.Valderrama
# Fecha: Febrero 3 de 2025
# R version: 4.31.


#---- Data ---------------------------------------------------------------------

# Primero debo realizar el análisis de una muestra scRNAseq para montar el 
# flujo de trabajo y posteriormente integrarlos con datos espaciales

# Establecer directorio de trabajo
setwd("D:/marval_windows/JR_MARVAL/est_INMEGEN/analyses/scrnaseq")
list.files()

# Instalar liberias
#install.packages("Seurat")
library(Seurat)
#install.packages("dylyr")
library(dplyr)

# Cargar los archivos
# Se genera una matriz genes*celulas con los valores de expresion
counts <- ReadMtx(mtx = "count_matrix_sparse.mtx",
                  features = "count_matrix_genes.tsv",
                  cells = "count_matrix_barcodes.tsv",
                  feature.column = 1) # La matriz de los genes solo tiene 1 columna

# Crear el objeto Seurat
seurat_object <- CreateSeuratObject(counts = counts,
                                    project = "scRNAseq_BreastCancer")

# Cargar el metadata
metadata <- read.csv("metadata.csv", row.names = 1)
seurat_object <- AddMetaData(seurat_object, metadata)

# Ver información del objeto
seurat_object


#---- Control de calidad -------------------------------------------------------

# Ver métricas básicas
# Detecta genes mitocondriales
seurat_object[["percent.mt"]] <- PercentageFeatureSet(seurat_object, 
                                                      pattern = "^MT-")  

# Visualizar QC con Violin Plots
VlnPlot(seurat_object, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"),
        ncol = 3)

# Mostrar estadísticos de QC
summary(seurat_object$nFeature_RNA)  # Genes detectados por célula
summary(seurat_object$nCount_RNA)    # UMI totales por célula
summary(seurat_object$percent.mt)    # Porcentaje de genes mitocondriales

# Aplicar filtros: definir umbrales (ajústalos según tu dataset)
dim(seurat_object)

seurat_object <- subset(seurat_object, 
                        subset = nFeature_RNA > 500 & nFeature_RNA < 4000 & 
                          percent.mt < 10)

# Revisar cuántas células quedan después del filtrado
dim(seurat_object)

# Visualizar QC con Violin Plots
VlnPlot(seurat_object, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"),
        ncol = 3)

# Relación entre genes detectados y UMIs por célula
FeatureScatter(seurat_object, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")

# Relación entre genes detectados y porcentaje mitocondrial
FeatureScatter(seurat_object, feature1 = "nCount_RNA", feature2 = "percent.mt")



#---- Normalización ------------------------------------------------------------

seurat_object <- NormalizeData(seurat_object,
                               normalization.method = "LogNormalize", 
                               scale.factor = 10000)

# Cuentas originales vs. normalizadas
# Cuentas originales log-transformadas
seurat_object$log_counts <- log1p(seurat_object$nCount_RNA)

# Visualización con gráficos de violín
VlnPlot(seurat_object, features = c("log_counts", "nCount_RNA"), 
        group.by = "orig.ident", 
        pt.size = 0.1, 
        split.by = "orig.ident")

# Comparar la distribución de características clave 
#(genes detectados, cuentas y log-transformed counts)
VlnPlot(seurat_object, features = c("nFeature_RNA", "nCount_RNA", "log_counts"),
        ncol = 3)
FeatureScatter(seurat_object, feature1 = "nCount_RNA", feature2 = "log_counts")


#---- Identifiación de genes altamente variables -------------------------------

seurat_object <- FindVariableFeatures(seurat_object, 
                                      selection.method = "vst", 
                                      nfeatures = 1000)

# Plot de variabilidad
VariableFeaturePlot(seurat_object)
# Identifiacion de genes top
top <- head(VariableFeatures(seurat_object),8)
LabelPoints(plot = VariableFeaturePlot(seurat_object), points = top, 
            repel = TRUE,
            max.overlaps = 20, xnudge = 0.01, ynudge = 0.01)


#---- Escalado -----------------------------------------------------------------

# Para todos los genes de alta variabilidad
seurat_object <- ScaleData(seurat_object, 
                           features = VariableFeatures(seurat_object))
DoHeatmap(seurat_object, features = head(VariableFeatures(seurat_object), 30))


# Para todos los genes del análisis
all.genes <- rownames(seurat_object)
seurat_object <- ScaleData(seurat_object, features = all.genes)
DoHeatmap(seurat_object, features = sample(all.genes, 30))


#---- Reducción de Dimensionalidad ---------------------------------------------

seurat_object <- RunPCA(seurat_object, 
                        features = VariableFeatures(seurat_object))
# Exploracion de los datos
print(seurat_object[["pca"]], dims = 1:5, nfeatures = 5)
# Plot de carga
VizDimLoadings(seurat_object, dims = 1:2, reduction = "pca")
# Vista del PCA
DimPlot(seurat_object, reduction = "pca")# + NoLegend()

# Como saber cuantos PC seleccionar?
DimHeatmap(seurat_object, dims = 1, cells = 500, balanced = TRUE)
DimHeatmap(seurat_object, dims = 1:5, cells = 500, balanced = TRUE)

# Para determinar cuántos PC usar
ElbowPlot(seurat_object, ndims = 50)  

# Cs significativos
seurat_object <- JackStraw(seurat_object, num.replicate = 100)
seurat_object <- ScoreJackStraw(seurat_object, dims = 1:20)
JackStrawPlot(seurat_object, dims = 1:20)

# Umbral de varianza
# Obtener la desviación estándar de cada PC
stdev <- seurat_object[["pca"]]@stdev
# Seleccionar PCs donde la varianza explicada sea mayor a 1
num_pcs <- sum(stdev > 1)
cat("Número óptimo de PCs:", num_pcs, "\n")

# Varianza explicada
# Visualiza la proporción de varianza explicada por cada componente principal
seurat_object[["pca"]]@stdev
# Calcula la proporción de varianza explicada
var_explained <- seurat_object[["pca"]]@stdev^2 / sum(seurat_object[["pca"]]@stdev^2)
# Visualiza la varianza explicada acumulada
cumsum(var_explained)
# Encuentra el número de PCs que explican el 90% de la varianza
pc_90 <- which(cumsum(var_explained) >= 0.90)[1]
pc_90


#---- Clustering ---------------------------------------------------------------

# Construye las redes de vecinos, considerando  el # de PC a retener
seurat_object <- FindNeighbors(seurat_object, dims = 1:50)
# Determina el número de clusters
seurat_object <- FindClusters(seurat_object, resolution = 0.5)
# Cluster IDs of the first 5 cells
head(Idents(seurat_object), 5)

# Realizar la reducción de dimensiones con UMAP
seurat_object <- RunUMAP(seurat_object, dims = 1:50)
# Visualizar los clusters en un gráfico UMAP
DimPlot(seurat_object, reduction = "umap")#, group.by = "seurat_clusters")

# Realizar la reducción de dimensiones con t-SNE
seurat_object <- RunTSNE(seurat_object, dims = 1:50)
# Visualizar los clusters en un gráfico t-SNE
DimPlot(seurat_object, reduction = "tsne")#, group.by = "seurat_clusters")

# Salvar resultado
saveRDS(seurat_object, file = "sc_rnaseq_test.rds")

#---- Finding Clusters Biomarkers ----------------------------------------------

# Encontrar marcadores para un cluster determinado
cluster2.markers <- FindMarkers(seurat_object, ident.1 = 2)
head(cluster2.markers, n = 5)

# Encontrar marcadores para un un grupo en determinados clusters
cluster5.markers <- FindMarkers(seurat_object, ident.1 = 5, ident.2 = c(0, 3))
head(cluster5.markers, n = 5)

# Encontrar marcadores para cada cluster comparando contra todas las céluas
pbmc.markers <- FindAllMarkers(seurat_object, only.pos = TRUE)
pbmc.markers %>%
  group_by(cluster) %>%
  dplyr::filter(avg_log2FC > 1)

# Se puede seleccionar la prueba estadistica a realizar para DEG
cluster0.markers <- FindMarkers(seurat_object, ident.1 = 0, 
                                logfc.threshold = 0.25, 
                                test.use = "roc", only.pos = TRUE)

# Visualizar los marcadores
VlnPlot(seurat_object, features = c("IL7R", "CCR7"))
# Plot raw counts
VlnPlot(seurat_object, features = c("IL7R", "CCR7"), slot = "counts", log = TRUE)
# Vusaluzar en el UMAP/TSNE
FeaturePlot(seurat_object, features = c("IL7R", "CCR7"))
# Heatmap
pbmc.markers %>%
  group_by(cluster) %>%
  dplyr::filter(avg_log2FC > 1) %>%
  slice_head(n = 3) %>%
  ungroup() -> top10

DoHeatmap(seurat_object, features = top10$gene) + NoLegend()

# Opcion 2

# Identificar los genes diferencialmente expresados en cada cluster
cluster_markers <- FindAllMarkers(seurat_object, 
                                  only.pos = TRUE, min.pct = 0.25, 
                                  logfc.threshold = 0.25)

# Visualizar los primeros marcadores encontrados
head(cluster_markers)

library(dplyr)

# Seleccionar los 10 genes más diferencialmente expresados en cada cluster
top10 <- cluster_markers %>% group_by(cluster) %>% top_n(10, avg_log2FC)
# Visualizar la tabla de genes seleccionados
top10

# Visualizacion 
DoHeatmap(seurat_object, features = top10$gene)# + scale_fill_viridis()
VlnPlot(seurat_object, features = c("CCR7"), pt.size = 0.1)
FeaturePlot(seurat_object, features = c("CCR7", "IL7R"))

#---- Asignacion de tipos celulares --------------------------------------------
new.cluster.ids <- c("Naive CD4 T", "CD14+ Mono", "Memory CD4 T", "B", 
                     "CD8 T", "FCGR3A+ Mono",
                     "NK", "DC", "Platelet",
                     rep("Unknown", 15))

names(new.cluster.ids) <- levels(seurat_object)
seurat_object <- RenameIdents(seurat_object, new.cluster.ids)
DimPlot(seurat_object, reduction = "umap", 
        label = TRUE, pt.size = 0.5) + NoLegend()

library(ggplot2)
plot <- DimPlot(seurat_object, reduction = "umap", label = TRUE, label.size = 4.5) 
  + xlab("UMAP 1") + ylab("UMAP 2") +
  theme(axis.title = element_text(size = 18), 
        legend.text = element_text(size = 18))
        + guides(colour = guide_legend(override.aes = list(size = 10)))
plot



# Salvar resultado
saveRDS(seurat_object, file = "sc_rnaseq_test.rds")

#Metodo supervisado
head(cluster_markers)

# Semisupervisado
BiocManager::install("SingleR")
BiocManager::install("celldex")
library(SingleR)
library(Celldex)

ref <- celldex::MonacoImmuneData()  # Base de datos de células inmunes
test <- GetAssayData(seurat_object, slot = "data")
pred <- SingleR(test = test, ref = ref, labels = ref$label.main)

# Asigna las etiquetas a los clusters
seurat_object$SingleR.labels <- pred$labels
DimPlot(seurat_object, reduction = "umap", group.by = "SingleR.labels", label = TRUE)


#---- Enriquecimiento de vias --------------------------------------------------

# Selecciona los genes con p_val ajustado < 0.05 y log2FC > 0.25
top_genes <- cluster_markers %>%
  filter(p_val_adj < 0.05 & avg_log2FC > 0.25) %>%
  pull(gene) # Extrae solo los nombres de los genes

# Instala clusterProfiler si no lo tienes
if (!requireNamespace("clusterProfiler", quietly = TRUE)) {
  install.packages("BiocManager")
  BiocManager::install("clusterProfiler")
}

# Cargar el paquete
library(clusterProfiler)

# Convertir nombres de genes a ENSEMBL si es necesario
library(org.Hs.eg.db) # Base de datos para humanos
genes_entrez <- bitr(top_genes, fromType="SYMBOL", toType="ENTREZID", OrgDb="org.Hs.eg.db")

# Hacer el enriquecimiento GO
go_enrichment <- enrichGO(gene         = genes_entrez$ENTREZID,
                          OrgDb        = org.Hs.eg.db,
                          keyType      = "ENTREZID",
                          ont          = "BP",  # Ontología: "BP" (procesos biológicos)
                          pAdjustMethod = "BH",
                          pvalueCutoff  = 0.05,
                          qvalueCutoff  = 0.05)

# Visualizar los resultados
dotplot(go_enrichment, showCategory = 10) # Muestra los 10 procesos más significativos

# Convertir genes a formato KEGG
kegg_enrichment <- enrichKEGG(gene         = genes_entrez$ENTREZID,
                              organism     = 'hsa', # Homo sapiens
                              pAdjustMethod = "BH",
                              pvalueCutoff  = 0.05)


#--- Referencias: --------------------------------------------------------------

# https://satijalab.org/seurat/articles/pbmc3k_tutorial#run-non-linear-dimensional-reduction-umaptsne
# https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSM5354531
# https://zenodo.org/records/4739739
# https://www.nature.com/articles/s41588-021-00911-1#Abs1

sessionInfo()










