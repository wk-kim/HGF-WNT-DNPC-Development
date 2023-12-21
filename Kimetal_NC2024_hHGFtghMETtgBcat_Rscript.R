#Add necessary tools to library
library(Seurat)
library(devtools)
library(dplyr)
library(Matrix)
library(cowplot)
library(ggplot2)
library(reticulate)
library(RColorBrewer)
library(ggpubr)
library(GGally)
library(Seurat)
library(devtools)
library(R.utils)
library(SeuratWrappers)
library(monocle3)
library(ggplot2)
library(dplyr)
library(BiocManager)
library(remotes) 

####Loading data####
HGFMETunfiltered.data <- Read10X("//isi-dcnl/user_data/zjsun/seq/201120_TGen/RUN_mm10_hMETtg_hHGFtg/count_39553_mm10_hMETtg_hHGFtg/outs/filtered_feature_bc_matrix")
HGFMETunfiltered <- CreateSeuratObject(counts = HGFMETunfiltered.data,  min.cells = 3, min.features = 200, project = "HGF-MET")
HGFMETunfiltered <- NormalizeData(HGFMETunfiltered)

####Initial processing, Filtering and Clustering####
#HGFMETunfiltered
HGFMETunfiltered[["percent.mt"]] <- PercentageFeatureSet(HGFMETunfiltered, pattern = "^mt-")
tiff(file = "HGFMET Pre-filter Vln.tiff", width = 6, height = 4, units = "in", compression = "lzw", res = 800)
VlnPlot(HGFMETunfiltered, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3, pt.size = 0.1)
dev.off()
tiff(file = "HGFMET Pre-filter nFeature Hist.tiff", width = 4, height = 4, units = "in", compression = "lzw", res = 800)
hist(HGFMETunfiltered@meta.data$nFeature_RNA, breaks = 100, col = "pink", xlab = "nFeature_RNA", main = "HGFMET Pre-filteration")
dev.off()
tiff(file = "HGFMET Pre-filter percentMT Hist.tiff", width = 4, height = 4, units = "in", compression = "lzw", res = 800)
hist(HGFMETunfiltered@meta.data$percent.mt, breaks = 100, col = "pink", xlab = "percent.mt", main = "HGFMET Pre-filteration")
dev.off()
#HGFMETfiltered
HGFMET <- subset(HGFMETunfiltered, subset = nFeature_RNA > 700 & nFeature_RNA < 7000 & percent.mt < 10)
tiff(file = "HGFMET Post-filter Vln.tiff", width = 6, height = 4, units = "in", compression = "lzw", res = 800)
VlnPlot(HGFMET, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3, pt.size = 0.1)
dev.off()
tiff(file = "HGFMET Post-filter nFeature Hist.tiff", width = 4, height = 4, units = "in", compression = "lzw", res = 800)
hist(HGFMET@meta.data$nFeature_RNA, breaks = 100, col = "lightblue", xlab = "nFeature_RNA", main = "HGFMET Post-filteration")
dev.off()
tiff(file = "HGFMET Post-filter percentMT Hist.tiff", width = 4, height = 4, units = "in", compression = "lzw", res = 800)
hist(HGFMET@meta.data$percent.mt, breaks = 100, col = "lightblue", xlab = "percent.mt", main = "HGFMET Post-filteration")
dev.off()

#Genes and UMI counts per cell
mean(HGFMET$nCount_RNA)
mean(HGFMET$nFeature_RNA)

#Clustering
HGFMET <- FindVariableFeatures(HGFMET, selection.method = "vst", nfeatures = 5000)
tiff(file = "HGFMET Variable Genes.tiff", width = 6, height = 4, units = "in", compression = "lzw", res = 800)
VariableFeaturePlot(HGFMET)
dev.off()

HGFMET <- ScaleData(HGFMET, verbose = FALSE)
HGFMET <- RunPCA(HGFMET, npcs = 50, verbose = FALSE)
tiff(file = "HGFMET ElbowPlot.tiff", width = 4, height = 4, units = "in", compression = "lzw", res = 800)
ElbowPlot(HGFMET, ndims = 50)
dev.off()

HGFMET <- FindNeighbors(HGFMET, reduction = "pca", dims = 1:20)
HGFMET <- FindClusters(HGFMET, resolution = 0.5)
HGFMET <- RunUMAP(HGFMET, reduction = "pca", dims = 1:20)
DimPlot(HGFMET, reduction = "umap", pt.size = 0.3) 

#Cell cycle assignment
mouse_cell_cycle_genes <- readRDS(file = "//isi-dcnl/user_data/zjsun/group/Won Kyung Kim/ARKOvCtrl_3timepoint/E18.5/mouse_cell_cycle_genes.rds")
s.genes <- mouse_cell_cycle_genes$s.genes
g2m.genes <- mouse_cell_cycle_genes$g2m.genes

DefaultAssay(HGFMET) <- "RNA"
all.genes <- rownames(HGFMET)
HGFMET <- ScaleData(HGFMET, features = all.genes)
HGFMET <- CellCycleScoring(HGFMET, s.features = s.genes, g2m.features = g2m.genes, set.ident = TRUE)

Idents(object = HGFMET) <- "Phase"
DimPlot(HGFMET, reduction = "umap")
tiff(file = "HGFMET Cell Cyle UMAP.tiff", width = 6.5, height = 6, units = "in", compression = "lzw", res = 800)
DimPlot(HGFMET, reduction = "umap", pt.size = 0.3, cols = c("goldenrod2", "grey75", "deepskyblue2"))
dev.off()

#SFig.2c
Idents(object = HGFMET) <- "stim"
tiff(file = "HGFMET grey UMAP.tiff", width = 6.5, height = 6, units = "in", compression = "lzw", res = 800)
DimPlot(HGFMET, reduction = "umap", pt.size = 0.3, cols = "grey")
dev.off()

#Cell Cycle regression
HGF_MET1 <- HGFMET
DefaultAssay(HGF_MET1) <- "RNA"
HGF_MET1 <- ScaleData(HGF_MET1, vars.to.regress = c("S.Score", "G2M.Score"), features = rownames(HGF_MET1))
HGF_MET1 <- RunPCA(HGF_MET1, features = VariableFeatures(HGF_MET1))
ElbowPlot(HGF_MET1, ndims = 50)

HGF_MET1 <- FindNeighbors(HGF_MET1, reduction = "pca", dims = 1:20)
HGF_MET1 <- FindClusters(HGF_MET1, resolution = 1.5)
HGF_MET1 <- RunUMAP(HGF_MET1, reduction = "pca", dims = 1:20)

#SFig.2c
Idents(object = HGF_MET1) <- "Phase"
tiff(file = "HGF_MET1 Cell Cyle UMAP after Cell Cycle Regression.tiff", width = 6.5, height = 6, units = "in", compression = "lzw", res = 800)
DimPlot(HGF_MET1, reduction = "umap", pt.size = 0.3, cols = c("goldenrod2", "grey75", "deepskyblue2"))
dev.off()

#FeaturePlots for cell markers
DefaultAssay(HGF_MET1)<-"RNA"
FeaturePlot(HGF_MET1, reduction = "umap", features = c("hHGFtg", "hMETtg", "Ar", "Pbsn",
                                                       "Krt5", "Trp63", "Krt8", "Cd24a", 
                                                       "Fbln1", "Myh11", "Plp1", "Pecam1",
                                                       "Rgs5", "Tyrobp", "Ccl5", "Pate4"), cols = c("light grey", "red"), pt.size = 0.3, max.cutoff = "q90")

#Cell Type Identification
Idents(object = HGF_MET1) <- "seurat_clusters"
HGF_MET1 <- RenameIdents(object = HGF_MET1, 
                         '23'="BE", '24' = "BE", '0'= "BE", '4'="BE",
                         '21'="LE", '10'="LE",'5'="LE",'1'="LE",
                         '3'="LE",'2'="LE",'20'="LE", '7'="LE", '16'="LE", '6'="LE", '11'="LE", '13'="LE",'19'="SV",
                         '15'="FB", '8'="FB",'17'="SM", '22'="Pericyte",
                         '25'="Glia",'14'="VE", '12'="Immune", '18'="Immune",
                         '9'="Immune")  
HGF_MET1[["CellTypes"]] <- Idents(object = HGF_MET1)

#SFig.2d
tiff(file = "HGF_MET1 CellTypes UMAP.tiff", width = 6.5, height = 6, units = "in", compression = "lzw", res = 800)
DimPlot(HGF_MET1, reduction = "umap", pt.size = 0.3, cols = c("chartreuse3", "salmon", "darkslategray3", "plum4", "brown3", "blueviolet", "blue", "bisque3", "steelblue1"))
dev.off()

#SFig.2e
DefaultAssay(HGF_MET1) <- "RNA"
tiff(file = "HGF_MET1 hHGFtg 3.0 expression plots.tiff", width = 6.5, height = 6, units = "in", compression = "lzw", res = 800)
FeaturePlot(HGF_MET1, reduction = "umap", features = c("hHGFtg"), cols = c("light grey", "red"), pt.size = 0.3, max.cutoff = 3.0)
dev.off()
tiff(file = "HGF_MET1 hMETtg 3.0 expression plots.tiff", width = 6.5, height = 6, units = "in", compression = "lzw", res = 800)
FeaturePlot(HGF_MET1, reduction = "umap", features = c("hMETtg"), cols = c("light grey", "red"), pt.size = 0.3, max.cutoff = 3.0)
dev.off()
tiff(file = "HGF_MET1 Pbsn 3.0 expression plots.tiff", width = 6.5, height = 6, units = "in", compression = "lzw", res = 800)
FeaturePlot(HGF_MET1, reduction = "umap", features = c("Pbsn"), cols = c("light grey", "red"), pt.size = 0.3, max.cutoff = 3.0)
dev.off()
tiff(file = "HGF_MET1 Ar 3.0 expression plots.tiff", width = 6.5, height = 6, units = "in", compression = "lzw", res = 800)
FeaturePlot(HGF_MET1, reduction = "umap", features = c("Ar"), cols = c("light grey", "red"), pt.size = 0.3, max.cutoff = 3.0)
dev.off()
tiff(file = "HGF_MET1 Epcam 3.0 expression plots.tiff", width = 6.5, height = 6, units = "in", compression = "lzw", res = 800)
FeaturePlot(HGF_MET1, reduction = "umap", features = c("Epcam"), cols = c("light grey", "red"), pt.size = 0.3, max.cutoff = 3.0)
dev.off()
tiff(file = "HGF_MET1 Vim 3.0 expression plots.tiff", width = 6.5, height = 6, units = "in", compression = "lzw", res = 800)
FeaturePlot(HGF_MET1, reduction = "umap", features = c("Vim"), cols = c("light grey", "red"), pt.size = 0.3, max.cutoff = 3.0)
dev.off()

#SFig.2f
Idents(object = HGF_MET1) <- "CellTypes"
table(Idents(HGF_MET1))

#Degs celltype clusters
DefaultAssay(HGF_MET1) <- "RNA"
Idents(object = HGF_MET1) <- "CellTypes"
all.genes <- rownames(HGF_MET1)
HGF_MET1 <- ScaleData(HGF_MET1, features = all.genes)
HGF_MET1.allmarkers <- FindAllMarkers(HGF_MET1, min.pct = 0.25, logfc.threshold = 0.25, only.pos=TRUE)
write.csv(HGF_MET1.allmarkers, file = "HGF_MET1.allmarkers.csv")

#SFig.2g
Idents(object = HGF_MET1) <- "CellTypes"
tiff(file = "HGF_MET1 markers DotPlot.tiff", width =12 , height = 4, units = "in", compression = "lzw", res = 800)
DotPlot(HGF_MET1, features = c("hMETtg", "hHGFtg",
                               "Krt14", "Lgals7", "Krt5", "Aqp3", "Col17a1",
                               "Gm5615", "Agr2", "5430419D17Rik", "Oit1", "Azgp1",
                               "Defb42", "Elf5", "Hoxd8", "Svs4", "Pate4",
                               "Apod", "Fbln1", "Dcn", "Crispld2", "Penk",
                               "1500015O10Rik", "Dkk2", "Itgbl1", "Mfap5", "Actg2",
                               "Rgs5", "Ndufa4l2", "Vtn", "Cox4i2", "Kcnj8",
                               "Plp1", "Fabp7", "Cdh19", "Kcna1", "Gfra3", 
                               "Flt1", "Plvap", "Aqp1", "Pecam1", "Cdh5", 
                               "Ccl5", "Rgs1", "Il1b", "Fcer1g", "Tyrobp"
),cols = c("light grey", "red")) + RotatedAxis()
dev.off()

####Re-clustering Epi####
#Fig.2e
Idents(object = HGF_MET1) <- "CellTypes"
tiff(file = "HGF_MET1 Epi Highlight UMAP.tiff", width = 6.5, height = 6, units = "in", compression = "lzw", res = 800)
DimPlot(HGF_MET1, reduction = "umap", pt.size = 0.3, cols = c("chartreuse3", "salmon", "lightgrey", "lightgrey", "lightgrey", "lightgrey"
                                                              , "lightgrey", "lightgrey", "lightgrey"))
dev.off()

#Subset epithelial cells
Idents(object = HGF_MET1) <- "CellTypes"
HGF_MET1.epi <- subset(HGF_MET1, idents = c("BE","LE"))

#Run the standard workflow for visualization and clustering
Idents(object = HGF_MET1.epi) <- "seurat_clusters"
DefaultAssay(HGF_MET1.epi) <- "RNA"
HGF_MET1.epi <- ScaleData(HGF_MET1.epi, verbose = FALSE)
HGF_MET1.epi <- RunPCA(HGF_MET1.epi, npcs = 50, verbose = FALSE)
ElbowPlot(HGF_MET1.epi, ndims = 50)

#Umap and Clustering
HGF_MET1.epi <- FindNeighbors(HGF_MET1.epi, reduction = "pca", dims = 1:18)
HGF_MET1.epi <- FindClusters(HGF_MET1.epi, resolution = 0.5)
HGF_MET1.epi <- RunUMAP(HGF_MET1.epi, reduction = "pca", dims = 1:18)
DimPlot(HGF_MET1.epi, reduction = "umap", pt.size = 0.3, label = TRUE)

#SFig.2h
Idents(object = HGF_MET1.epi) <- "stim"
tiff(file = "HGF_MET1.epi grey UMAP.tiff", width = 6.5, height = 6, units = "in", compression = "lzw", res = 800)
DimPlot(HGF_MET1.epi, reduction = "umap", pt.size = 0.3, cols = "grey")
dev.off()

#Cell cycle assignment
DefaultAssay(HGF_MET1.epi) <- "RNA"
all.genes <- rownames(HGF_MET1.epi)
HGF_MET1.epi <- ScaleData(HGF_MET1.epi, features = all.genes)
HGF_MET1.epi <- CellCycleScoring(HGF_MET1.epi, s.features = s.genes, g2m.features = g2m.genes, set.ident = TRUE)
Idents(object = HGF_MET1.epi) <- "Phase"
DimPlot(HGF_MET1.epi, reduction = "umap", pt.size = 0.3, cols = c("goldenrod2", "grey75", "deepskyblue2"))

#Cell Cycle Regression
HGF_MET1.epi1 <- HGF_MET1.epi
DefaultAssay(HGF_MET1.epi1) <- "RNA"
HGF_MET1.epi1 <- ScaleData(HGF_MET1.epi1, vars.to.regress = c("S.Score", "G2M.Score"), features = rownames(HGF_MET1.epi1))
HGF_MET1.epi1 <- RunPCA(HGF_MET1.epi1, features = VariableFeatures(HGF_MET1.epi1))
ElbowPlot(HGF_MET1.epi1, ndims = 30)

HGF_MET1.epi1 <- FindNeighbors(HGF_MET1.epi1, reduction = "pca", dims = 1:15)
HGF_MET1.epi1 <- FindClusters(HGF_MET1.epi1, resolution = 0.3)
HGF_MET1.epi1 <- RunUMAP(HGF_MET1.epi1, reduction = "pca", dims = 1:15)

#SFig.2h
Idents(object = HGF_MET1.epi1) <- "Phase"
tiff(file = "HGF_MET1.epi1 Cell Cyle UMAP after regression.tiff", width = 6.5, height = 6, units = "in", compression = "lzw", res = 800)
DimPlot(HGF_MET1.epi1, reduction = "umap", pt.size = 0.3, cols = c("goldenrod2", "grey75", "deepskyblue2"))
dev.off()

#Rename
Idents(object = HGF_MET1.epi1) <- "seurat_clusters"
DimPlot(HGF_MET1.epi1, reduction = "umap", pt.size = 0.3, label = TRUE)
HGF_MET1.epi1 <- RenameIdents(object = HGF_MET1.epi1, '0' = "BE", '5' = "LE1", '4' = "LE2", '3' = "LE3", '11' = "LE4", 
                              '2' = "LE5", '10' = "LE5", '1' = "LE6", '8' = "LE7", '6' = "LE8", '7' = "UrLE", '9' = "OE")
HGF_MET1.epi1[["EpiCellTypes"]] <- Idents(object = HGF_MET1.epi1)

#SFig.2i
Idents(object = HGF_MET1.epi1) <- "EpiCellTypes"
tiff(file = "HGF_MET1.epi1 EpiCellTypes 0.5 UMAP.tiff", width = 6.5, height = 6, units = "in", compression = "lzw", res = 800)
DimPlot(HGF_MET1.epi1, reduction = "umap", pt.size = 0.5, cols = c("red", "#FF9933", "#1D762E", "purple", "#FFD966", "#28CC44", "#3399FF", "#3333FF", "#51CCC9","#FF5CFF", "grey"))
dev.off()

#SFig.2j
Idents(object = HGF_MET1.epi1) <- "EpiCellTypes"
table(Idents(HGF_MET1.epi1))

#DEGs_allclusters
DefaultAssay(HGF_MET1.epi1) <- "RNA"
Idents(object = HGF_MET1.epi1) <- "EpiCellTypes"
HGF_MET1.epi1 <- ScaleData(HGF_MET1.epi1, features = rownames(HGF_MET1.epi1))
HGF_MET1.epi1.allMarkers <- FindAllMarkers(HGF_MET1.epi1, min.pct = 0.25, logfc.threshold = 0.25, only.pos = TRUE)
write.csv(HGF_MET1.epi1.allMarkers, "HGF_MET1.epi1.allMarkers.csv")

#SFig.2k
Idents(object = HGF_MET1.epi1) <- "EpiCellTypes"
tiff(file = "HGF_MET1.epi1 EpiCellType markers DotPlot.tiff", width =12 , height = 4, units = "in", compression = "lzw", res = 800)
DotPlot(HGF_MET1.epi1, features = c("hMETtg", "hHGFtg", "Krt14", "Krt17", "Lgals7", "Krt5", "Col17a1", 
                                    "Gm42418", "Gm26917", "Lars2", "Atf3", "Slc7a5", 
                                    "Rps26", "Serf2", "Sec61g", "Cox8a", "Pfdn5",
                                    "Tac1", "Ank", "Kcnk3", "Crabp1", "Spink5",
                                    "Laptm5", "Tnfrsf9", "Ptprc", "Coro1a", "Rac2",
                                    "Msmb", "Chodl", "Cmbl", "Serpinb11", "Reg3g",
                                    "Areg", "Btc", "Ctse", "Podxl", "Nlrp10",
                                    "Cmpk2", "Apol9a", "Mx1", "Gbp2", "Igtp", 
                                    "Fgl1", "Chn2", "Glb1l3", "Gsdma", "Pgm2l1",
                                    "Gsdmc2", "Gsdmc3", "Naip5", "Ces1d", "Cxcl15", 
                                    "Igfbp6", "Col1a2", "Sparcl1", "Serping1", "Bgn"
),cols = c("light grey", "red")) + RotatedAxis()
dev.off()

#Fig.2f
DefaultAssay(HGF_MET1.epi1) <- "RNA"
tiff(file = "HGF_MET1.epi1 hHGFtg 3.0 expression plots size 0.5.tiff", width = 6.5, height = 6, units = "in", compression = "lzw", res = 800)
FeaturePlot(HGF_MET1.epi1, reduction = "umap", features = c("hHGFtg"), cols = c("light grey", "red"), pt.size = 0.5, max.cutoff = 3.0)
dev.off()
tiff(file = "HGF_MET1.epi1 hMETtg 3.0 expression plots 0.5.tiff", width = 6.5, height = 6, units = "in", compression = "lzw", res = 800)
FeaturePlot(HGF_MET1.epi1, reduction = "umap", features = c("hMETtg"), cols = c("light grey", "red"), pt.size = 0.5, max.cutoff = 3.0)
dev.off()
tiff(file = "HGF_MET1.epi1 Pbsn 3.0 expression plots 0.5.tiff", width = 6.5, height = 6, units = "in", compression = "lzw", res = 800)
FeaturePlot(HGF_MET1.epi1, reduction = "umap", features = c("Pbsn"), cols = c("light grey", "red"), pt.size = 0.5, max.cutoff = 3.0)
dev.off()
tiff(file = "HGF_MET1.epi1 Ar 3.0 expression plots 0.5.tiff", width = 6.5, height = 6, units = "in", compression = "lzw", res = 800)
FeaturePlot(HGF_MET1.epi1, reduction = "umap", features = c("Ar"), cols = c("light grey", "red"), pt.size = 0.5, max.cutoff = 3.0)
dev.off()
tiff(file = "HGF_MET1.epi1 Krt5 3.0 expression plots 0.5.tiff", width = 6.5, height = 6, units = "in", compression = "lzw", res = 800)
FeaturePlot(HGF_MET1.epi1, reduction = "umap", features = c("Krt5"), cols = c("light grey", "red"), pt.size = 0.5, max.cutoff = 3.0)
dev.off()
tiff(file = "HGF_MET1.epi1 Krt8 3.0 expression plots 0.5.tiff", width = 6.5, height = 6, units = "in", compression = "lzw", res = 800)
FeaturePlot(HGF_MET1.epi1, reduction = "umap", features = c("Krt8"), cols = c("light grey", "red"), pt.size = 0.5, max.cutoff = 3.0)
dev.off()

####hMETtg+ vs hMETtg-####

#Add hMETtg info
DefaultAssay(HGF_MET1.epi1) <- "RNA"
HGF_MET1.epi1.hMETtgPos <- subset(x=HGF_MET1.epi1, subset = hMETtg > 0)
HGF_MET1.epi1.hMETtgNeg <- subset(x=HGF_MET1.epi1, subset = hMETtg == 0)
Idents(object = HGF_MET1.epi1.hMETtgPos) <- "hMETtgPos"
Idents(object = HGF_MET1.epi1.hMETtgNeg) <- "hMETtgNeg"
HGF_MET1.epi1.hMETtgPos[["hMETtgExp"]] <- Idents(object = HGF_MET1.epi1.hMETtgPos)
HGF_MET1.epi1.hMETtgNeg[["hMETtgExp"]] <- Idents(object = HGF_MET1.epi1.hMETtgNeg)
HGF_MET1.epi1.hMETtg <- merge(x = HGF_MET1.epi1.hMETtgPos, y = HGF_MET1.epi1.hMETtgNeg)
Idents(object = HGF_MET1.epi1.hMETtg) <- "hMETtgExp"
HGF_MET1.epi1$hMETtgExp <- Idents(object = HGF_MET1.epi1.hMETtg)

#subset LEs
Idents(object = HGF_MET1.epi1) <- "EpiCellTypes"
HGF_MET1.LE <- subset(HGF_MET1.epi1, idents = c("LE1", "LE2", "LE3", "LE4", "LE5", "LE6", "LE7", "LE8"))

#DEGs
DefaultAssay(HGF_MET1.LE) <- "RNA"
Idents(object = HGF_MET1.LE) <- "hMETtgExp"
HGF_MET1.LE <- ScaleData(HGF_MET1.LE, features = rownames(HGF_MET1.LE))
HGF_MET1.LE.0.Markers <- FindMarkers(HGF_MET1.LE, ident.1 = "hMETtgPos", ident.2 = "hMETtgNeg", min.pct = 0, logfc.threshold = 0)
write.csv(HGF_MET1.LE.0.Markers, "HGF_MET1.LE.0.Markers.csv")

#p.adjust
DEG_hMETtgPosvhMETtgNeg <- read.csv("HGF_MET1.LE.0.Markers.csv") 
DEG_hMETtgPosvhMETtgNeg_pvalue <- DEG_hMETtgPosvhMETtgNeg$p_val
DEG_hMETtgPosvhMETtgNeg_pvalue=as.numeric(DEG_hMETtgPosvhMETtgNeg_pvalue)
DEG_hMETtgPosvhMETtgNeg_BH = p.adjust(DEG_hMETtgPosvhMETtgNeg_pvalue, "BH")
write.csv(DEG_hMETtgPosvhMETtgNeg_BH, "DEG_hMETtgPosvhMETtgNeg_BH-1.csv")

#Fig.2g
DefaultAssay(HGF_MET1.LE) <- "RNA"
Idents(object = HGF_MET1.LE) <- "hMETtgExp"
HGF_MET1.LE <- ScaleData(HGF_MET1.LE, features = rownames(HGF_MET1.LE))
HGF_MET1.LE.all.markers <- FindAllMarkers(HGF_MET1.LE, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
HGF_MET1.LE.all.markers.Top50 <- HGF_MET1.LE.all.markers %>% group_by(cluster) %>% top_n(n = 50, wt = avg_log2FC)
tiff(file = "HGF_MET1.LE Heatmap Top50 purple.tiff", width = 6, height = 10, units = "in", compression = "lzw", res = 200)
DoHeatmap(HGF_MET1.LE, features = c(HGF_MET1.LE.all.markers.Top50$gene), draw.lines = TRUE) + scale_fill_gradientn(colors = c("purple", "black", "yellow"))
dev.off()

#Fig.2i
HGF_MET1.LE <- RenameIdents(object = HGF_MET1.LE, 'hMETtgNeg' = "hMETtgNeg", 'hMETtgPos' = "hMETtgPos")
HGF_MET1.LE[["hMETtgExp"]] <- Idents(object = HGF_MET1.LE)
Idents(object = HGF_MET1.LE) <- "hMETtgExp"
tiff(file = "HGF_MET1.LE hMETtg Vln.tiff", width = 5, height = 4, units = "in", compression = "lzw", res = 800)
VlnPlot(HGF_MET1.LE, features = "hMETtg", pt.size = 0.3, cols = c("#3399FF", "#E06666"))
dev.off()
tiff(file = "HGF_MET1.LE Plaur Vln.tiff", width = 5, height = 4, units = "in", compression = "lzw", res = 800)
VlnPlot(HGF_MET1.LE, features = "Plaur", pt.size = 0.3, cols = c("#3399FF", "#E06666"))
dev.off()
tiff(file = "HGF_MET1.LE Sox9 Vln.tiff", width = 5, height = 4, units = "in", compression = "lzw", res = 800)
VlnPlot(HGF_MET1.LE, features = "Sox9", pt.size = 0.3, cols = c("#3399FF", "#E06666"))
dev.off()
tiff(file = "HGF_MET1.LE Mmp7 Vln.tiff", width = 5, height = 4, units = "in", compression = "lzw", res = 800)
VlnPlot(HGF_MET1.LE, features = "Mmp7", pt.size = 0.3, cols = c("#3399FF", "#E06666"))
dev.off()
tiff(file = "HGF_MET1.LE Cd44 Vln.tiff", width = 5, height = 4, units = "in", compression = "lzw", res = 800)
VlnPlot(HGF_MET1.LE, features = "Cd44", pt.size = 0.3, cols = c("#3399FF", "#E06666"))
dev.off()
tiff(file = "HGF_MET1.LE Tcf7l2 Vln.tiff", width = 5, height = 4, units = "in", compression = "lzw", res = 800)
VlnPlot(HGF_MET1.LE, features = "Tcf7l2", pt.size = 0.3, cols = c("#3399FF", "#E06666"))
dev.off()

#Gene-Gene Spearman correlation
Idents(object = HGF_MET1.LE) <- "hMETtgExp"
GOI <- c('hMETtg', 'Plaur','Sox9', 'Mmp7', 'Cd44', 'Tcf7l2')  
GOI_index <- is.element(rownames(HGF_MET1.LE),GOI)
Cell_index <- is.element(Idents(HGF_MET1.LE), c('hMETtgPos','hMETtgNeg'))
expr_GOI <- HGF_MET1.LE@assays$RNA@data[GOI_index,Cell_index] 
expr_GOI <- HGF_MET1.LE@assays$RNA@counts[GOI_index,Cell_index]
ggcorr(t(expr_GOI), method = c("pairwise", "spearman"), label = TRUE)

#Fig.2j
tiff(file = "HGF_MET1.LE spearman correlation NoLabel.tiff", width = 5, height = 5, units = "in", compression = "lzw", res = 800)
ggcorr(t(expr_GOI), method = c("pairwise", "spearman"))
dev.off()

####Load data####
setwd("//isi-dcnl/user_data/zjsun/group/Lab Members/Won Kyung Kim/hHGFtg-hMETtg-b-cat-Pb/LabelledMulti")
list.files()
scRNA <- readRDS("scRNA.rds")

#subset primary
Idents(object = HGFMETBcat_Multi1) <- "orig.ident"
DimPlot(HGFMETBcat_Multi1, reduction = "umap")
HGFMETBcat_Primary <- subset(HGFMETBcat_Multi1, idents = c("Primary"))
HGFMETBcat_Primary <- FindVariableFeatures(HGFMETBcat_Primary, selection.method = "vst", nfeatures = 5000)

#Clustering
HGFMETBcat_Primary <- ScaleData(HGFMETBcat_Primary, verbose = FALSE)
HGFMETBcat_Primary <- RunPCA(HGFMETBcat_Primary, npcs = 50, verbose = FALSE)
ElbowPlot(HGFMETBcat_Primary, ndims = 50)
HGFMETBcat_Primary <- FindNeighbors(HGFMETBcat_Primary, reduction = "pca", dims = 1:20)
HGFMETBcat_Primary <- FindClusters(HGFMETBcat_Primary, resolution = 0.5)
HGFMETBcat_Primary <- RunTSNE(HGFMETBcat_Primary, reduction = "pca", dims = 1:20)
HGFMETBcat_Primary <- RunUMAP(HGFMETBcat_Primary, reduction = "pca", dims = 1:20)

#Fig.4a
Idents(object = HGF_MET1) <- "stim"
tiff(file = "HGF_MET1 gray UMAP.tiff", width = 6, height = 6, units = "in", compression = "lzw", res = 800)
DimPlot(HGF_MET1, reduction = "umap", pt.size = 0.3, cols = c("gray")) + NoLegend()
dev.off()
tiff(file = "HGFMETBcat_Primary darkblue UMAP.tiff", width = 6, height = 6, units = "in", compression = "lzw", res = 800)
DimPlot(HGFMETBcat_Primary, reduction = "umap", pt.size = 0.3, group.by = "stim", cols = c("darkblue")) + NoLegend()
dev.off()

####Merge Data####
#Stash old idents
HGFMETBcat_Primary[["orig.clusters"]] <- Idents(object = HGFMETBcat_Primary)
HGF_MET1[["orig.clusters"]] <- Idents(object = HGF_MET)

#Set Current idents
Idents(object = HGFMETBcat_Primary) <- "seurat_clusters"
Idents(object = HGF_MET1) <- "seurat_clusters"
HGFMETBcat_Primary$stim <- "Triple"
HGF_MET1$stim <- "Double"
TriplevDouble.anchors <- FindIntegrationAnchors(object.list = list(HGFMETBcat_Primary, HGF_MET1), dims = 1:20)
TriplevDouble.combined <- IntegrateData(anchorset = TriplevDouble.anchors, dims = 1:20)
DefaultAssay(TriplevDouble.combined) <- "integrated"

#clustering
TriplevDouble.combined <- ScaleData(TriplevDouble.combined, verbose = FALSE)
TriplevDouble.combined <- RunPCA(TriplevDouble.combined, npcs = 30, verbose = FALSE)
ElbowPlot(TriplevDouble.combined, ndims = 50)

#Umap and Clustering
TriplevDouble.combined <- FindNeighbors(TriplevDouble.combined, reduction = "pca", dims = 1:18)
TriplevDouble.combined <- FindClusters(TriplevDouble.combined, resolution = 0.5)
TriplevDouble.combined <- RunUMAP(TriplevDouble.combined, reduction = "pca", dims = 1:18)
Idents(object = TriplevDouble.combined) <- "stim"
TriplevDouble.combined <- RenameIdents(object = TriplevDouble.combined, 'Double' = "Double", 'Triple' = "Triple")
TriplevDouble.combined[["stim"]] <- Idents(object = TriplevDouble.combined)
DimPlot(TriplevDouble.combined, reduction = "umap", pt.size = 0.3, label = TRUE)

#Cell cycle assignment
mouse_cell_cycle_genes <- readRDS(file = "//isi-dcnl/user_data/zjsun/group/Won Kyung Kim/ARKOvCtrl_3timepoint/E18.5/mouse_cell_cycle_genes.rds")
s.genes <- mouse_cell_cycle_genes$s.genes
g2m.genes <- mouse_cell_cycle_genes$g2m.genes
DefaultAssay(TriplevDouble.combined) <- "RNA"
all.genes <- rownames(TriplevDouble.combined)
TriplevDouble.combined <- ScaleData(TriplevDouble.combined, features = all.genes)
TriplevDouble.combined <- CellCycleScoring(TriplevDouble.combined, s.features = s.genes, g2m.features = g2m.genes, set.ident = TRUE)

#Cell Cycle regression
TriplevDouble.combined1 <- TriplevDouble.combined
DefaultAssay(TriplevDouble.combined1) <- "integrated"
TriplevDouble.combined1 <- ScaleData(TriplevDouble.combined1, vars.to.regress = c("S.Score", "G2M.Score"), features = rownames(TriplevDouble.combined1))
TriplevDouble.combined1 <- RunPCA(TriplevDouble.combined1, features = VariableFeatures(TriplevDouble.combined1))
ElbowPlot(TriplevDouble.combined1, ndims = 50)
TriplevDouble.combined1 <- FindNeighbors(TriplevDouble.combined1, reduction = "pca", dims = 1:26)
TriplevDouble.combined1 <- FindClusters(TriplevDouble.combined1, resolution = 1.5)
TriplevDouble.combined1 <- RunUMAP(TriplevDouble.combined1, reduction = "pca", dims = 1:26)

#SFig.4a
Idents(object = TriplevDouble.combined) <- "Phase"
tiff(file = "TriplevDouble.combined Cell Cyle UMAP.tiff", width = 6.5, height = 6, units = "in", compression = "lzw", res = 800)
DimPlot(TriplevDouble.combined, reduction = "umap", pt.size = 0.05, cols = c("goldenrod2", "grey75", "deepskyblue2"))
dev.off()
Idents(object = TriplevDouble.combined1) <- "Phase"
tiff(file = "TriplevDouble.combined1 Cell Cyle UMAP.tiff", width = 6.5, height = 6, units = "in", compression = "lzw", res = 800)
DimPlot(TriplevDouble.combined1, reduction = "umap", pt.size = 0.05, cols = c("goldenrod2", "grey75", "deepskyblue2"))
dev.off()

#Cell type identification
DefaultAssay(TriplevDouble.combined1)<-"RNA"
FeaturePlot(TriplevDouble.combined1, reduction = "umap", features = c("Krt5", "Krt4", "Krt8", "Cd24a", "Plp1",
                                                                      "Fbln1", "Myh11", "Pecam1",
                                                                      "Vim", "Tyrobp", "Ccl5", "Pate4"), cols = c("light grey", "red"), pt.size = 0.3, max.cutoff = "q90")

#Rename CellTypes
Idents(object = TriplevDouble.combined1) <- "seurat_clusters"
DimPlot(TriplevDouble.combined1, reduction = "umap", pt.size = 0.3, label=TRUE)
TriplevDouble.combined1 <- RenameIdents(object = TriplevDouble.combined1, 
                                        '1'="BE",'3'="BE",'20'="BE",
                                        '6'="LE", '7'="LE", "0"="LE", '19'="LE",
                                        '2'="LE",'12'="LE",'5'="LE", '8'="LE",'11'="LE",
                                        '22'="LE",'21'="LE",'18'="LE",'13'="LE",'4'="LE",'25'="SV",
                                        '26'="SV", '16'="SV",
                                        '10'="FB",'14'="FB", '23'="SM",'27'="Pericyte", '17'="VE",
                                        '9'="Immune",'24'="Immune", '15'="Immune")  
TriplevDouble.combined1[["CellTypes"]] <- Idents(object = TriplevDouble.combined1)

#SFig.4b
Idents(object = TriplevDouble.combined1) <- "CellTypes"
tiff(file = "TriplevDouble.combined1 CellTypes UMAP.tiff", width = 6.5, height = 6, units = "in", compression = "lzw", res = 800)
DimPlot(TriplevDouble.combined1, reduction = "umap", pt.size = 0.3, cols = c("chartreuse3", "salmon", "darkslategray3", "plum4", "brown3", "blueviolet", "blue", "bisque3", "steelblue1"))
dev.off()

#SFig.4c
DefaultAssay(TriplevDouble.combined1) <- "RNA"
tiff(file = "TriplevDouble.combined1 hMETtg split expression plots.tiff", width = 10, height = 5, units = "in", compression = "lzw", res = 800)
FeaturePlot(TriplevDouble.combined1, reduction = "umap", features = c("hMETtg"), split.by = "stim", cols = c("light grey", "red"), pt.size = 0.3, max.cutoff = "q90", keep.scale = "all")
dev.off()
tiff(file = "TriplevDouble.combined1 Ar split expression plots.tiff", width = 10, height = 5, units = "in", compression = "lzw", res = 800)
FeaturePlot(TriplevDouble.combined1, reduction = "umap", features = c("Ar"), split.by = "stim", cols = c("light grey", "red"), pt.size = 0.3, max.cutoff = "q90", keep.scale = "all")
dev.off()
tiff(file = "TriplevDouble.combined1 Pbsn split expression plots.tiff", width = 10, height = 5, units = "in", compression = "lzw", res = 800)
FeaturePlot(TriplevDouble.combined1, reduction = "umap", features = c("Pbsn"), split.by = "stim", cols = c("light grey", "red"), pt.size = 0.3, max.cutoff = "q90", keep.scale = "all")
dev.off()
tiff(file = "TriplevDouble.combined1 Krt5 split expression plots.tiff", width = 10, height = 5, units = "in", compression = "lzw", res = 800)
FeaturePlot(TriplevDouble.combined1, reduction = "umap", features = c("Krt5"), split.by = "stim", cols = c("light grey", "red"), pt.size = 0.3, max.cutoff = "q90", keep.scale = "all")
dev.off()
tiff(file = "TriplevDouble.combined1 Krt8 split expression plots.tiff", width = 10, height = 5, units = "in", compression = "lzw", res = 800)
FeaturePlot(TriplevDouble.combined1, reduction = "umap", features = c("Krt8"), split.by = "stim", cols = c("light grey", "red"), pt.size = 0.3, min.cutoff = "q5", max.cutoff = "q90", keep.scale = "all")
dev.off()

#DEGs
DefaultAssay(TriplevDouble.combined1) <- "RNA"
Idents(object = TriplevDouble.combined1) <- "CellTypes"
TriplevDouble.combined1 <- ScaleData(TriplevDouble.combined1, features = rownames(TriplevDouble.combined1))
TriplevDouble.combined1.allMarkers <- FindAllMarkers(TriplevDouble.combined1, min.pct = 0.25, logfc.threshold = 0.25, only.pos = TRUE)
write.csv(TriplevDouble.combined1.allMarkers, "TriplevDouble.combined1.CellTypes.Markers.csv")

#SFig.4d
Idents(object = TriplevDouble.combined1) <- "CellTypes"
tiff(file = "TriplevDouble.combined1 CellTypes markers DotPlot.tiff", width =12 , height = 3.1, units = "in", compression = "lzw", res = 800)
DotPlot(TriplevDouble.combined1, features = c("Krt15", "Krt14", "Krt5", "Aqp3", "Col17a1", 
                                              "Prr9", "Pigr","Slc12a2", "Arl14", "Bex1",
                                              "Svs4",  "Pate4", "A630095E13Rik", "D730048I06Rik", "Sprr2f", 
                                              "Apod", "Bgn", "Igfbp6", "Serping1", "Col1a2",
                                              "Ndufa4l2", "Adamts4", "Myh11", "Pdgfrb", "Des",
                                              "Plp1", "Kcna1", "Cdh19", "S100b", "Fxyd1",
                                              "Aqp1", "Plvap", "Cdh5", "Pecam1", "Cd93",
                                              "H2-Eb1", "C1qa", "Rgs1", "Tyrobp", "Fcer1g"
),cols = c("light grey", "red")) + RotatedAxis()
dev.off()

#SFig.4g
Idents(object = TriplevDouble.combined1) <- "CellTypes"
TriplevDouble.combined1$stim.CellTypes <- paste(Idents(TriplevDouble.combined1), TriplevDouble.combined1$stim, sep = "_")
Idents(object = TriplevDouble.combined1) <- "stim.CellTypes"
table(Idents(TriplevDouble.combined1))

####Subcluster epi TriplevDouble####
Idents(object = TriplevDouble.combined1) <- "CellTypes"
TriplevDouble.combined1.epi <- subset(TriplevDouble.combined1, idents = c("BE","LE"))
DimPlot(TriplevDouble.combined1.epi, reduction = "umap", pt.size = 0.3, label = TRUE)

#Run the standard workflow for visualization and clustering
DefaultAssay(TriplevDouble.combined1.epi) <- "integrated"
TriplevDouble.combined1.epi <- ScaleData(TriplevDouble.combined1.epi, verbose = FALSE)
TriplevDouble.combined1.epi <- RunPCA(TriplevDouble.combined1.epi, npcs = 50, verbose = FALSE)
ElbowPlot(TriplevDouble.combined1.epi, ndims = 50)

#Umap and Clustering
TriplevDouble.combined1.epi <- FindNeighbors(TriplevDouble.combined1.epi, reduction = "pca", dims = 1:20)
TriplevDouble.combined1.epi <- FindClusters(TriplevDouble.combined1.epi, resolution = 0.7)
TriplevDouble.combined1.epi <- RunUMAP(TriplevDouble.combined1.epi, reduction = "pca", dims = 1:20)
DimPlot(TriplevDouble.combined1.epi, reduction = "umap", pt.size = 0.3, label = TRUE)

#Cell cycle assignment epi
DefaultAssay(TriplevDouble.combined1.epi) <- "RNA"
all.genes <- rownames(TriplevDouble.combined1.epi)
TriplevDouble.combined1.epi <- ScaleData(TriplevDouble.combined1.epi, features = all.genes)
TriplevDouble.combined1.epi <- CellCycleScoring(TriplevDouble.combined1.epi, s.features = s.genes, g2m.features = g2m.genes, set.ident = TRUE)

#Cell Cycle Regression epi
TriplevDouble.combined1.epi2 <- TriplevDouble.combined1.epi
DefaultAssay(TriplevDouble.combined1.epi2) <- "integrated"
TriplevDouble.combined1.epi2 <- ScaleData(TriplevDouble.combined1.epi2, vars.to.regress = c("S.Score", "G2M.Score"), features = rownames(TriplevDouble.combined1.epi2))
TriplevDouble.combined1.epi2 <- RunPCA(TriplevDouble.combined1.epi2, features = VariableFeatures(TriplevDouble.combined1.epi2))
ElbowPlot(TriplevDouble.combined1.epi2, ndims = 50)
TriplevDouble.combined1.epi2 <- FindNeighbors(TriplevDouble.combined1.epi2, reduction = "pca", dims = 1:22)
TriplevDouble.combined1.epi2 <- FindClusters(TriplevDouble.combined1.epi2, resolution = 2.5)
TriplevDouble.combined1.epi2 <- RunUMAP(TriplevDouble.combined1.epi2, reduction = "pca", dims = 1:22)
DimPlot(TriplevDouble.combined1.epi2, reduction = "umap", pt.size = 0.3, label = TRUE, split.by = "stim")

#SFig.4e
Idents(object = TriplevDouble.combined1.epi) <- "Phase"
tiff(file = "TriplevDouble.combined1.epi Cell Cyle UMAP.tiff", width = 6.5, height = 6, units = "in", compression = "lzw", res = 800)
DimPlot(TriplevDouble.combined1.epi, reduction = "umap", pt.size = 0.05, cols = c("goldenrod2", "grey75", "deepskyblue2"))
dev.off()
tiff(file = "TriplevDouble.combined1.epi2 Cell Cyle UMAP.tiff", width = 6.5, height = 6, units = "in", compression = "lzw", res = 800)
DimPlot(TriplevDouble.combined1.epi2, reduction = "umap", pt.size = 0.05, cols = c("goldenrod2", "grey75", "deepskyblue2"))
dev.off()

#DEGs
DefaultAssay(TriplevDouble.combined1.epi2) <- "RNA"
Idents(object = TriplevDouble.combined1.epi2) <- "seurat_clusters"
TriplevDouble.combined1.epi2 <- ScaleData(TriplevDouble.combined1.epi2, features = rownames(TriplevDouble.combined1.epi2))
TriplevDouble.combined1.epi2.seurat.Markers <- FindAllMarkers(TriplevDouble.combined1.epi2, min.pct = 0.25, logfc.threshold = 0.25, only.pos = TRUE)
write.csv(TriplevDouble.combined1.epi2.seurat.Markers, "TriplevDouble.combined1.epi2.seurat.Markers.csv")

#Rename Clusters
Idents(object = TriplevDouble.combined1.epi2) <- "seurat_clusters"
TriplevDouble.combined1.epi2 <- RenameIdents(object = TriplevDouble.combined1.epi2, 
                                             '33'="BE1", 
                                             '6'="BE2",'16'="BE2",'29'="BE2",
                                             '15'="BE3", '5'="BE3", 
                                             '23' = "BE4", '30'="BE4", '28' = "BE4", '34'="BE4", 
                                             '35'="LE1",'21'="LE1", '12'="LE1", '13'="LE1",
                                             '25'="LE2", '10'="LE2",'17'="LE2",
                                             '18' = "LE3", '14' ="LE3",'32' = "LE3", '31' ="LE3",'4' ="LE3", '20'="LE3",
                                             '11' = "LE4", '1'="LE4",
                                             '0'="LE5",'9'="LE5",
                                             '2'="LE6",'36'="LE6",
                                             '3'="LE7", 
                                             '22'="LE8",'8'="LE8", '19'="LE8",
                                             '7'="LE9", '24'="LE9",
                                             '26'="UrLE", 
                                             '27'="OE")  
TriplevDouble.combined1.epi2[["EpiCellTypes"]] <- Idents(object = TriplevDouble.combined1.epi2)

#SFig.4f
tiff(file = "TriplevDouble.combined1.epi2 EpiCellTypes UMAP.tiff", width = 5.5, height = 5, units = "in", compression = "lzw", res = 800)
DimPlot(TriplevDouble.combined1.epi2, reduction = "umap", pt.size = 0.3, label.size = 6, 
        cols = c("salmon", "skyblue1", "olivedrab2", "brown3", "deeppink1",  "blue", "darkorange", "red", "turquoise3",
                 "bisque3", "slategray3", "mediumorchid3", 
                 "yellow2", "green4", "black"))
dev.off()

#Fig.4b
tiff(file = "TriplevDouble.combined1.epi2 EpiCellTypes split UMAP.tiff", width = 10, height = 5, units = "in", compression = "lzw", res = 800)
DimPlot(TriplevDouble.combined1.epi2, reduction = "umap", pt.size = 0.3, label.size = 6, split.by = "stim",
        cols = c("salmon", "skyblue1", "olivedrab2", "brown3", "deeppink1",  "blue", "darkorange", "red", "turquoise3",
                 "bisque3", "slategray3", "mediumorchid3", 
                 "yellow2", "green4", "black"))
dev.off()

#Fig.4b
DefaultAssay(TriplevDouble.combined1.epi2)<-"RNA"
Idents(object = TriplevDouble.combined1.epi2) <- "EpiCellTypes"
tiff(file = "TriplevDouble.combined1.epi2 Ar split plots.tiff", width = 10, height = 5, units = "in", compression = "lzw", res = 800)
FeaturePlot(TriplevDouble.combined1.epi2, reduction = "umap", features = c("Ar"), split.by = "stim", cols = c("light grey", "red"), pt.size = 0.3, max.cutoff = "q95")
dev.off()
tiff(file = "TriplevDouble.combined1.epi2 Pbsn split plots.tiff", width = 10, height = 5, units = "in", compression = "lzw", res = 800)
FeaturePlot(TriplevDouble.combined1.epi2, reduction = "umap", features = c("Pbsn"), split.by = "stim", cols = c("light grey", "red"), pt.size = 0.3, max.cutoff = "q95")
dev.off()
tiff(file = "TriplevDouble.combined1.epi2 Fkbp5 plots.tiff", width = 10, height = 5, units = "in", compression = "lzw", res = 800)
FeaturePlot(TriplevDouble.combined1.epi2, reduction = "umap", features = c("Fkbp5"), split.by = "stim", cols = c("light grey", "red"), pt.size = 0.3, max.cutoff = "q95", keep.scale = "all")
dev.off()
tiff(file = "TriplevDouble.combined1.epi2 Tcf4 split plots.tiff", width = 10, height = 5, units = "in", compression = "lzw", res = 800)
FeaturePlot(TriplevDouble.combined1.epi2, reduction = "umap", features = c("Tcf4"), split.by = "stim", cols = c("light grey", "red"), pt.size = 0.3, max.cutoff = "q90", keep.scale = "all")
dev.off()
tiff(file = "TriplevDouble.combined1.epi2 Axin2 split plots.tiff", width = 10, height = 5, units = "in", compression = "lzw", res = 800)
FeaturePlot(TriplevDouble.combined1.epi2, reduction = "umap", features = c("Axin2"), split.by = "stim", cols = c("light grey", "red"), pt.size = 0.3, max.cutoff = "q90", keep.scale = "all")
dev.off()

#SFig.4h
Idents(object = TriplevDouble.combined1.epi2) <- "EpiCellTypes"
tiff(file = "TriplevDouble.combined1.epi2 EpiCellTypes markers DotPlot.tiff", width =17.5 , height = 4.8, units = "in", compression = "lzw", res = 800)
DotPlot(TriplevDouble.combined1.epi2, features = c("Tmem171", "Egfl6", "Ncam1", "Clca3a2", "Adm", 
                                                   "Krt15", "Palld", "Ctsl", "Tubb6", "Tpm1",
                                                   "Aqp3", "Lgals7", "Col17a1", "Lamb3", "Pvrl1",
                                                   "Sncg", "Ifi202b", "Gpnmb", "Dapl1", "Gpr87", 
                                                   "Lars2", "AY036118", "Gm42418", "Gm26917", "Hbb-bs",
                                                   "Defa21", "Gm15293", "Defa5", "Cryba4", "Ltf", 
                                                   "Rnf149", "Car2", "Lap3", "Plat", "Tm4sf1", 
                                                   "Coch", "Tgfb2", "Dkk2", "Zeb2", "Apoc1",
                                                   "Crip1", "Tspan8", "Ly6c1", "Btc", "B2m", 
                                                   "Tgm4", "9530053A07Rik", "Gm5615", "Man1a", "Spink8", 
                                                   "Msmb", "Mme", "Apof", "Pcp4", "Agtr1a",
                                                   "Pigr", "Tspan1", "Tnfrsf21", "Dcxr", "Cldn3",
                                                   "Spink1", "Sbpl", "Crabp1", "Col6a3", 'Gucy2g', 
                                                   "Gsdmc2", "Gsdmc3", "Barx2", "Cxcl15", "Krt4", 
                                                   "Serping1", "Igfbp6", "Fbln1", "Serpinf1", "Col1a2"
),cols = c("light grey", "red")) + RotatedAxis()
dev.off()

#SFig.4j
Idents(object = TriplevDouble.combined1.epi2) <- "EpiCellTypes"
TriplevDouble.combined1.epi2$stim.EpiCellTypes <- paste(Idents(TriplevDouble.combined1.epi2), TriplevDouble.combined1.epi2$stim, sep = "_")
Idents(object = TriplevDouble.combined1.epi2) <- "stim.EpiCellTypes"
table(Idents(TriplevDouble.combined1.epi2))

####hMETtgPos vs hMETtgNeg####
#Add hMETtg information
DefaultAssay(TriplevDouble.combined1.epi2) <- "RNA"
TriplevDouble.combined1.epi2METPos <- subset(x=TriplevDouble.combined1.epi2,  subset = `hMETtg` > 0)
TriplevDouble.combined1.epi2METNeg <- subset(x=TriplevDouble.combined1.epi2,  subset = `hMETtg` == 0)
Idents(object = TriplevDouble.combined1.epi2METPos) <- "METPos"
Idents(object = TriplevDouble.combined1.epi2METNeg) <- "METNeg"
TriplevDouble.combined1.epi2METPos[["METExp"]] <- Idents(object = TriplevDouble.combined1.epi2METPos)
TriplevDouble.combined1.epi2METNeg[["METExp"]] <- Idents(object = TriplevDouble.combined1.epi2METNeg)
TriplevDouble.combined1.epi2MET <- merge(x = TriplevDouble.combined1.epi2METPos, y = TriplevDouble.combined1.epi2METNeg)
TriplevDouble.combined1.epi2$METExp <- Idents(object = TriplevDouble.combined1.epi2MET)
Idents(object = TriplevDouble.combined1.epi2) <- "METExp"

#Subset hMETPos BELE
Idents(object = TriplevDouble.combined1.epi2) <- "METExp"
TriplevDouble.combined1.epi2METPos <- subset(TriplevDouble.combined1.epi2, idents = c("METPos"))
Idents(object = TriplevDouble.combined1.epi2METPos) <- "EpiCellTypes"
TriplevDouble.combined1.BELEMETPos <- subset(TriplevDouble.combined1.epi2METPos, idents = c("BE1","BE2", "BE3", "BE4",
                                                                                            "LE1", "LE2", "LE3", "LE4",
                                                                                            "LE5", "LE6", "LE7", "LE8", "LE9"))

#Fig.4c & SFig.4i
DefaultAssay(TriplevDouble.combined1.BELEMETPos) <- "RNA"
Idents(object = TriplevDouble.combined1.BELEMETPos) <- "stim"
all.genes <- rownames(TriplevDouble.combined1.BELEMETPos)
TriplevDouble.combined1.BELEMETPos <- ScaleData(TriplevDouble.combined1.BELEMETPos, features = all.genes)
hMETtgPos_BELE_TriplevDouble.0.Markers <- FindMarkers(TriplevDouble.combined1.BELEMETPos, ident.1 = c("Triple"), 
                                                      ident.2 = c("Double"), min.pct = 0, logfc.threshold = 0)
write.csv(hMETtgPos_BELE_TriplevDouble.0.Markers, "hMETtgPos_BELE_TriplevDouble.0.Markers.csv")

#p.adjust
DEG_hMETtgPos_BELE_TriplevDouble <- read.csv("hMETtgPos_BELE_TriplevDouble.0.Markers.csv") 
DEG_hMETtgPos_BELE_TriplevDouble_pvalue <- DEG_hMETtgPos_BELE_TriplevDouble$p_val
DEG_hMETtgPos_BELE_TriplevDouble_pvalue=as.numeric(DEG_hMETtgPos_BELE_TriplevDouble_pvalue)
DEG_hMETtgPos_BELE_TriplevDouble_BH = p.adjust(DEG_hMETtgPos_BELE_TriplevDouble_pvalue, "BH")
write.csv(DEG_hMETtgPos_BELE_TriplevDouble_BH, "DEG_hMETtgPos_BELE_TriplevDouble_BH.csv")

#qval
library(qvalue)
SupplementaryData2 <- read_csv("Desktop/TripleTg Revision/SupplementaryData2.csv")
pvalues <- SupplementaryData2$pval
qobj <- qvalue_truncp(p = pvalues)
qvalues <- qobj$qvalues
write.csv(qvalues, file = "SupplementaryData2_final.csv")

#Fig.4d
data_summary <- function(x) {
  m <- mean(x)
  ymin <- m-sd(x)
  ymax <- m+sd(x)
  return(c(y=m,ymin=ymin,ymax=ymax))
}
DefaultAssay(TriplevDouble.combined1.BELEMETPos) <- "RNA"
Idents(object = TriplevDouble.combined1.BELEMETPos) <- "stim"
tiff(file = "TriplevDouble.combined1.BELEMETPos Ar Vln.tiff", width = 3.7, height = 4.5, units = "in", compression = "lzw", res = 800)
VlnPlot(TriplevDouble.combined1.BELEMETPos, features = "Ar", pt.size = 0, cols = c("#3399FF",   "#E06666")) + NoLegend() +
  stat_summary(fun.data=data_summary) + stat_summary(fun = mean, geom='point', size = 25, colour = "green2", shape = 95)
dev.off()
tiff(file = "TriplevDouble.combined1.BELEMETPos Pbsn Vln.tiff", width = 3.7, height = 4.5, units = "in", compression = "lzw", res = 800)
VlnPlot(TriplevDouble.combined1.BELEMETPos, features = "Pbsn", pt.size = 0, cols = c("#3399FF",   "#E06666")) + NoLegend() +
  stat_summary(fun.data=data_summary) + stat_summary(fun = mean, geom='point', size = 25, colour = "green2", shape = 95)
dev.off()
tiff(file = "TriplevDouble.combined1.BELEMETPos Fkbp5 Vln.tiff", width = 3.7, height = 4.5, units = "in", compression = "lzw", res = 800)
VlnPlot(TriplevDouble.combined1.BELEMETPos, features = "Fkbp5", pt.size = 0, cols = c("#3399FF",   "#E06666")) + NoLegend() +
  stat_summary(fun.data=data_summary) + stat_summary(fun = mean, geom='point', size = 25, colour = "green2", shape = 95)
dev.off()
tiff(file = "TriplevDouble.combined1.BELEMETPos Axin2 Vln.tiff", width = 3.7, height = 4.5, units = "in", compression = "lzw", res = 800)
VlnPlot(TriplevDouble.combined1.BELEMETPos, features = "Axin2", pt.size = 0, cols = c("#3399FF",   "#E06666")) + NoLegend() +
  stat_summary(fun.data=data_summary) + stat_summary(fun = mean, geom='point', size = 25, colour = "green2", shape = 95)
dev.off()
tiff(file = "TriplevDouble.combined1.BELEMETPos Tcf4 Vln.tiff", width = 3.7, height = 4.5, units = "in", compression = "lzw", res = 800)
VlnPlot(TriplevDouble.combined1.BELEMETPos, features = "Tcf4", pt.size = 0, cols = c("#3399FF",   "#E06666")) + NoLegend() +
  stat_summary(fun.data=data_summary) + stat_summary(fun = mean, geom='point', size = 25, colour = "green2", shape = 95)
dev.off()
tiff(file = "TriplevDouble.combined1.BELEMETPos Xpo1 Vln.tiff", width = 3.7, height = 4.5, units = "in", compression = "lzw", res = 800)
VlnPlot(TriplevDouble.combined1.BELEMETPos, features = "Xpo1", pt.size = 0, cols = c("#3399FF",   "#E06666"), y.max = 1) + NoLegend() +
  stat_summary(fun.data=data_summary) + stat_summary(fun = mean, geom='point', size = 25, colour = "green2", shape = 95)
dev.off()
tiff(file = "TriplevDouble.combined1.BELEMETPos Rpl12 Vln.tiff", width = 3.7, height = 4.5, units = "in", compression = "lzw", res = 800)
VlnPlot(TriplevDouble.combined1.BELEMETPos, features = "Rpl12", pt.size = 0, cols = c("#3399FF",   "#E06666")) + NoLegend() +
  stat_summary(fun.data=data_summary) + stat_summary(fun = mean, geom='point', size = 25, colour = "green2", shape = 95)
dev.off()
tiff(file = "TriplevDouble.combined1.BELEMETPos Rps16 Vln.tiff", width = 3.7, height = 4.5, units = "in", compression = "lzw", res = 800)
VlnPlot(TriplevDouble.combined1.BELEMETPos, features = "Rps16", pt.size = 0, cols = c("#3399FF",   "#E06666")) + NoLegend() +
  stat_summary(fun.data=data_summary) + stat_summary(fun = mean, geom='point', size = 25, colour = "green2", shape = 95)
dev.off()

####Heatmap####
library(tidyverse)
library(pheatmap)
library(ggplot2)
library(ggplotify)
#scale.data
Epi_df <- as.data.frame(t(TriplevDouble.combined1.epi2$RNA@scale.data)) 
Epi_df_selected <- cbind( Epi_df$Pbsn, Epi_df$`Nkx3-1`, Epi_df$Fkbp5, Epi_df$Azgp1, 
                          Epi_df$Hpn, Epi_df$Mtor, Epi_df$Uba52, Epi_df$Mapk8, 
                          Epi_df$Axin2, Epi_df$Cd44, Epi_df$Dkk2, Epi_df$Tcf4, 
                          Epi_df$Rpl3, Epi_df$Rpl12,Epi_df$Rps5, Epi_df$Rps16,
                          Epi_df$Pcna, Epi_df$Cenpf, Epi_df$Pttg1, Epi_df$Cdk4 )
Epi_df_selected <- as.data.frame(Epi_df_selected)
colnames(Epi_df_selected) <- c( "Pbsn",	'Nkx3-1', "Fkbp5", "Azgp1",
                                "Hpn", "Mtor", "Uba52", "Mapk8", 
                                "Axin2", "Cd44",	"Dkk2", "Tcf4", 
                                "Rpl3", "Rpl12", "Rps5", "Rps16",
                                "Pcna",	"Cenpf", "Pttg1", "Cdk4")
rownames(Epi_df_selected) <- row.names(Epi_df)
write.csv(Epi_df_selected, file = "Epi_df_selected_scaledata.csv")

#meta.data
write.csv(TriplevDouble.combined1.epi2@meta.data, file = "TriplevDouble.combined1.epi2_metadata.csv")

#Env
df <- read.csv("Heatmap_Triple_and_Double_Epi_EpiCellTypes.csv", header = TRUE, sep = ",")
colnames(df)[1] <- "Clusters"
df <- as.data.frame(df)
df <- column_to_rownames(df, var = "Clusters")
df <- as.matrix(df)
d <- pheatmap::pheatmap(df, cluster_cols = F, cluster_rows = F)

#Fig.4f
tiff(file = "TriplevDouble.combined1.Epi EpiCellTypes Heatmap G2M AR MET WNT AR Ribosome.tiff", width = 6, height = 4, units = "in", compression = "lzw", res = 800)
pheatmap::pheatmap(log10(df+10),scale="column",
                   fontsize = 10,color = colorRampPalette(c("purple","grey0", "yellow"))(101),
                   cluster_cols = F, cluster_rows = F)
dev.off()

####hMETtgPos Monocle3####
Idents(object = TriplevDouble.combined1.BELEMETPos) <- "stim"
Double.combined1.BELEMETPos <- subset(TriplevDouble.combined1.BELEMETPos, idents = c("Double"))
Triple.combined1.BELEMETPos <- subset(TriplevDouble.combined1.BELEMETPos, idents = c("Triple"))

###Triple
##Triple.BELEMETPos re-clustering
DefaultAssay(Triple.combined1.BELEMETPos) <- "RNA"
Triple.combined1.BELEMETPos <- FindVariableFeatures(Triple.combined1.BELEMETPos, selection.method = "vst", nfeatures = 5000)
Triple.combined1.BELEMETPos <- ScaleData(Triple.combined1.BELEMETPos, verbose = FALSE)
Triple.combined1.BELEMETPos <- RunPCA(Triple.combined1.BELEMETPos, npcs = 50, verbose = FALSE)
ElbowPlot(Triple.combined1.BELEMETPos, ndims = 50)
Triple.combined1.BELEMETPos <- FindNeighbors(Triple.combined1.BELEMETPos, reduction = "pca", dims = 1:14)
Triple.combined1.BELEMETPos <- FindClusters(Triple.combined1.BELEMETPos, resolution = 0.5)
Triple.combined1.BELEMETPos <- RunUMAP(Triple.combined1.BELEMETPos, reduction = "pca", dims = 1:14)

#Cell cycle scoring
DefaultAssay(Triple.combined1.BELEMETPos) <- "RNA"
all.genes <- rownames(Triple.combined1.BELEMETPos)
Triple.combined1.BELEMETPos <- ScaleData(Triple.combined1.BELEMETPos, features = all.genes)
Triple.combined1.BELEMETPos <- CellCycleScoring(Triple.combined1.BELEMETPos, s.features = s.genes, g2m.features = g2m.genes, set.ident = TRUE)

#Cell Cycle regression
Triple.combined1.BELEMETPos1 <- Triple.combined1.BELEMETPos
DefaultAssay(Triple.combined1.BELEMETPos1) <- "integrated"
Triple.combined1.BELEMETPos1 <- ScaleData(Triple.combined1.BELEMETPos1, vars.to.regress = c("S.Score", "G2M.Score"), features = rownames(Triple.combined1.BELEMETPos1))
Triple.combined1.BELEMETPos1 <- RunPCA(Triple.combined1.BELEMETPos1, features = VariableFeatures(Triple.combined1.BELEMETPos1))
ElbowPlot(Triple.combined1.BELEMETPos1, ndims = 50)
Triple.combined1.BELEMETPos1 <- FindNeighbors(Triple.combined1.BELEMETPos1, reduction = "pca", dims = 1:16)
Triple.combined1.BELEMETPos1 <- FindClusters(Triple.combined1.BELEMETPos1, resolution = 0.5)
Triple.combined1.BELEMETPos1 <- RunUMAP(Triple.combined1.BELEMETPos1, reduction = "pca", dims = 1:16)

#Convert Seurat to Monocle3 cell data set class
DefaultAssay(Triple.combined1.BELEMETPos1) <- "RNA"
cds <- as.cell_data_set(Triple.combined1.BELEMETPos1)

#Calculate size factors using built-in function in monocle3
cds <- estimate_size_factors(cds)

#Add gene names into CDS
cds@rowRanges@elementMetadata@listData[["gene_short_name"]] <- rownames(Triple.combined1.BELEMETPos1)

#Construction single cell trajectories
cds <- cluster_cells(cds = cds, reduction_method = "UMAP") 

#plot single cell trajectory
cds <- learn_graph(cds, use_partition = FALSE)
plt <- plot_cells(cds,
                  color_cells_by = "EpiCellTypes",
                  show_trajectory_graph = TRUE,
                  label_branch_points = FALSE,
                  trajectory_graph_color = "black",
                  trajectory_graph_segment_size = 1.8,
                  label_cell_groups = FALSE,
                  label_leaves=FALSE,
                  label_roots = FALSE,
                  cell_size = 0.6,
                  alpha = 1) 

#Fig.4g
tiff(file = "Triple.combined1.BELEMETPos1 re-clustering trajectory UMAP.tiff", width = 6, height = 5, units = "in", compression = "lzw", res = 800)
plt + scale_color_manual(values = c("salmon", "skyblue1", "olivedrab2", "brown3", "deeppink1",  "blue", "darkorange", "red", "turquoise3",
                                    "bisque3", "slategray3", "mediumorchid3", 
                                    "yellow2", "green4", "black"))
dev.off()

#plot single cell trajectory by pseudotime
cds <- order_cells(cds, reduction_method = "UMAP")
plt <- plot_cells(cds,
                  color_cells_by = "pseudotime",
                  show_trajectory_graph = TRUE,
                  label_branch_points = FALSE,
                  trajectory_graph_color = "black",
                  trajectory_graph_segment_size = 1.8,
                  label_cell_groups = FALSE,
                  label_leaves=FALSE,
                  label_roots = FALSE,
                  cell_size = 0.6,
                  alpha = 1)

#Fig.4g
tiff(file = "Triple.combined1.BELEMETPos1 re-clustering pseudotime UMAP.tiff", width = 6, height = 5, units = "in", compression = "lzw", res = 800)
plt
dev.off()

#Fig.4h
Solid_genes <- c("Pbsn")
Solid_lineage_cds <- cds[rowData(cds)$gene_short_name %in% Solid_genes,
                         colData(cds)$stim %in% c("Triple")]

tiff(file = "Triple.combined1.BELEMETPos4 Pbsn in pseudotime.tiff", width = 5, height = 4, units = "in", compression = "lzw", res = 1000)
plot_genes_in_pseudotime(Solid_lineage_cds, cell_size = 1.5,
                         color_cells_by="EpiCellTypes"
) + scale_color_manual(values = c("salmon", "skyblue1", "olivedrab2", "brown3", "deeppink1",  "blue", "darkorange", "red", "turquoise3",
                                  "bisque3", "slategray3", "mediumorchid3", 
                                  "yellow2", "green4", "black"))
dev.off()

###Double
##Double.combined1.BELEMETPos re-clustering
DefaultAssay(Double.combined1.BELEMETPos) <- "RNA"
Double.combined1.BELEMETPos <- FindVariableFeatures(Double.combined1.BELEMETPos, selection.method = "vst", nfeatures = 5000)
Double.combined1.BELEMETPos <- ScaleData(Double.combined1.BELEMETPos, verbose = FALSE)
Double.combined1.BELEMETPos <- RunPCA(Double.combined1.BELEMETPos, npcs = 50, verbose = FALSE)
ElbowPlot(Double.combined1.BELEMETPos, ndims = 50)
Double.combined1.BELEMETPos <- FindNeighbors(Double.combined1.BELEMETPos, reduction = "pca", dims = 1:14)
Double.combined1.BELEMETPos <- FindClusters(Double.combined1.BELEMETPos, resolution = 0.5)
Double.combined1.BELEMETPos <- RunUMAP(Double.combined1.BELEMETPos, reduction = "pca", dims = 1:14)

#Cell cycle scoring
DefaultAssay(Double.combined1.BELEMETPos) <- "RNA"
all.genes <- rownames(Double.combined1.BELEMETPos)
Double.combined1.BELEMETPos <- ScaleData(Double.combined1.BELEMETPos, features = all.genes)
Double.combined1.BELEMETPos <- CellCycleScoring(Double.combined1.BELEMETPos, s.features = s.genes, g2m.features = g2m.genes, set.ident = TRUE)

#Cell Cycle regression
Double.combined1.BELEMETPos1 <- Double.combined1.BELEMETPos
DefaultAssay(Double.combined1.BELEMETPos1) <- "integrated"
Double.combined1.BELEMETPos1 <- ScaleData(Double.combined1.BELEMETPos1, vars.to.regress = c("S.Score", "G2M.Score"), features = rownames(Triple.combined1.BELEMETPos1))
Double.combined1.BELEMETPos1 <- RunPCA(Double.combined1.BELEMETPos1, features = VariableFeatures(Double.combined1.BELEMETPos1))
ElbowPlot(Double.combined1.BELEMETPos1, ndims = 50)
Double.combined1.BELEMETPos1 <- FindNeighbors(Double.combined1.BELEMETPos1, reduction = "pca", dims = 1:16)
Double.combined1.BELEMETPos1 <- FindClusters(Double.combined1.BELEMETPos1, resolution = 0.5)
Double.combined1.BELEMETPos1 <- RunUMAP(Double.combined1.BELEMETPos1, reduction = "pca", dims = 1:16)

#Convert Seurat to Monocle3 cell data set class
cds <- as.cell_data_set(Double.combined1.BELEMETPos1)

#Calculate size factors using built-in function in monocle3
cds <- estimate_size_factors(cds)

#Add gene names into CDS
cds@rowRanges@elementMetadata@listData[["gene_short_name"]] <- rownames(Double.combined1.BELEMETPos1)

#Construction single cell trajectories
cds <- cluster_cells(cds = cds, reduction_method = "UMAP") 

#plot single cell trajectory
cds <- learn_graph(cds, use_partition = FALSE)

#
plt <- plot_cells(cds,
                  color_cells_by = "EpiCellTypes",
                  show_trajectory_graph = TRUE,
                  label_branch_points = FALSE,
                  trajectory_graph_color = "black",
                  trajectory_graph_segment_size = 1.8,
                  label_cell_groups = FALSE,
                  label_leaves=FALSE,
                  label_roots = FALSE,
                  cell_size = 0.6,
                  alpha = 1) 
#SFig.4l
tiff(file = "Double.combined1.BELEMETPos1 re-clustering trajectory UMAP.tiff", width = 6, height = 5, units = "in", compression = "lzw", res = 800)
plt + scale_color_manual(values = c("salmon", "skyblue1", "olivedrab2", "brown3", "deeppink1",  "blue", "darkorange", "red", "turquoise3",
                                    "bisque3", "slategray3", "mediumorchid3", 
                                    "yellow2", "green4", "black"))
dev.off()

#plot single cell trajectory by pseudotime
cds <- order_cells(cds, reduction_method = "UMAP")
plt <- plot_cells(cds,
                  color_cells_by = "pseudotime",
                  show_trajectory_graph = TRUE,
                  label_branch_points = FALSE,
                  trajectory_graph_color = "black",
                  trajectory_graph_segment_size = 1.8,
                  label_cell_groups = FALSE,
                  label_leaves=FALSE,
                  label_roots = FALSE,
                  cell_size = 0.6,
                  alpha = 1)

#SFig.4l
tiff(file = "Double.combined1.BELEMETPos1 re-clustering pseudotime UMAP.tiff", width = 6, height = 5, units = "in", compression = "lzw", res = 800)
plt
dev.off()

#plot genes in pseudotime
Solid_genes <- c("Pbsn")
Solid_lineage_cds <- cds[rowData(cds)$gene_short_name %in% Solid_genes,
                         colData(cds)$stim %in% c("Double")]

#SFig.4m
tiff(file = "Double.combined1.BELEMETPos4 Pbsn in pseudotime.tiff", width = 5, height = 4, units = "in", compression = "lzw", res = 1000)
plot_genes_in_pseudotime(Solid_lineage_cds, cell_size = 1.5,
                         color_cells_by="EpiCellTypes", min_expr = 0.1,
) + scale_color_manual(values = c("salmon", "skyblue1", "olivedrab2", "brown3", "deeppink1",  "blue", "darkorange", "red", "turquoise3",
                                  "bisque3", "slategray3", "mediumorchid3", 
                                  "yellow2", "green4", "black"))
dev.off()

####RNAseq####
library(edgeR)
library(reshape2)
library(vegan)
library(rgl)
library(gplots)
library(grid)
library(gridExtra)
library(GenomicFeatures)
library(ggplot2)
library(statmod)
library(ComplexHeatmap)
library(circlize)
library(pheatmap)

setwd("//isi-dcnl/user_data/zjsun/group/Lab Members/Won Kyung Kim/hHGFtg-hMETtg-b-cat-Pb/RNAseq/CasVsPrimary")
countdata <- read.csv("countdata.csv")

##CPM
wk_RAW <- cbind(countdata$COHP_48624, countdata$COHP_48625, countdata$COHP_48626, countdata$COHP_48627, countdata$COHP_48628)
wk_RAW <- as.data.frame(wk_RAW)
colnames(wk_RAW) <- c("48625_Castration", "48626_Castration", "48627_Primary", "48628_Primary")
rownames(wk_RAW) <- row.names(countdata)
View(wk_RAW)

wk_RAW_DGEList <- DGEList(counts=wk_RAW[,1:4], group=c("Castration", "Castration", "Primary", "Primary"), genes=data.frame(Symbol=countdata$symbol, Length=countdata$total_exon_length))
colnames(wk_RAW_DGEList) <- c("48625_Castration", "48626_Castration", "48627_Primary", "48628_Primary")
names(wk_RAW_DGEList)
names(wk_RAW)
head(wk_RAW_DGEList$counts)

wk_RAW_DGEList$samples
head(wk_RAW_DGEList$genes)
w <- calcNormFactors(wk_RAW_DGEList, method = "TMM")
w$samples
CPM_w <- cpm(w)
dim(CPM_w)
View(CPM_w)

##TPM
normalized_RPKM <- CPM_w
normalized_RPKM[,1] <- normalized_RPKM[,1]*1000/countdata$total_exon_length
normalized_RPKM[,2] <- normalized_RPKM[,2]*1000/countdata$total_exon_length
normalized_RPKM[,3] <- normalized_RPKM[,3]*1000/countdata$total_exon_length
normalized_RPKM[,4] <- normalized_RPKM[,4]*1000/countdata$total_exon_length
keep_wk <- rowSums(normalized_RPKM >= 1) >= 1.4
sum(keep_wk)
w1 <- w[keep_wk, , keep.lib.sizes=FALSE]
w1$samples
colSums(w1$counts)
plotMDS(w1)
normalized_RPKM.keep <- normalized_RPKM[keep_wk,]
normalized_RPKM.keep.log2 <- log2(normalized_RPKM.keep+1)
write.table(normalized_RPKM.keep.log2, "normalized_RPKM.keep.log2.txt", quote=FALSE, sep="\t", col.names=NA)

##Further analysis
class(normalized_RPKM.keep.log2)
pca.values <- prcomp(na.omit(data.matrix(normalized_RPKM.keep.log2)))
pc.values <- data.frame(pca.values$rotation)
pc.values
variance.explained <- (pca.values $sdev)^2 / sum(pca.values $sdev^2)
pca.table <- data.frame(PC = 1:length(variance.explained), percent.variation = variance.explained, t(pc.values))
pca.text.file = paste("Group","_pca_values.txt",sep="")
write.table(pca.table, pca.text.file, quote=F, row.names=F, sep="\t")

qc.grp <- c("48625_Castration", "48626_Castration", "48627_Primary", "48628_Primary")
qc.grp
qc.grp <- as.character(w1$samples$group)
qc.grp
groups = levels(as.factor(as.character(qc.grp)))
groups
num.sample.types = length(groups)
num.sample.types
fixed.color.palatte = c("orange","purple","cyan","pink", colors())
class(fixed.color.palatte)
length(fixed.color.palatte)
color.palette <- fixed.color.palatte[1:num.sample.types]
color.palette
labelColors = rep("black",times=ncol(normalized_RPKM.keep.log2))
labelColors
for (i in 1:num.sample.types){
  labelColors[qc.grp == as.character(groups[i])] = color.palette[i]
}
labelColors
pca.file = paste("pca_by_","Group",".png",sep="")
png(file=pca.file)
plot(pc.values$PC1, pc.values$PC2, col = labelColors, xlab = paste("PC1 (",round(100* variance.explained[1] , digits = 2),"%)", sep = ""),ylab = paste("PC2 (",round(100* variance.explained[2] , digits = 2),"%)", sep = ""), pch=19)
text(pc.values$PC1, pc.values$PC2, c("48625_Castration", "48626_Castration", "48627_Primary", "48628_Primary"), pos = 1)
legend("bottomright", legend=groups, col=color.palette, pch=19)
dev.off()
groups

class(groups)
treatmentR3 <- factor(as.character(qc.grp))
treatmentR3
comp_10 <- relevel(treatmentR3, ref="Primary")
comp_10
comp_10_design <- model.matrix(~comp_10)
comp_10_design
rownames(comp_10_design) <- colnames(w1)
comp_10_design
w1_comp_10 <- estimateDisp(w1, comp_10_design, robust = TRUE)
w1_comp_10$common.dispersion
plotBCV(w1_comp_10)
fit_comp_10 <- glmQLFit(w1_comp_10, comp_10_design, robust=TRUE)
plotQLDisp(fit_comp_10)
names(fit_comp_10)
head(fit_comp_10$coefficients)

fit_comp_10$design
qlf_Cas_vs_Primary <- glmQLFTest(fit_comp_10, coef=2)
topTags(qlf_Cas_vs_Primary)

qlf_Cas_vs_Primary$table$normRPKM_48625_Castration <- normalized_RPKM[row.names(qlf_Cas_vs_Primary$table),1]
qlf_Cas_vs_Primary$table$normRPKM_48626_Castration <- normalized_RPKM[row.names(qlf_Cas_vs_Primary$table),2]
qlf_Cas_vs_Primary$table$normRPKM_48627_Primary <- normalized_RPKM[row.names(qlf_Cas_vs_Primary$table),3]
qlf_Cas_vs_Primary$table$normRPKM_48628_Primary <- normalized_RPKM[row.names(qlf_Cas_vs_Primary$table),4]

qlf_Cas_vs_Primary$table$RAW_48625_Castration <- w1_comp_10$counts[row.names(qlf_Cas_vs_Primary$table),1]
qlf_Cas_vs_Primary$table$RAW_48626_Castration <- w1_comp_10$counts[row.names(qlf_Cas_vs_Primary$table),2]
qlf_Cas_vs_Primary$table$RAW_48627_Primary <- w1_comp_10$counts[row.names(qlf_Cas_vs_Primary$table),3]
qlf_Cas_vs_Primary$table$RAW_48628_Primary <- w1_comp_10$counts[row.names(qlf_Cas_vs_Primary$table),4]

qlf_Cas_vs_Primary_data_frame <- as.data.frame(topTags(qlf_Cas_vs_Primary, n=30000))
sum(ifelse((qlf_Cas_vs_Primary_data_frame$logFC >= 1)&(qlf_Cas_vs_Primary_data_frame$FDR < 0.05),1,0))
sum(ifelse((qlf_Cas_vs_Primary_data_frame$logFC <= -1)&(qlf_Cas_vs_Primary_data_frame$FDR < 0.05),-1,0))

sum(ifelse((qlf_Cas_vs_Primary_data_frame$logFC >= log2(1.5))&(qlf_Cas_vs_Primary_data_frame$PValue < 0.05),1,0))
sum(ifelse((qlf_Cas_vs_Primary_data_frame$logFC <= -log2(1.5))&(qlf_Cas_vs_Primary_data_frame$PValue < 0.05),-1,0))

qlf_Cas_vs_Primary_data_frame$UP <- ifelse((qlf_Cas_vs_Primary_data_frame$logFC >= log2(1.5))&(qlf_Cas_vs_Primary_data_frame$PValue < 0.05),1,0)
qlf_Cas_vs_Primary_data_frame$DOWN <- ifelse((qlf_Cas_vs_Primary_data_frame$logFC <= -log2(1.5))&(qlf_Cas_vs_Primary_data_frame$PValue < 0.05),-1,0)

sum(qlf_Cas_vs_Primary_data_frame$UP)
sum(qlf_Cas_vs_Primary_data_frame$DOWN)
qlf_Cas_vs_Primary_data_frame$UP_OR_DOWN <- qlf_Cas_vs_Primary_data_frame$UP + qlf_Cas_vs_Primary_data_frame$DOWN
write.table(qlf_Cas_vs_Primary_data_frame, file="qlf_Cas_vs_Primary_data_frame_final.txt", quote=F, sep="\t", col.names=NA)


