# ********************************************************************* #
# Re-analyze HiRes dataset from Strunz. et al. figure 4 onwards
# ********************************************************************* #
# Data is obtained from publication: Strunz, Maximilian, et al. 
# "Alveolar regeneration through a Krt8+ transitional stem cell state 
# that persists in human lung fibrosis." Nature communications 11.1 
# (2020): 1-20.
# 

set.seed(0206)

# ----------------------------------------------------------------------
# *** Load Necessary Packages ***
# ----------------------------------------------------------------------
packages <- c("Seurat", "tidyverse", "ggplot2","Matrix", "ade4")
lapply(packages, library, character.only = TRUE)
# ----------------------------------------------------------------------

# Load existing Seurat object included in publication.
# This object needed further analysis.
HiRes <- readRDS("data/EpiHiRes_seurat.RDS")

# extract raw data
counts <- HiRes@assays$RNA@counts
ncol(counts)
features <- rownames(counts)
barcodes <-colnames(counts)
metadata <- HiRes@meta.data

# remove the HiRes Object to save memory
rm(HiRes)
gc()

# Creat New Seurat Object

HiRes <- CreateSeuratObject(counts = counts,
                                    project = "sc_HiRes",
                                    min.features = 200)

# apply metadata from previous
HiRes@meta.data$timepoint <- metadata$grouping
HiRes@meta.data$name <- metadata$name

## this data has already been cleared from mitochodrial cut offs. etc.
#VlnPlot(HiRes, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)

#plot1 <- FeatureScatter(HiRes, feature1 = "nCount_RNA", feature2 = "percent.mito")
#plot2 <- FeatureScatter(HiRes, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
#plot1 + plot2

HiRes <- SCTransform(object = HiRes, verbose = T)

HiRes <- FindVariableFeatures(HiRes, verbose = T, nfeatures = 3000)

HiRes <- ScaleData(HiRes, features = row.names(HiRes@assays$SCT@data))
HiRes <- RunPCA(HiRes)

# function has been created to perform Mantel test and elbow plots, see Functions.R
HiResPCs <- evaluatePCs(HiRes, 1000, 8, 20)
saveRDS(HiResPCs$mantel, file = "HiRes_mantel1.rds")
saveRDS(HiResPCs$elbow, file = "HiRes_elbow1.rds")

###############################################################################

set.seed(0206)

# Use PC chosen based on Mantel and elbow  to run analysis: PC12 was chosen
HiRes <- RunUMAP(object = HiRes, dims = 1:50, verbose = T)
HiRes <- FindNeighbors(object = HiRes, dims = 1:30, verbose = T)
HiRes <- FindClusters(object = HiRes, resolution = 1, verbose = T)

pdf("Dimplot_post_umap.pdf")
DimPlot(HiRes)
dev.off()



#### depending on the seed you have used, the number of clusters may be altered. 
# Here we use DotPlot function to examine expression of matched markers from 
# each cluster to what has been published with the dataset. 
# Consult the webtool published by original authors to select markers for clusters:
# https://theislab.github.io/LungInjuryRegeneration/

Idents(HiRes,cells = WhichCells(HiRes, idents = c(0))) <- "Activated ATII"
Idents(HiRes,cells = WhichCells(HiRes, idents = c(1,3))) <- "ATII cells"
Idents(HiRes,cells = WhichCells(HiRes, idents = c(6))) <- "ATI cells"
Idents(HiRes,cells = WhichCells(HiRes, idents = c(17))) <- "Basal Cells"
Idents(HiRes,cells = WhichCells(HiRes, idents = c(5))) <- "Krt8+ ADI"
Idents(HiRes,cells = WhichCells(HiRes, idents = c(10,4,12,16))) <- "Ciliated Cells"
Idents(HiRes,cells = WhichCells(HiRes, idents = c(7, 13, 14))) <- "Activated Ciliated Cells"

Idents(HiRes,cells = WhichCells(HiRes, idents = c(2))) <- "Club cells"
Idents(HiRes,cells = WhichCells(HiRes, idents = c(15))) <- "Activated Club cells"
Idents(HiRes,cells = WhichCells(HiRes, idents = c(8))) <- "Other Activated Club cells"
Idents(HiRes,cells = WhichCells(HiRes, idents = c(9))) <- "Goblet cells"

Idents(HiRes,cells = WhichCells(HiRes, idents = c(18))) <- "Neuroendocrine Cells"

Idents(HiRes,cells = WhichCells(HiRes, idents = c(11))) <- "Mki67+ Proliferation"

x <- FindAllMarkers(HiRes)

x <- split(x, x$cluster)

HiRes@meta.data$population <- Idents(HiRes)

DimPlot(HiRes)



# mark the PBS population

HiRes@meta.data$treatment <- "Bleomycin"
HiRes@meta.data[HiRes@meta.data$timepoint=="d14_PBS",]$treatment <- "PBS"


#combine activated cells with its own cell type for simplicity.
p1 <- DimPlot(HiRes, pt.size = 1) + NoAxes()
p2 <- DimPlot(HiRes, group.by = "treatment", cols = c("lightgrey", "blue"),pt.size = 1) + ggtitle("") + NoAxes()

g1 <- ggplotGrob(p1)
g2 <- ggplotGrob(p2)

g2$widths = g1$widths

grid.arrange(g1, g2, ncol = 2)

### Figure 4B.
ggsave(p1, filename = "Fig4B.png")

### Subset PBS for use with deconvolution
HiResPBS <- subset(HiRes, subset = timepoint == "d14_PBS")







