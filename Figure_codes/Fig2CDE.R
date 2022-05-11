# ********************************************************************* #
# Figure 2C. PCA of proximal progenitors before and after ALI
# ********************************************************************* #
# NOTE:
# Some functions used in this file are defined in Functions.R. 
# MAKE SURE TO RUN ALL CODE IN Functions.R before running anything else.
# Load packages

pckgs <- c("ggplot2", "plotly","ggbiplot", "tidyverse")
lapply(pckgs, library, character.only=TRUE)


# ********************************************************************* #
# Load TPM data for MTECs pellet and ALI
# ********************************************************************* #

### TPMs generated previously.

tpm <- readRDS("data/TPMs_gene_annotated.rds") %>% column_to_rownames("ens_gene")

#load sample
s2c <- read.table("data/s2c.csv", header = T, row.names = 1, sep = ",")

#ensure correct nomenclature
s2c[s2c$cell == "MTEC",]$cell <- "Proximal Progenitors"

tpm[,colnames(tpm) %in% c(s2c[s2c$condition %in% c("pellet_MTEC", "ALI_MTEC"),]$sample,"ext_gene")] -> mtecs


# ********************************************************************* #
# Perform PCA using all genes
# ********************************************************************* #

keeping1 <- rowSums(mtecs[,1:6])>=1
mtecs <- mtecs[keeping1,]

PCA.data <- t(mtecs[,1:6])
colnames(PCA.data)<-mtecs$ext_gene

labels <- s2c[s2c$condition %in% c("pellet_MTEC", "ALI_MTEC"),]$sampleName
groups <- s2c[s2c$condition %in% c("pellet_MTEC", "ALI_MTEC"),]$condition
groups <- rep(c("ALI_Dy28", "MTEC_Day0"),each=3)

matrix.pca <- prcomp(PCA.data, center = TRUE, scale. = T) 

g <- ggbiplot(matrix.pca, scale = 1, pc.biplot = TRUE, obs.scale = 1, var.scale = 1, groups = groups, ellipse = TRUE, circle = TRUE, var.axes = F, labels = NULL, labels.size = 1 , alpha = 0) 

g <- g + geom_point(aes(color = groups), size = 5, alpha = 0.8)
g <- g + scale_color_manual(values = c("goldenrod", "purple") )
g <- g + theme(legend.direction = 'horizontal', legend.position = 'top')
g <- g + theme(axis.line = element_line(colour = "black"))
g <- g + theme(panel.background = element_rect(fill = "white", size = 1, color = "black"))
g <- g + theme(panel.grid.major = element_line(colour = "gray90")) 
g <- g + theme(legend.key = element_rect(fill="white"))
g <- g + geom_vline(xintercept = 0, linetype = "dotdash")
g <- g + geom_hline(yintercept = 0, linetype = "dotdash") + 
  theme(text = element_text(size=20))

ggsave("Fig2C.png",g)




# ********************************************************************* #
# Figure 2D. Differential expression proximal progenitors before and after ALI
# ********************************************************************* #
# NOTE:
# Some functions used in this file are defined in Functions.R. 
# MAKE SURE TO RUN ALL CODE IN Functions.R before running anything else.

# Load packages

pckgs <- c("DESeq2","EnhancedVolcano",  "Biobase", "ggplot2", "plotly","Matrix",
           "tidyverse", "Seurat", "ade4")
lapply(pckgs, library, character.only=TRUE)

# ********************************************************************* #
# Load Data
# ********************************************************************* #

#load sample metadata
s2c <- read.table("data/s2c.csv", header = T, row.names = 1, sep = ",")


#load counts 
counts <- readRDS("data/feature_counts_annotated.rds")
head(counts)

#filter counts. Remove all rows where the sum of all samples is less than counts
counts <- counts[!(rowSums(counts[,2:19]) < 5),]

#check duplicate ext genes
counts <- counts[!(duplicated(counts$ext_gene)),]

rownames(counts) <- counts$ext_gene

counts %>% select(-contains("gene")) ->counts

#order column names same as metadata
counts <- counts[,s2c$sample]


# ********************************************************************* #
# Differential expression with DEseq
# ********************************************************************* #

s2c$condition <- as.factor(s2c$condition)

# this last line is to set the reference condition
s2c$condition <- relevel(s2c$condition, ref = "pellet_MTEC")

dds <- DESeqDataSetFromMatrix(counts, colData = s2c, design = ~condition)
dds <- estimateSizeFactors( dds )
sizeFactors( dds )

dds<-estimateDispersions(dds)
keep <- rowSums(counts(dds))>=1
dds <- dds[keep,]
dds <- DESeq(dds)

### Function DE is defined in Functions.R. MAKE SURE TO RUN ALL CODE IN Functions.R before running this.
results_MTEC1 <- DE(dds, conditions = c("ALI_MTEC", "pellet_MTEC"), xlim = c(-15,25),labSize=6)

results_MTEC1$volcano
ggsave(filename = "Fig2D.png", plot = results_MTEC1$volcano)

# ********************************************************************* #
# Figure 2E. Heatmap of Proximal epithelial markers
# ********************************************************************* #


genes <- c("Cckar", "Clic5","Ccdc78", "Cyp2f2",  "Dapl1","Foxj1", "Muc1",  "Muc5ac", "Muc5b", "Krt14", "Krt5","Trp63","Trp73", "Scgb3a2")

#Function plotHeatmap is defined in Functions.R
plotHeatmap(genes, f=c("pellet_MTEC","ALI_MTEC"), data=tpm, labs=c("condition"), h=8, w=6, fil="Fig2E.png")


