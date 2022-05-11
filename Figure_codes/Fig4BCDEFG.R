# ********************************************************************* #
# Figure 4 Deconvolution of bulk RNAseq of distal progenitor pellets 
# predicts higher proportions of contaminating proximal airway cells 
# in the classic isolation.
# ********************************************************************* #
# NOTE: YOU NEED TO RUN ON A HIGH PERFORMANCE COMPUTATION CLUSTER
# ********************************************************************* #
# Fig4B. Preparation of the sc Dataset
# ********************************************************************* #

# please refer to R files in relation to the dataset:
# Reanalysis_of_HiRes_scRNAseq.R

# ********************************************************************* #
### load necessary packages
pckgs <- c("BisqueRNA", "Biobase", "ggplot2", "plotly","Matrix", "tidyverse", "Seurat", "ade4")
lapply(pckgs, library, character.only=TRUE)

# ********************************************************************* #
# 1. Prepartion of bulk ESET 
# ********************************************************************* #

## Prepare the bulkRNAseq expressionSet

bulk_counts <- readRDS("data/feature_counts_clean.rds")
s2c <- read.csv("data/s2c.csv")
s2c$X <- NULL


### select samples of interest
smpls <- s2c[s2c$cell=="ATII" & s2c$sampleType=="pellet",]
bulk_counts <- bulk_counts[,colnames(bulk_counts) %in% smpls$sample]

### remove zero rows 
bulk_counts <- bulk_counts[rowSums(bulk_counts)>0,]

### to be considered expressed, a gene has to be present in 2 out of 3 samples
groups <- as.character(unique(smpls$condition))
expressed <- bulk_counts[,0]

for (i in 1: length(groups)){
  ex <- bulk_counts[,colnames(bulk_counts)%in%smpls[smpls$condition==groups[i],]$sample]
  zeros <- apply(ex, 1, function(x) length(which(x==0)))
  flamingo <- rep(1,nrow(bulk_counts))
  flamingo[zeros>1]<-0
  table(flamingo)
  expressed <- cbind(expressed,flamingo)
  colnames(expressed)[i] <- groups[i]
}
head(expressed)
head(bulk_counts)

## remove the genes with the profile given
bulk_counts <- bulk_counts[rowSums(expressed)>0,]

## create bulk matrix
bulk.matrix <- round(as.matrix(bulk_counts))

head(bulk.matrix)
### create an expression set:
bulk.eset <- Biobase::ExpressionSet(assayData = bulk.matrix)

# ********************************************************************* #
#  generate annotation table
# ********************************************************************* #

library(biomaRt)
mart <- biomaRt::useMart(biomart = "ensembl", dataset = "mmusculus_gene_ensembl")
t2g <- biomaRt::getBM(attributes = c("ensembl_transcript_id", "ensembl_gene_id",
                                     "external_gene_name", "transcript_biotype"), mart = mart)
t2g <- dplyr::rename(t2g, target_id = ensembl_transcript_id,
                     ens_gene = ensembl_gene_id, ext_gene = external_gene_name, transcript_type = transcript_biotype)

# ********************************************************************* #
#  Prepare single_cell dataset for Bisque
# ********************************************************************* #

set.seed(0406)
raw.sc <- GetAssayData(HiResPBS, slot = "counts")
data.frame("ext_gene"=rownames(raw.sc)) -> sc.genes
mart <- biomaRt::useMart(biomart = "ensembl", dataset = "mmusculus_gene_ensembl")
e2e <- biomaRt::getBM(attributes = c("ensembl_gene_id",
                                     "external_gene_name"),
                      values = sc.genes$ext_gene, 
                      filters = c("external_gene_name"), 
                      mart = mart,uniqueRows = T)
e2e <- dplyr::rename(e2e, ens_gene = ensembl_gene_id, ext_gene = external_gene_name)

e2e %>% distinct(ext_gene, .keep_all = T) -> e2e2
raw.sc <- raw.sc[rownames(raw.sc) %in% e2e2$ext_gene,]
e2e2 <- e2e2[e2e2$ext_gene %in% rownames(raw.sc),]
e2e2 <- e2e2[match(rownames(raw.sc), e2e2$ext_gene),]

rownames(raw.sc) <- e2e2$ens_gene
raw.sc <- as.matrix(raw.sc)

sc.pheno <- data.frame(check.names=F, check.rows=F,
                       stringsAsFactors=F,
                       row.names=rownames(HiResPBS@meta.data),
                       SubjectName=HiResPBS@meta.data$orig.ident,
                       cellType=HiResPBS@meta.data$population)

sc.meta <- data.frame(labelDescription=c("SubjectName",
                                         "cellType"),
                      row.names=c("SubjectName",
                                  "cellType"))
sc.pdata <- new("AnnotatedDataFrame",
                data=sc.pheno,
                varMetadata=sc.meta)


sc.eset <- Biobase::ExpressionSet(assayData=raw.sc,
                                  phenoData=sc.pdata)


### RUN Bisque

res <- BisqueRNA::ReferenceBasedDecomposition(bulk.eset = bulk.eset
                                              ,sc.eset = sc.eset 
                                              ,use.overlap = FALSE
                                              ,old.cpm=F)
res$genes.used

#special coloring
special.colors <- c("purple", "goldenrod")
smpls <- s2c[s2c$cell=="ATII" & s2c$sampleType=="pellet",]

# ********************************************************************* #
# Fig4C.
# ********************************************************************* #

library(ggplot2)
library(reshape2)
library(tidyverse)


smpls <- s2c[s2c$cell=="ATII" & s2c$sampleType=="pellet",]
RBE <- as.data.frame(res$bulk.props)

#combine activated cells with their corresponding cell type for simplicity
RBE[rownames(RBE)=="Ciliated Cells",] <- RBE[rownames(RBE)=="Ciliated Cells",]+RBE[rownames(RBE)=="Activated Ciliated Cells",]
RBE <- RBE[!(rownames(RBE)=="Activated Ciliated Cells"),]

RBE[rownames(RBE)=="Club cells",] <- RBE[rownames(RBE)=="Club cells",]+RBE[rownames(RBE)=="Activated Club cells",]+RBE[rownames(RBE)=="Other Activated Club cells",]
RBE <- RBE[!(rownames(RBE)=="Activated Club cells"),]
RBE <- RBE[!(rownames(RBE)=="Other Activated Club cells"),]

RBE[rownames(RBE)=="ATII cells",] <- RBE[rownames(RBE)=="ATII cells",]+RBE[rownames(RBE)=="Activated ATII",]
RBE <- RBE[!(rownames(RBE)=="Activated ATII"),]


#convert to dataframe for ggplot2
RBE.df <- data.frame(RBE)
#transpose dataframe for use
RBE.df.t <- as.data.frame(t(RBE.df))

RBE.df.t$sample <- rownames(RBE.df.t)
RBE.df.t$condition <- factor(smpls$condition, levels = c("pellet_Classic","pellet_3DLD"))

RBE.df.t %>% gather(key="celltype", value = "value", -sample, -condition) -> RBE.df.t.g
RBE.df.t.g$value <- as.numeric(RBE.df.t.g$value)

RBE.df.t.g$label <- factor(paste0(gsub("pellet_","", RBE.df.t.g$condition), "_", rep(c(1,2,3),18)), levels = unique(paste0(gsub("pellet_","", RBE.df.t.g$condition), "_", rep(c(1,2,3),18))))


library(RColorBrewer)
n <- 9
qual_col_pals <- brewer.pal.info[brewer.pal.info$category == 'qual',]
col_vector = unlist(mapply(brewer.pal, qual_col_pals$maxcolors, rownames(qual_col_pals)))[1:n]
col_vector[4] <- "#00ACC1"
c <- ggplot(RBE.df.t.g, aes(fill=celltype, y=100*value, x=label)) + 
  ylab("Cell Proportion (%)")+ xlab("")+
  geom_bar(position = "stack",stat = "identity") + 
  theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
  theme(axis.text.x = element_text(color = c("red", "green", "blue", "brown", "orange", "black"))) +
  theme(panel.background = element_rect(fill = "white", size = 0.5, color = "black")) +
  theme(panel.grid.major = element_line(colour = "gray90")) + 
  theme(panel.grid.minor = element_line(colour = "gray90")) +
  scale_fill_manual(values=col_vector) + 
  theme(text = element_text(size=25))
c
ggsave("Fig4C.png",plot = c, width = 7, height = 8)

# ********************************************************************* #
# Fig4D.
# ********************************************************************* #

b <- ggplot(RBE.df.t.g, aes(color=sample,fill=condition, y=100*value, x=condition)) + 
  geom_violin(aes(colour=condition, fill=condition), alpha=0.3)+
  geom_jitter(aes(fill=sample),stat="identity", width = 0.25, size=3) + ylab("Cell Proportions (%)") + xlab("") +
  facet_wrap(~celltype, scales = scales, ncol = 2) +
  theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
  guides(fill=FALSE) +
  expand_limits(y=0) +
  theme(panel.background = element_rect(fill = "white", size = 0.5, color = "black")) +
  theme(panel.grid.major = element_line(colour = "white")) + 
  theme(panel.grid.minor = element_line(colour = "grey")) + theme(legend.position="none")+
  scale_color_manual(values = cols) + scale_fill_manual(values=cols) + theme(text = element_text(size=25))
b
ggsave("Fig4D.png",plot = b, width = 8, height = 17)


# ********************************************************************* #
# Fig4E.
# ********************************************************************* #

tpm <- readRDS("data/TPMs_gene_annotated.rds")
head(tpm)

s2c <- read.csv("data/s2c.csv")
s2c$X <- NULL


### select samples of interest
smpls <- s2c[s2c$cell=="ATII",]
tpm[,colnames(tpm) %in% c(smpls$sample,"ens_gene")] %>% column_to_rownames("ens_gene") -> ATII_tpm
ATII_tpm %>% head()

### remove zero rows 
ATII_tpm <- ATII_tpm[rowSums(ATII_tpm)>0,]

### to be considered expressed, a gene has to be present in 2 out of 3 samples
groups <- as.character(unique(smpls$condition))
expressed <- ATII_tpm[,0]

for (i in 1: length(groups)){
  ex <- ATII_tpm[,colnames(ATII_tpm)%in%smpls[smpls$condition==groups[i],]$sample]
  zeros <- apply(ex, 1, function(x) length(which(x==0)))
  flamingo <- rep(1,nrow(ATII_tpm))
  flamingo[zeros>1]<-0
  table(flamingo)
  expressed <- cbind(expressed,flamingo)
  colnames(expressed)[i] <- groups[i]
}

## remove the genes with the profile given
ATII_tpm <- ATII_tpm[rowSums(expressed)>0,]


# ********************************************************************* #
#  generate annotation table
# ********************************************************************* #

library(biomaRt)
mart <- biomaRt::useMart(biomart = "ensembl", dataset = "mmusculus_gene_ensembl")
t2g <- biomaRt::getBM(attributes = c("ensembl_transcript_id", "ensembl_gene_id",
                                     "external_gene_name", "transcript_biotype"), mart = mart)
t2g <- dplyr::rename(t2g, target_id = ensembl_transcript_id,
                     ens_gene = ensembl_gene_id, ext_gene = external_gene_name, transcript_type = transcript_biotype)

# ********************************************************************* #
#  Prepare single_cell dataset - See Reanalysis_of_HiRes_scRNAseq.R
# ********************************************************************* #
#  Prepare bisque markers
# --------------------------------------------------------------------- #
#generate gene lists for each cluster
x1 <- FindAllMarkers(HiRes,only.pos = T)
# Filter with avg logFC of 0.5
y1 <- x1[x1$avg_log2FC>0.5,]

data.frame("ext_gene"=y1$gene,"cellType"=y1$cluster) %>% 
  merge(., t2g, by = "ext_gene") %>% 
  distinct(ext_gene, .keep_all = T) %>% arrange(cellType) -> HR_markers

HR <- summary(HR_markers$cellType)

#write.table(HR_markers, "HiRes_markers_0.5.csv")


ATII_filtered <- ATII_tpm[rownames(ATII_tpm)%in%HR_markers$ens_gene,]

head(ATII_filtered)

#convert annotation to gene names
ATII_filtered$ens_gene <- rownames(ATII_filtered)

HR_markers %>% select(contains("gene")) %>% 
  merge(.,ATII_filtered, by="ens_gene") %>% 
  column_to_rownames("ext_gene") %>% select(-contains("gene")) -> ATII_f_annotated


######################################################################
# generate heatmap
library(pheatmap)

atii_filtered_log <- log(ATII_f_annotated)
atii_filtered_log[atii_filtered_log==-Inf]<-NA


library(pheatmap)
library(RColorBrewer)

rownames(s2c) <- s2c$sample


annot.col <- s2c[s2c$cell %in% c("ATII"),] %>% select(contains("condition")) %>% arrange(condition)

annot.row <- HR_markers %>% column_to_rownames("ext_gene") %>% arrange(cellType) %>% select(cellType)

annot.row1 <- annot.row

unique(annot.row1$cellType)

annot.row1$cellType <- gsub("Other ", "", annot.row1$cellType)
annot.row1 %>% arrange(cellType) -> annot.row1

atii_filtered_log <- atii_filtered_log[rownames(annot.row1),]
atii_filtered_log <- atii_filtered_log[,rownames(annot.col)]

library(heatmaply)


#### Calcuate Z score for heatmap

atIImax <- apply(ATII_f_annotated, 1, max)
atIImin <- apply(ATII_f_annotated, 1, min)

atIImean <- apply(ATII_f_annotated, 1, mean)
atIIsd <- apply(ATII_f_annotated, 1, sd)

ATII <- (ATII_f_annotated - atIImean)/(atIIsd)

ATII <- ATII[rownames(annot.row1),]
ATII <- ATII[,rownames(annot.col)]

ATIIFIXED <- ATII

#### limit the range between -2, 2 for the z-score.

ATIIFIXED[ATIIFIXED< -2]<- -2
ATIIFIXED[ATIIFIXED> 2]<- 2


##Figure 4E
a <-heatmaply(ATIIFIXED,
              row_side_colors = annot.row1,
              col_side_colors = annot.col,
              dendrogram = F, showticklabels = F, 
              row_dend_left = T
)
a[["x"]][["layout"]][["annotations"]][[2]][["text"]] <- ""
a


# ********************************************************************* #
# Fig4F
# ********************************************************************* #

results_orgs_3dld <- DE(dds, conditions = c("organoids_3DLD", "pellet_3DLD"), xlim = c(-20,20), labSize = 6.0)
results_orgs_3dld$volcano
results_orgs_3dld$res

orgs3dld <- results_orgs_3dld$res
orgs3dld <- na.omit(orgs3dld)
orgs3dld$delabel <- NA
orgs3dld$diffexpressed <- "NO"
orgs3dld[orgs3dld$pvalue<0.05 & abs(orgs3dld$log2FoldChange)> 1,]$diffexpressed <- "YES"

orgs3dld$delabel[orgs3dld$diffexpressed != "NO"] <- rownames(orgs3dld[orgs3dld$diffexpressed != "NO",])


markers<- read.csv("HiRes_markers_0.5.csv", header = T, stringsAsFactors = F, sep = " ")

markers$cellType <- gsub("Activated ", "", markers$cellType)
markers$cellType <- gsub("ATII cells", "ATII", markers$cellType)
markers$cellType <- gsub("ATII", "ATII cells", markers$cellType)
markers$cellType <- gsub("Other ", "", markers$cellType)



orgs3dld$ext_gene <- rownames(orgs3dld)

orgs3dld1 <- merge(orgs3dld, markers, by="ext_gene",all.x=T)


orgs3dld1$label <- NA


orgs3dld1[orgs3dld1$diffexpressed=="YES",]$label <- "Differentially Expressed"

orgs3dld1[orgs3dld1$diffexpressed=="YES" & is.na(orgs3dld1$cellType)==FALSE,]$label <- "Cell Marker"

orgs3dld1

summary(as.factor(orgs3dld1$label))

orgs3dld1$labelCells <- NA
orgs3dld1[orgs3dld1$diffexpressed=="YES",]$labelCells <- "Differentially Expressed"
orgs3dld1[orgs3dld1$diffexpressed=="YES" & is.na(orgs3dld1$cellType)==FALSE,]$labelCells <- orgs3dld1[orgs3dld1$diffexpressed=="YES" & is.na(orgs3dld1$cellType)==FALSE,]$cellType

orgs3dld2 <- orgs3dld1[orgs3dld1$diffexpressed=="YES" & is.na(orgs3dld1$cellType)==FALSE,]

TPM <- readRDS("data/TPMs_gene_annotated.rds")
head(counts)

TPM <- TPM[!(rowSums(TPM[,2:19]) < 5),]

#check duplicate ext genes
TPM <- TPM[!(duplicated(TPM$ext_gene)),]
rownames(TPM) <- TPM$ext_gene
TPM %>% select(-contains("gene")) ->TPM

orgs3dldTPM <- TPM[, colnames(TPM)%in%s2c[s2c$condition=="organoids_3DLD",]$sample]
orgs3dldTPM$mean <- rowMeans(orgs3dldTPM)

orgs3dldTPM %>% select(contains("mean")) %>% rownames_to_column("ext_gene") -> orgs3dldTPM

orgs3dld3 <- merge(orgs3dld1, orgs3dldTPM, by="ext_gene")
orgs3dld4 <- orgs3dld3[orgs3dld3$diffexpressed=="YES" & is.na(orgs3dld3$cellType)==FALSE,]

orgs3dld4b <- orgs3dld3[is.na(orgs3dld3$cellType)==TRUE,]


###################### FIGURE 4F
ggplot(data=orgs3dld4b, aes(x=log10(mean), y=log2FoldChange)) + 
  geom_point(color="grey85") +
  geom_point(data = orgs3dld4,aes(x=log10(mean), y=log2FoldChange, col=labelCells))+
  geom_text_repel(data = orgs3dld4, mapping = aes(label=ext_gene, color=labelCells),show_guide = F)+
  xlab("log(TPM[organoids_average])")+
  ylab("log2FC (organoids vs. pellet))")+
  theme_minimal()

######### SAME CODE IS APPLIED AS ABOVE FOR CLASSIC ISOLATION> 
######### JUST CHANGE THE DE TO COMPARE c("organoids_Classic", "pellet_Classic")



# ********************************************************************* #
# Fig4G
# ********************************************************************* #
# Load count Data and sample information

#sample information
s2c <- read.table("data/s2c.csv", header = T, stringsAsFactors = F, sep = ",", row.names = 1)

#load counts 
counts <- readRDS("data/feature_counts_annotated.rds")
head(counts)

counts <- counts[!(rowSums(counts[,2:19]) < 5),]

#check duplicate ext genes
counts <- counts[!(duplicated(counts$ext_gene)),]

rownames(counts) <- counts$ext_gene

counts %>% select(-contains("gene")) ->counts

#order column names same as metadata
counts <- counts[,s2c$sample]


# DEseq

s2c$condition <- as.factor(s2c$condition)
# this last line is to set the reference condition
s2c$condition <- relevel(s2c$condition, ref = "pellet_Classic")

dds <- DESeqDataSetFromMatrix(counts, colData = s2c, design = ~condition)
dds <- estimateSizeFactors( dds )
sizeFactors( dds )

dds<-estimateDispersions(dds)
keep <- rowSums(counts(dds))>=1
dds <- dds[keep,]

dds$condition <- relevel(dds$condition, ref="pellet_3DLD") 

dds <- DESeq(dds)


results_pell_newVsold <- DE(dds, conditions = c("pellet_Classic", "pellet_3DLD"), xlim = c(-13,13), labSize = 6.0)
results_pell_newVsold$volcano


results_pell_newVsold$res

res_pellets <- results_pell_newVsold$res


#### added on 2022-04-13
classic_pel_up <- rownames(res_pellets[res_pellets$log2FoldChange >= 1 & res_pellets$pvalue < 0.05,])
classic_pel_up <- classic_pel_up[-grep("NA", classic_pel_up)]
classic_pel_down <- rownames(res_pellets[res_pellets$log2FoldChange <= -1 & res_pellets$pvalue < 0.05,])
classic_pel_down <- classic_pel_down[-grep("NA", classic_pel_down)]



#### added on 2022-04-13
results_orgs_classic <- DE(dds, conditions = c("organoids_Classic","pellet_Classic"), xlim = c(-13,13), labSize = 6.0)
results_orgs_classic$volcano
results_orgs_classic$res

res_orgsClassic <- results_orgs_classic$res


#### added on 2022-04-13
orgclassic_up <- rownames(res_orgsClassic[res_orgsClassic$log2FoldChange > 1 & res_orgsClassic$pvalue < 0.05,])
orgclassic_up <- orgclassic_up[-grep("NA", orgclassic_up)]
orgclassic_down <- rownames(res_orgsClassic[res_orgsClassic$log2FoldChange < -1 & res_orgsClassic$pvalue < 0.05,])
orgclassic_down <- orgclassic_down[-grep("NA", orgclassic_down)]






results_orgs_3dld <- DE(dds, conditions = c("organoids_3DLD", "pellet_3DLD"), xlim = c(-13,13), labSize = 6.0, ylim = c(0,150))
results_orgs_3dld$volcano
results_orgs_3dld$res

res_orgs3dld <- results_orgs_3dld$res

grep("NA", rownames(res_orgs3dld))

orgs3dld_up <- rownames(res_orgs3dld[res_orgs3dld$log2FoldChange > 1 & res_orgs3dld$pvalue < 0.05,])
orgs3dld_up <- orgs3dld_up[-grep("NA", orgs3dld_up)]


orgs3dld_down <- rownames(res_orgs3dld[res_orgs3dld$log2FoldChange < -1 & res_orgs3dld$pvalue < 0.05,])
orgs3dld_down <- orgs3dld_down[-grep("NA", orgs3dld_down)]

x <- list("Up_Classic_Organoids"=orgclassic_up, "Down_Classic_Organoids"= orgclassic_down, "Up_3DLD_Organoids"=orgs3dld_up, "Down_3DLD_organoids"= orgs3dld_down)

a <- ggvenn(
  x, 
  fill_color = c("#0073C2FF", "#EFC000FF", "#868686FF", "#CD534CFF"),
  stroke_size = 0.75, show_percentage = F, text_size = 6
)

#### this file was saved and formated properly on illustrator to adjust labels
ggsave(filename = "VennDiagram_DEclassicVs3DLDorgs.svg", plot = a)



