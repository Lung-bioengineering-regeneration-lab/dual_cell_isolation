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
# DESeq analysis required for figures 4F and 4G.
# ********************************************************************* #
# Load count Data and sample information
#sample information
s2c <- read.table("data/s2c.csv", header = T, stringsAsFactors = F, sep = ",", row.names = 1)
#load counts 
counts <- readRDS("data/feature_counts_annotated.rds")

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

# ********************************************************************* #
# Fig4F
# ********************************************************************* #
# DO NOT FORGET TO LOAD ALL FUNCTIONS in Functions.R into the workspace environment.
# Run the differential exression for the 3DLD pellets vs. Organoids.
# 
results_orgs_3dld <- DE(dds, conditions = c("organoids_3DLD", "pellet_3DLD"), xlim = c(-20,20), labSize = 6.0)
orgs3dld <- results_orgs_3dld$res %>% na.omit()
orgs3dld$diffexpressed <- "NO"
orgs3dld[orgs3dld$pvalue<0.05 & abs(orgs3dld$log2FoldChange)> 1,]$diffexpressed <- "YES"

# This is just to correct the way the cells are annotated. And to combine what was refered to as activated cells from
# the original paper of the single cell dataset.
markers<- read.csv("data/HiRes_markers_0.5.csv", header = T, stringsAsFactors = F, sep = " ")
markers$cellType <- gsub("Activated ", "", markers$cellType)
markers$cellType <- gsub("ATII cells", "ATII", markers$cellType)
markers$cellType <- gsub("ATII", "ATII cells", markers$cellType)
markers$cellType <- gsub("Other ", "", markers$cellType)

## Combine DE list with labels for makers.
orgs3dld$ext_gene <- rownames(orgs3dld)
orgs3dld1 <- merge(orgs3dld, markers, by="ext_gene",all.x=T)


TPM <- readRDS("data/TPMs_gene_annotated.rds")
#filter TPMS with row sumns below 5 transcripts per milion in all samples combined.
TPM <- TPM[!(rowSums(TPM[,2:19]) < 5),]
#check duplicate ext genes
TPM <- TPM[!(duplicated(TPM$ext_gene)),] %>% remove_rownames() %>% column_to_rownames("ext_gene")

# filter out the TPMS for organoids only
orgs3dldTPM <- TPM[, colnames(TPM)%in%s2c[s2c$condition=="organoids_3DLD",]$sample]
# filter out the pellet TPMS only
pel3dldTPM <- TPM[, colnames(TPM)%in%s2c[s2c$condition=="pellet_3DLD",]$sample]
#calculate means for both.
orgs3dldTPM_means <- data.frame(row.names = rownames(TPM), "meanPel"= rowMeans(pel3dldTPM), "meanOrgs"= rowMeans(orgs3dldTPM))
orgs3dldTPM_means %>% rownames_to_column("ext_gene") -> orgs3dldTPM_means
# Combine differential expression data with the means of TPMS.
orgs3dld3 <- merge(orgs3dld1, orgs3dldTPM_means, by="ext_gene")
# generate a ggplot friendly data fram of the sig. genes that are NOT cell markers
orgs3dld3 %>% filter(diffexpressed=="YES"&is.na(cellType)==T) %>%
  select(contains(c("ext_gene", "cellType", "meanPel", "meanOrgs", "log2F"))) %>% 
  gather(key="differentiation", value = "TPM", -ext_gene, -cellType, -log2FoldChange) %>% 
  arrange(ext_gene) %>% arrange(cellType) -> orgs_TPMs_gray
# Factorize the gene names to maintain the order made in the previous line.
orgs_TPMs_gray$ext_gene <- factor(orgs_TPMs_gray$ext_gene, levels = unique(orgs_TPMs_gray$ext_gene))

## make dataframe with all significant genes that are ALSO cell marker genes.
orgs3dld3 %>% filter(diffexpressed=="YES"&is.na(cellType)==F) %>% 
  select(contains(c("ext_gene", "cellType", "meanPel", "meanOrgs", "log2F"))) %>% 
  gather(key="differentiation", value = "TPM", -ext_gene, -cellType, -log2FoldChange) %>% 
  arrange(ext_gene) %>% arrange(cellType) -> orgs_TPMs
orgs_TPMs$ext_gene <- factor(orgs_TPMs$ext_gene, levels = unique(orgs_TPMs$ext_gene))

######################## LETS DO the same for Classic organoids

results_orgs <- DE(dds, conditions = c("organoids_Classic", "pellet_Classic"), xlim = c(-20,20), labSize = 6.0)
orgs <- results_orgs$res %>% na.omit()
orgs$diffexpressed <- "NO"
orgs[orgs$pvalue<0.05 & abs(orgs$log2FoldChange)> 1,]$diffexpressed <- "YES"

# This is just to correct the way the cells are annotated. And to combine what was refered to as activated cells from
# the original paper of the single cell dataset.
orgs$ext_gene <- rownames(orgs)
orgs1 <- merge(orgs, markers, by="ext_gene",all.x=T)
# filter out the TPMS for organoids only
orgsTPM <- TPM[, colnames(TPM)%in%s2c[s2c$condition=="organoids_Classic",]$sample]
# filter out the pellet TPMS only
pelTPM <- TPM[, colnames(TPM)%in%s2c[s2c$condition=="pellet_Classic",]$sample]
#calculate means for both.
orgsTPM_means <- data.frame(row.names = rownames(TPM), "meanPel"= rowMeans(pelTPM), "meanOrgs"= rowMeans(orgsTPM))
orgsTPM_means %>% rownames_to_column("ext_gene") -> orgsTPM_means
# Combine differential expression data with the means of TPMS.
orgs3 <- merge(orgs1, orgsTPM_means, by="ext_gene")

# make dataframe with all significant genes that are ALSO cell marker genes.
orgs3 %>% filter(diffexpressed=="YES"&is.na(cellType)==F) %>% 
  select(contains(c("ext_gene", "cellType", "meanPel", "meanOrgs", "log2F"))) %>% 
  gather(key="differentiation", value = "TPM", -ext_gene, -cellType, -log2FoldChange) %>% 
  arrange(ext_gene) %>% arrange(cellType) -> orgs_TPMs_classic
orgs_TPMs_classic$ext_gene <- factor(orgs_TPMs_classic$ext_gene, levels = unique(orgs_TPMs_classic$ext_gene))

# generate a ggplot friendly data fram of the sig. genes that are NOT cell markers
orgs3 %>% filter(diffexpressed=="YES"&is.na(cellType)==T) %>%
  select(contains(c("ext_gene", "cellType", "meanPel", "meanOrgs", "log2F"))) %>% 
  gather(key="differentiation", value = "TPM", -ext_gene, -cellType, -log2FoldChange) %>% 
  arrange(ext_gene) %>% arrange(cellType) -> orgs_TPMs_classic_gray
# Factorize the gene names to maintain the order made in the previous line.
orgs_TPMs_classic_gray$ext_gene <- factor(orgs_TPMs_classic_gray$ext_gene, levels = unique(orgs_TPMs_classic_gray$ext_gene))

################################################# Lets prepare to comined both dataframes for 3dld and classic into 1.
colnames(orgs_TPMs) <- c("ext_gene", "cellType", "log2FC_3DLD", "differentiation", "TPM")
colnames(orgs_TPMs_gray) <- c("ext_gene", "cellType", "log2FC_3DLD", "differentiation", "TPM")

colnames(orgs_TPMs_classic) <- c("ext_gene", "cellType", "log2FC_Classic", "differentiation", "TPM")
colnames(orgs_TPMs_classic_gray) <- c("ext_gene", "cellType", "log2FC_Classic", "differentiation", "TPM")

## combine the two data frames.
orgsDE <- merge(orgs_TPMs_classic, orgs_TPMs, by="ext_gene",all=T) %>% mutate(cellType=coalesce(cellType.x, cellType.y))
orgsDE_gray <- merge(orgs_TPMs_classic_gray, orgs_TPMs_gray, by="ext_gene",all=T) %>% mutate(cellType=coalesce(cellType.x, cellType.y))

## remove duplicate genes.
orgsDE %>% filter(differentiation.x=="meanOrgs") %>% distinct(ext_gene, .keep_all = T) -> orgs_DE_d
orgsDE_gray %>% filter(differentiation.x=="meanOrgs") %>% distinct(ext_gene, .keep_all = T) -> orgs_DE_gray_d

orgs_DE_gray_d$TPM_mean <- rowMeans(data.frame("x"=orgs_DE_gray_d$TPM.x,"y"=orgs_DE_gray_d$TPM.y))
orgs_DE_d$TPM_mean <- rowMeans(data.frame("x"=orgs_DE_d$TPM.x,"y"=orgs_DE_d$TPM.y))

#### this filtration step gets rid of the NA in only one of both.
orgs_DE_gray_d %>% filter(TPM_mean >= 1) -> orgs_DE_gray_d_f
orgs_DE_d %>% filter(TPM_mean >= 1) -> orgs_DE_d_f

orgs_DE_d$NA.x <-NULL
orgs_DE_d$NA.y <- NULL
orgs_DE_all <- rbind(orgs_DE_d, orgs_DE_gray_d)
orgs_DE_all %>% filter(TPM_mean > 1) -> orgs_DE_all_f

## add this part only if you would like to label the dots in the 4th quadrant.
orgs_DE_all_f$lab <- NA
orgs_DE_all_f[orgs_DE_all_f$log2FC_Classic>0&orgs_DE_all_f$log2FC_3DLD<0,]$lab <- as.character(orgs_DE_all_f[orgs_DE_all_f$log2FC_Classic>0&orgs_DE_all_f$log2FC_3DLD<0,]$ext_gene)

ggplot(orgs_DE_gray_d_f, aes(x=log2FC_Classic, y=log2FC_3DLD)) + 
  geom_point(aes(size=log(TPM_mean)),color="grey85") +
  geom_point(data = orgs_DE_d_f,aes(x=log2FC_Classic, y=log2FC_3DLD, col=cellType,size=log(TPM_mean)))+
  #Add the following line to label dots in the 4th quadrant.
  #geom_text_repel(data = orgs_DE_all_f, mapping = aes(label=lab), color="red",show.legend = F, max.overlaps = 20)+
  theme_bw()+
  theme(text = element_text(size=20))+
  geom_vline(xintercept=c(0), size=1, linetype=1, color="black")+
  geom_hline(yintercept = c(0), size=1, linetype=1, color="black") + 
  guides(colour = guide_legend(override.aes = list(size=4)))-> orgsDEplot

orgsDEplot

ggsave(filename = "figure4G_legend_colvector.png",orgsDEplot,width = 8, height = 6)


# ********************************************************************* #
# Fig4G
# ********************************************************************* #
res_orgsClassic <- results_orgs$res

## Get lists for classic orgs vs pell
orgclassic_up <- rownames(res_orgsClassic[res_orgsClassic$log2FoldChange > 1 & res_orgsClassic$pvalue < 0.05,])
orgclassic_up <- orgclassic_up[-grep("NA", orgclassic_up)]
orgclassic_down <- rownames(res_orgsClassic[res_orgsClassic$log2FoldChange < -1 & res_orgsClassic$pvalue < 0.05,])
orgclassic_down <- orgclassic_down[-grep("NA", orgclassic_down)]

res_orgs3dld <- results_orgs_3dld$res
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



