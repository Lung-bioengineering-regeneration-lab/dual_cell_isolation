# ********************************************************************* #
# Figure 3 Load and prepare organoid data
# ********************************************************************* #

#load all required packages
pckgs <- c("ggplot2", "ggsignif","ggbiplot", "tidyverse")
lapply(pckgs, library, character.only=TRUE)

organoids <- read.csv("data/atii_organoid_sizes.csv", header = T, stringsAsFactors = F, sep = ";", dec = ".")

organoids[organoids$isolation=="new",]$isolation <- "3DLD"
organoids[organoids$isolation=="old",]$isolation <- "Classic"


head(organoids)
organoids$experiment <- as.factor(organoids$experiment)

organoids$isolation <- factor(organoids$isolation, levels = c("Classic","3DLD"))
#convert area in um

organoids$area_um <- organoids$Area * 1000000



# ********************************************************************* #
# Figure 3C. Organoid CFE
# ********************************************************************* #
# draw rectangle of known organoid CFE range (0.5 - 2%)
rect <- data.frame(xmin=-Inf, xmax=Inf, ymin=0.5, ymax=2)


per_well <- organoids %>% group_by(well_id) %>% select(contains("area_um")) %>% summarize_all(funs(mean,sd, length))

colnames(per_well) <- c("well_id", "Organoid_size_average", "organoid_size_sd", "Organoids_per_well")

per_well$isolation <- c(rep("3DLD",2),rep("Classic",2),rep("3DLD",3),rep("Classic",3))


per_well$isolation <- factor(per_well$isolation,levels = c("Classic","3DLD"))


seeded=20000
per_well$CFE <- (per_well$Organoids_per_well/seeded)*100

p4<-per_well %>% ggplot(aes(x=isolation, y=CFE, fill=isolation)) + 
                  geom_boxplot(width=0.3) + 
                  geom_jitter(width = 0.3)+ expand_limits(y=0) + 
                  ylab("Colony Formation Effeciency (%)") + 
                  scale_fill_manual(values = c("purple","goldenrod1")) +
                  theme_bw()+xlab("") + 
                  geom_rect(data=rect, aes(xmin=xmin, xmax=xmax, ymin=ymin, ymax=ymax),
                            fill="green",
                            alpha=0.1,
                            inherit.aes = FALSE)+ theme(text = element_text(size=20))

p5 <- p4 + geom_signif(comparisons = list(c("Classic", "3DLD")), map_signif_level = T, y_position = 1.25, textsize = 6)

p5

ggsave(plot = p5, filename = "Fig2C.png", width = 4 , height = 5)



# ********************************************************************* #
# Figure 3D. Organoid diamter
# ********************************************************************* #

organoids$well_id <- paste0(organoids$experiment,organoids$isolation, organoids$Well)
size_per_well <- aggregate(area_um ~ well_id+isolation, organoids, FUN="mean")

organoids$diameter <- 2*sqrt(organoids$area_um / pi)


ppp <- ggplot(organoids, aes(x=isolation, y=diameter, fill= isolation)) + geom_boxplot(width=0.5) + geom_jitter(width=0.25,alpha=0.1) + ylab("Organoid Diameter (µm)")+ scale_fill_manual(values = c("darkorchid2","goldenrod1")) + xlab("Isolation")  +theme_bw() + theme(text = element_text(size=20)) + expand_limits(y=2500)


pppp <- ppp + geom_signif(comparisons = list(c("Classic", "3DLD")), map_signif_level = T, y_position = 2250,textsize = 6 )


ggsave(plot = pppp, filename = "Fig3D.png", width = 4 , height = 5)


# ********************************************************************* #
# Figure 3E. Organoid Area Distributions
# ********************************************************************* #

### Calculate the number of bins to be used:
### USe struge's law: K = 1 + 3. 322 logN

N1<- nrow(organoids[organoids$isolation=="Classic",])
N2<- nrow(organoids[organoids$isolation=="3DLD",])

K1 = 1 + 3.322*log(N1)

K2 = 1 + 3.322*log(N2)

bins=round(mean(K1, K2))

p1 <-ggplot(organoids, aes(x=log(area_um), fill=isolation)) +
  #facet_wrap(~isolation,nrow = 2)+
  geom_histogram(aes(y =..density..),alpha=0.5, position="dodge", bins = bins)+
  geom_density(alpha=0.2)+
  theme(legend.position="top") + scale_fill_manual(values = c("darkorchid2","goldenrod1")) +ylab("Frequency")+ xlab("Organoid Size: log(Area(µm^2))") +theme_bw()+ theme(text = element_text(size=20))
p1


ggsave(plot = p1, filename = "organoid_size_dist_bin.png", width = 8, height = 4)

# ********************************************************************* #
# Figure 3G. Organoid Area Distributions
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
s2c$condition <- relevel(s2c$condition, ref = "pellet_Classic")

dds <- DESeqDataSetFromMatrix(counts, colData = s2c, design = ~condition)
dds <- estimateSizeFactors( dds )
sizeFactors( dds )

dds<-estimateDispersions(dds)
keep <- rowSums(counts(dds))>=1
dds <- dds[keep,]
dds <- DESeq(dds)

### Function DE is defined in Functions.R. MAKE SURE TO RUN ALL CODE IN Functions.R before running this.
results_pell_newVsold <- DE(dds, conditions = c("pellet_Classic", "pellet_3DLD"), xlim = c(-13,13), labSize = 6.0, ylim = c(0,12.5))
results_pell_newVsold$volcano

