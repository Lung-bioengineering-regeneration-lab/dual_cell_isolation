# ********************************************************************* #
# Functions used in publication
# ********************************************************************* #

DE <- function(dds=NA, conditions, coefs=NA, xlim,ylim, labSize){
  require(DESeq2)
  require(EnhancedVolcano)
  if (is.na(dds)){
    print("NO Deseq object inputed")
  }
  conds <- c("condition", conditions)
  res <- results(dds, contrast = conds)
  if(is.na(coefs)==FALSE){
    res <- lfcShrink(coef = coefs, dds = dds)
  }  
  print(summary(res))
  
  res <- res[order(res$pvalue, decreasing = F),]
  write.csv(res, file = paste0(c("R_output/", conditions , ".csv"), collapse = ""))
  
  res <- as.data.frame(res)
  
  volcano <- EnhancedVolcano::EnhancedVolcano(res
                                              ,lab = rownames(res)
                                              , x = "log2FoldChange"
                                              , y = "pvalue"
                                              , pCutoff = 0.01, xlim = xlim, ylim=ylim, labSize = labSize
  )
  
  r <- list("res"=res, "volcano"=volcano)
  
  return(r)
}


plotHeatmap <- function(genes, c, f, data, dataType=NA, fil=NA, labs=c("condition"), h, w, cc=T, cr=T){
  
  gdata <- data[tolower(data$ext_gene) %in% tolower(genes),]
  hmap <- gdata[,2:19]
  annot.row <- as.data.frame(cbind(gdata$ens_gene, gdata$ext_gene, gdata$transcript_type))
  colnames(annot.row) <- c("ens_gene", "ext_gene")
  rownames(annot.row) <- annot.row$ens_gene
  
  s2c <- s2c[s2c$condition %in% f,] 
  hmap <- hmap[, colnames(hmap)%in%rownames(s2c)]
  
  
  if(is.na(dataType)==FALSE){
    hmap <- log(hmap)
    hmap[hmap==-Inf]<-0
    hmap[hmap==Inf]<-0    
  }
  annot.col <- as.data.frame(s2c[,labs])
  rownames(annot.col) <- rownames(s2c)
  colnames(annot.col) <- labs
  print(hmap)
  
  b <- pheatmap(hmap, cluster_rows = T, cluster_cols = cc, legend = cr,
                cellwidth = w, cellheight = h, 
                annotation_col = annot.col,
                annotation_row = annot.row,
                labels_row = gdata$ext_gene,
                angle_col = 0,
                filename = fil
  )
  
}




evaluatePCs <- function(so, nC, lowPC, highPC){
  
  colNum <- highPC-lowPC+1
  flamingo <- matrix(ncol = colNum, nrow = ncol(so))
  flamingo2 <- matrix(ncol = colNum, nrow = ncol(so))
  
  j <- 0
  for(i in lowPC:highPC) {
    j <- j + 1
    blossom <- RunUMAP(so, dims = 1:i, verbose = F)
    flamingo[, j] <- blossom@reductions$umap@cell.embeddings[, 1]
    flamingo2[, j] <- blossom@reductions$umap@cell.embeddings[, 2]
    print(paste((100*j/(highPC - lowPC +1)),"% of UMAP optimization done !!!"))
  }
  
  colnames(flamingo) <- lowPC:highPC
  colnames(flamingo2) <- lowPC:highPC
  
  mrand_obs <- NULL
  flamingo_new <- flamingo[sample(1:ncol(so), nC),]
  flamingo2_new <- flamingo2[sample(1:ncol(so), nC),]
  
  for(i in  1:c(colNum - 1)) {
    blossom <- dist(flamingo_new[, i])
    blossom2 <- dist(flamingo_new[, i + 1])
    blossom3 <- mantel.randtest(blossom, blossom2)$obs
    blossom <- dist(flamingo2_new[, i])
    blossom2 <- dist(flamingo2_new[, i + 1])
    blossom4 <- mantel.randtest(blossom, blossom2)$obs
    mrand_obs <- rbind(mrand_obs, cbind(blossom3, blossom4))
    print(paste((100*i/(colNum -1)),"% of Mantel test done !!!"))
  }
  data.frame("x"=c(lowPC:c(highPC-1)), "Mantel1"=mrand_obs[, 1], "Mantel2"=mrand_obs[, 2]) %>% 
    gather(key="line", value="observation",-x) %>% 
    ggplot(aes(x=x, y=observation, col=line)) + geom_point() + geom_smooth(method = "loess", se = FALSE,span = 0.25, linetype=2, size = 0.5)+ theme_bw() -> mantel_plot
  
  ElbowPlot(so, ndims = 30) -> elbow_plot
  
  return(list("mantel"=mantel_plot, "elbow"=elbow_plot))
}


