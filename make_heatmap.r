make.heatmap <- function(dds, lib_index=1:6, title) {
  require(pheatmap)
  require(RColorBrewer)
  require(DESeq2)
  dds <- estimateSizeFactors(dds)
  idx <- rowSums( counts(dds, normalized=T)[,lib_index] > 0 ) >=2
  table(idx)
  dds <- dds[idx]
  vsd <- vst(dds) 
  sampleDist <- dist(t(assay(vsd)))
  sampleDistMatrix <- as.matrix(sampleDist)
  rownames(sampleDistMatrix) <- paste(vsd$genotype, vsd$time, vsd$block, vsd$tissue, sep="_")
  colnames(sampleDistMatrix) <- paste(vsd$genotype, vsd$time, vsd$block, vsd$tissue, sep="_")
  colors <- brewer.pal(9, 'Spectral')
  pheatmap(sampleDistMatrix, main=title, cluster_cols = F,
           clustering_distance_rows=sampleDist,
           clustering_distance_cols=sampleDist,
           col=colors)
  
}