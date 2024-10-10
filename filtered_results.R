filtered_results <- function(dds=dds_leaf.le, lib_index=1:6, name_contrast, contr=T){
  dds <- estimateSizeFactors(dds)
  idx <- rowSums( counts(dds, normalized=T)[,lib_index] > 0 ) >=2
  table(idx)
  dds <- dds[idx]
  des <- DESeq(dds, test = 'Wald', betaPrior = F, parallel = T)
  if (contr==T){
    res <- results(des, contrast = name_contrast,
                   filter = rowVars(counts(des, normalized=T)[, lib_index]),
                   theta = seq(0, .80, .01))}
  else {
    res <- results(des, name = name_contrast,
                   filter = rowVars(counts(des, normalized=T)[, lib_index]),
                   theta = seq(0, .80, .01))
  }
  return(res)
  
}