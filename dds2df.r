dds2df <- function(res, name){
  df <- as.data.frame(mutate(as.data.frame(res), sig=ifelse(res$padj<0.05, 'FDR<0.05', 'Not Sig')), row.names = rownames(res))
  df <- df[order(df$padj),]
  df$geneid <- rownames(df)
  df$comp <- name
  df$annot <- cr[match(df$geneid, crv3.2), c('cr')]
  df$annot2 <- cr[match(df$geneid, crv3.2), c('cr')]
  df$araport <- cr[match(df$geneid, crv3.2), c('araport')]
  df$araport2 <- cr[match(df$geneid, crv3.2), c('araport2')]
  return(df)
}
