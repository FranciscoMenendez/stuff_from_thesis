# Make GSEA enrichement and print curves 

gsea_enrich <- function(ranked_list_file='leaf_24h_0h_gsea_list.txt', gmt_file = sl, seed=23, save2file=F, byNES=T) {
  require('fgsea')
  require(ggplot2)
  require(gggsea)
  require(dplyr)
  #print(ranked)
  
  print(1)
  
  # Read and prep ranked list
  #print(ranked_list_file)
  ranked_list <- read.table(ranked_list_file)
  rl <- ranked_list$V2
  names(rl) <- ranked_list$V1
  print(head(rl))
  print(2)
  #calculate enrichement
  set.seed(seed)
  
  gsea_test <- fgsea::fgseaMultilevel(sl,  rl)
  print(3)
  if(byNES==T) {
    gsea_test <- gsea_test[order(abs(gsea_test$NES), decreasing = T), ]
    gsea_test <- gsea_test[gsea_test$padj<0.05,]
    
  }else{
    gsea_test <- gsea_test[order(abs(gsea_test$padj), decreasing = T), ]
    gsea_test <- gsea_test[gsea_test$padj<0.05,]
  }
  gsea_test$comparison <- ranked_list_file
  
  print(4)
  #make curve DF and print to file
  
  print(5)
  
  if (save2file==T) {
    curve_gsea <- gseaCurve(rl = rl, setlist = sl[gsea_test$pathway], gsea = gsea_test)
    curve_gsea$set <- factor(curve_gsea$set, unique(curve_gsea$set), labels=go_terms[go_terms$ids %in% unique(curve_gsea$set), "go_names"])
    gcurve <- ggplot()+gggsea::geom_gsea(curve_gsea)
    
    ggsave(filename = paste(ranked_list_file, '_curves', sep = ''), plot = gcurve, device = 'pdf', path = '.', dpi = 300)
    
  }else{
    return(gsea_test)
  }
  
}

setwd('res_root/GSEA_lists/')
test1 <- gsea_enrich('root_pa121_tsh660_0h_0h_gsea_list.txt')
test2 <- gsea_enrich('root_pa121_tsh660_48h_48h_gsea_list.txt')

setwd('../../res_leaf/GSEA_lists/')
test3 <- gsea_enrich('leaf_pa121_tsh660_0h_0h_gsea_list.txt')
test4 <- gsea_enrich('leaf_pa121_tsh660_24h_24h_gsea_list.txt')
test5 <- gsea_enrich('leaf_pa121_tsh660_48h_48h_gsea_list.txt')

go_terms[grep('GO:0080022', go_terms$ids),]

go_terms[go_terms$ids %in% unique(gcurve_root_pa121_48_0$set), "go_names"]

go_terms[go_terms$ids %in% unique(gcurve_root_pa121_48_0$set), ]

df_fg <- gseaCurve(rl = ranklist, setlist = sl[c("GO:0048046", "GO:0005507", "GO:0009505")], gsea = gsea_test ) #make curve with gggsea
#fgsea::plotEnrichment(pathway = sl[["GO:0005507"]], ranked_root_48)

df_fg$set <- factor(df_fg$set, unique(df_fg$set), labels=c('Apoplast', 'Copper ion binding', 'Plant-type cell wall'))
gsea_root_48h <- ggplot()+
  gggsea::geom_gsea(df_fg)
ggsave(filename = 'enrichment_score_gsea_root_48_48h', plot = gsea_root_48h, device = 'pdf', width = 21, height = 18, units = 'cm', dpi = 300)
