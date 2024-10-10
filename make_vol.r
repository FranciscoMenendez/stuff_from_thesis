make.vol <- function(res, res_name='noname', out=F, 
                     file_name=NULL, file_path=NULL, 
                     file_width = 7, file_height=7,
                     file_units='in', file_dpi=300){
  require(dplyr)
  require(ggplot2)
  df <- as.data.frame(mutate(as.data.frame(res), sig=ifelse(res$padj<0.05, 'FDR<0.05', 'Not Sig')), row.names = rownames(res))
  df <- df[order(df$padj),]
  vol <- print(ggplot(df, aes(log2FoldChange, -log10(padj)), environment = environment()) +
                 geom_point(aes(col = sig, alpha=0.5), size=0.1) +
                 scale_color_manual(values = c("#18A558", "#050A30")) +
                 ggtitle(res_name)+
                 ggrepel::geom_text_repel(data=df[1:10, ], aes(label=rownames(df)[1:10]), size=3)+
                 scale_alpha(guide='none')+
                 theme_minimal()+
                 #ylim(0,51)+
                 theme(legend.key.size = unit(1, "cm"), legend.title = element_blank(),legend.position = 'none'))
  
  if (out==F) {
    return(vol) 
  } else{
    ggsave(filename = file_name,
           path=file_path,
           plot = vol,
           device = 'png',
           width = file_width,
           height = file_height,
           units = file_units,
           dpi = file_dpi)
    
  }
}
