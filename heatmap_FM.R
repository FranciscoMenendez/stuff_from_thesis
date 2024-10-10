make_zscr_heatmap <- function(data = counts(dd_leaf.le),
                              file_name, 
                              ngene=150, 
                              dendo='row', 
                              row_ordered=T,
                              colsidecolors=NULL,
                              colsp=NULL, 
                              rowsp=NULL, 
                              named_genes=NULL,
                              brks=NULL, 
                              named_cat=gene_cat,
                              ncat=3, 
                              HEIGHT=26,
                              WIDTH=12,
                              MARGINS = c(14, 12)) {
  
  # Set the plot dimensions.
  require(RColorBrewer)
  # Command line argument.
  args = commandArgs(trailingOnly=TRUE)  
  
  # Optional argument for FDR cutoff.
  if (length(args) > 0) {
    LIMIT = as.integer(args[1])
  } else {
    # Default FDR cutoff
    LIMIT = 5
  }
  
  # Input values should be in percent!
  LIMIT = LIMIT/100
  
  # Set the margins
  
  
  # Relative heights of the rows in the plot.
  LHEI = c(.45, 5)
  
  # Load the library.
  suppressPackageStartupMessages(require(gplots))
  
  # Read normalized counts from the standard input.
  #data = read.csv('stdin', header=T, as.is=TRUE)
  print(1)
  # Subset data for values under a threshold.
  ##data = subset(data, data$FDR <= LIMIT)
  print(2)
  # The heatmap row names will be taken from the first column.
  #row_names = data[, 1]
  
  # The normalized data starts past the rightmost of these columns.
  #idx = which(colnames(data) == "falsePos") + 1
  
  # The normalized counts are on the right size.
  #counts = data[, idx : ncol(data)]
  
  # Load the data from the second column on.
  values = as.matrix(data)
  
  # Adds a little noise to each element to avoid the
  # clustering function failure on zero variance rows.
  values = jitter(values, factor = 1, amount = 0.00001)
  print(3)
  # Normalize each row to a z-score
  #zscores = apply(X = values, MARGIN = 1, FUN = function(row) (row - mean(row)) / sd(row), simplify = T)
  zscores = t(scale(t(values)))
  #zscores2 = NULL
  
  #for (i in 1 : nrow(values)) {
  #  row = values[i,]
  #  zrow = (row - mean(row)) / sd(row)
  #  zscores2 = rbind(zscores2, zrow)
  #}
    
  # Set the row names on the zscores.

  #row.names(zscores2) <- rownames(data)
  
  # Turn the data into a matrix for heatmap2.
  zscores = as.matrix(zscores)
  zscores <- zscores[order(rowVars(zscores)), ]
  print(4)
  # Open the drawing device.
  pdf(file_name, width = WIDTH, height = HEIGHT)
  
  # Set the color palette.
  col = greenred
  print(5)
  print(sum(rownames(zscores) %in% named_genes))
  df_mat <- zscores[rownames(zscores) %in% named_genes,]  #matrix with named_genes
  f <- named_cat[!duplicated(named_cat$gene_id), ]
  
  if (sum(duplicated(named_cat$gene_id)!=0)){
    dup <- named_cat[duplicated(named_cat$gene_id), 'gene_id']
    mypal <- brewer.pal((ncat+1), "Dark2")
    f[f$gene_id %in% dup, 'Category'] <- 'Multiple'
  } else{mypal <- brewer.pal(ncat, "Dark2")}#vector with duplicated genes in named_genes
  
  f <- factor(f$Category)
  print(length(f))
  print(dim(df_mat))
  
  # Draw the heatmap.
  set.seed(3)
  if (is.null(named_genes)){
    heatmap.2(zscores[1:ngene,], col=col, density.info="none", Colv=NULL, colsep = colsp, cexCol = 20, cexRow = 20,
              dendrogram=dendo, trace="none", margins=MARGINS, lhei=LHEI, srtCol = 45)  
  } else {
    

    heatmap.2(df_mat, col=col, density.info="none", Colv=NULL, colsep = colsp, Rowv = row_ordered,
              dendrogram=dendo, trace="none", margins=MARGINS, lhei=LHEI, srtCol = 45, rowsep = rowsp, breaks = brks, RowSideColors = mypal[f])
    legend(x=c('top'), ncol=2, legend=unique(f), col=mypal[], pch=15)
    
  }
  # Turn off the device.
  dev.off()
  
}


#!/usr/bin/env Rscript
#
# Draws a heatmap from the output that contains a normalized matrix
#

#make_zscr_heatmap(data = counts(dd_leaf.le), file_name = 'leaf_libraries_top150.pdf', ngene = 150, colsp = c(3, 6, 9, 12, 15, 18))
#make_zscr_heatmap(data = counts(dd_root_nl.le), file_name = 'root_libraries_top150.pdf', colsp = c(3, 6, 9))
