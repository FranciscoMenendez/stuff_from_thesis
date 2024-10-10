#res=res_leaf_geno$leaf_pa121_cd_no_cd
#des=dd_leaf.le
#index=index_leaf_geno$leaf_pa121_cd_no_cd





baseMean_add <- function(res, des=dd_leaf.le, index, uneq=F, index_a, index_b){
  des=estimateSizeFactors(des)
  res$baseMeanA <- 1
  res$baseMeanB <- 1
  
  print(1)
  
  normed = counts(des, normalized=TRUE) 
  
  print(rownames(res)[1:10])
  print(rownames(normed)[1:10])
  
  normed <- normed[rownames(normed) %in% rownames(res),]
  
  print(dim(normed))
  
  print(2)
  
  if(uneq==T){
    
    print(colnames(normed[,index_a]))
    
    res$baseMeanA = rowMeans(normed[,index_a])
    
    print(colnames(normed[,index_b]))
    
    res$baseMeanB = rowMeans(normed[,index_b])
    
    print(3)
    
    return(res)
    
  } else if (length(index) %% 2 == 0) {
    
    print(4)
    
    print(colnames(normed[, 1:(length(index)/2)]))
    
    res$baseMeanA <- rowMeans(normed[, 1:(length(index)/2)])
    
    print(5)
    
    print(colnames(normed[, (length(index)/2+1):length(index)]))
    
    res$baseMeanB <- rowMeans(normed[, (length(index)/2+1):length(index)])
    
    print(6)
    
    return(res)
    
  }
}


