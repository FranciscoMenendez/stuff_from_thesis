
hm_data2 <- function(gene_ids){ #heatmap data
  
  # Root & LEAF DATA
  
  #CREATE EMPTY DATA TABLES
  
  leaf_48h_0h_l2f <- data.table()
  leaf_tsh660_48h_0h_l2f <- data.table()
  leaf_pa121_48h_0h_l2f <- data.table()
  
  leaf_48h_0h_padj <- data.table()
  leaf_tsh660_48h_0h_padj <- data.table()
  leaf_pa121_48h_0h_padj <- data.table()
  
  root_48h_0h_l2f <- data.table()
  root_tsh660_48h_0h_l2f <- data.table()
  root_pa121_48h_0h_l2f <- data.table()
  
  root_48h_0h_padj <- data.table()
  root_tsh660_48h_0h_padj <- data.table()
  root_pa121_48h_0h_padj <- data.table()
  
  # LOOP AROSS GENES 
  for (x in 1:length(gene_ids)) {
    print(x)
    # Extract l2f from each data table
    root_48h_0h_l2f <- rbind(root_48h_0h_l2f, df_res$root_48h_0h[grep(gene_ids[x], geneid), .(orthogroup, geneid, log2FoldChange), with = TRUE])
    root_tsh660_48h_0h_l2f <- rbind(root_tsh660_48h_0h_l2f, df_res$root_tsh660_48h_0h[grep(gene_ids[x], geneid), .(geneid, log2FoldChange), with = TRUE])
    root_pa121_48h_0h_l2f <- rbind(root_pa121_48h_0h_l2f, df_res$root_pa121_48h_0h[grep(gene_ids[x], geneid), .(geneid, log2FoldChange), with = TRUE])
    
    
    leaf_48h_0h_l2f <- rbind(leaf_48h_0h_l2f, df_res$leaf_48h_0h[grep(gene_ids[x], geneid), .(orthogroup, geneid, log2FoldChange), with = TRUE])
    leaf_tsh660_48h_0h_l2f <- rbind(leaf_tsh660_48h_0h_l2f, df_res$leaf_tsh660_48h_0h[grep(gene_ids[x], geneid), .(geneid, log2FoldChange), with = TRUE])
    leaf_pa121_48h_0h_l2f <- rbind(leaf_pa121_48h_0h_l2f, df_res$leaf_pa121_48h_0h[grep(gene_ids[x], geneid), .(geneid, log2FoldChange), with = TRUE])
    
    # Extract padj from each data table
    root_48h_0h_padj <- rbind(root_48h_0h_padj, df_res$root_48h_0h[grep(gene_ids[x], geneid), .(orthogroup, geneid, padj), with = TRUE])
    root_tsh660_48h_0h_padj <- rbind(root_tsh660_48h_0h_padj, df_res$root_tsh660_48h_0h[grep(gene_ids[x], geneid), .(geneid, padj), with = TRUE])
    root_pa121_48h_0h_padj <- rbind(root_pa121_48h_0h_padj, df_res$root_pa121_48h_0h[grep(gene_ids[x], geneid), .(geneid, padj), with = TRUE])
    
    leaf_48h_0h_padj <- rbind(leaf_48h_0h_padj, df_res$leaf_48h_0h[grep(gene_ids[x], geneid), .(orthogroup, geneid, padj), with = TRUE])
    leaf_tsh660_48h_0h_padj <- rbind(leaf_tsh660_48h_0h_padj, df_res$leaf_tsh660_48h_0h[grep(gene_ids[x], geneid), .(geneid, padj), with = TRUE])
    leaf_pa121_48h_0h_padj <- rbind(leaf_pa121_48h_0h_padj, df_res$leaf_pa121_48h_0h[grep(gene_ids[x], geneid), .(geneid, padj), with = TRUE])
    
  }
  
  #MERGE DATA TABLES, WITH = ALL
  leaf_l2f <- merge(x = leaf_48h_0h_l2f, y = leaf_tsh660_48h_0h_l2f, by.x = 'geneid', by.y='geneid', all = TRUE)
  leaf_l2f <- merge(x = leaf_l2f, y = leaf_pa121_48h_0h_l2f, by.x = 'geneid', by.y = 'geneid', all = TRUE)
  
  leaf_padj <- merge(leaf_48h_0h_padj, leaf_tsh660_48h_0h_padj, by.x = 'geneid', by.y='geneid', all = TRUE)
  leaf_padj <- merge(leaf_padj, leaf_pa121_48h_0h_padj, by.x = 'geneid', by.y='geneid', all = TRUE)
  
  root_l2f <- merge(root_48h_0h_l2f, root_tsh660_48h_0h_l2f, by = 'geneid', all = TRUE)
  root_l2f <- merge(root_l2f, root_pa121_48h_0h_l2f, by.x = 'geneid', by.y='geneid', all = TRUE)
  
  root_padj <- merge(root_48h_0h_padj, root_tsh660_48h_0h_padj, by = 'geneid', all = TRUE)
  root_padj <- merge(root_padj, root_pa121_48h_0h_padj, by.x = 'geneid', by.y='geneid', all = TRUE)
  
  #ORDER BY ORTHOGROUP 
  
  leaf_l2f <- leaf_l2f[order(orthogroup)]
  leaf_padj <- leaf_padj[match(leaf_l2f$geneid, leaf_padj$geneid),]
  colnames(leaf_l2f) <- c('geneid', 'gfam', 'BOTH', "TSH660", 'PA121')
  colnames(leaf_padj) <- c('geneid', 'gfam', 'BOTH', "TSH660", 'PA121')
  
  root_l2f <- root_l2f[order(orthogroup)]
  root_padj <- root_padj[match(root_l2f$geneid, root_padj$geneid), ]
  colnames(root_l2f) <- c('geneid', 'gfam', 'BOTH', "TSH660", 'PA121')
  colnames(root_padj) <- c('geneid', 'gfam', 'BOTH', "TSH660", 'PA121')
  
  
  
  if(length(!is.na(root_l2f$geneid))>length(!is.na(leaf_l2f$geneid))){
    leaf_l2f <- leaf_l2f[match(root_l2f[,geneid],  leaf_l2f[, geneid]),]
    leaf_l2f$geneid <- root_l2f$geneid
    leaf_padj <- leaf_padj[match(root_padj[,geneid],  leaf_padj[, geneid]),]
    leaf_padj$geneid <- root_padj$geneid
  } else if(length(!is.na(root_l2f$geneid))<length(!is.na(leaf_l2f$geneid))){
    root_l2f <- root_l2f[match(leaf_l2f[,geneid],  root_l2f[, geneid]),]
    root_l2f$geneid <- leaf_l2f$geneid
    root_padj <- root_padj[match(leaf_padj[,geneid],  root_padj[, geneid]),]
    root_padj$geneid <- leaf_padj$geneid
  }
  
  lis <- list(root_l2f, root_padj, leaf_l2f, leaf_padj)
  names(lis) <- c('Rtl2f', 'Rtpadj', 'Lfl2f', 'Lfpadj')
  
  return(lis)
}


add_labels_to_factor <- function(data_table, column_name, labels, level_order = NULL) {
  unique_values <- unique(data_table[[column_name]])
  unique_values <- unique_values[!is.na(unique_values)]
  factor_column <- factor(data_table[[column_name]], levels = unique_values)
  
  # Relabel the factor levels using the new labels
  levels(factor_column) <- labels[match(levels(factor_column), unique_values)]
  
  # Reorder the factor levels if level_order is provided
  if (!is.null(level_order)) {
    factor_column <- factor(factor_column, levels = level_order)
  }
  
  data_table[[column_name]] <- factor_column
  
  return(data_table)
}


load('df_res.rda')
load('cr.rda')
library(data.table)

# just the heatmap

# 1. Choose gene families,

# 2. Make hheatmap data hm_data2
# 3. make hm_data list 
# 4. Name gene families (Add_gene_labels)
# 5. Make Heatmap (only sig )

#define the metal transporters

hma_ids <- c('Tc01v2_t005580',
             'Tc01v2_t018910',
             'Tc02v2_t000500',
             'Tc02v2_t011530',
             'Tc03v2_t023170', 
             'Tc09v2_t019080',
             'Tc09v2_t015490',
             'Tc09v2_t015500')

hma_tc <- unique(cr[unlist(lapply(hma_ids, grep, lookup)),crv3.2])

hma_hd <- hm_data2(gene_ids = unique(cr[unlist(lapply(hma_ids, grep, lookup)),crv3.2]))

hma_hd <- lapply(hma_hd, add_labels_to_factor, column_name='gfam', labels=c('HMA'))

hma_hd <- lapply(hma_hd, function(x, column_name, label) {
  x[[column_name]] <- label  # Assign the label to the specified column
  return(x)  # Return the modified data table
}, column_name = 'gfam', label = 'HMA')  # Specify the column_name and label arguments

cr[ncbi %in% c('XM_007047334.2', 'XM_018116084.1', 'XM_007040138.2', 'XM_018127237.1'), geneid]

hma_hd <-  hm_data2(gene_ids = unique(cr[ncbi %in% c('rna-XM_007047334.2', 'rna-XM_018116084.1', 'rna-XM_007040138.2', 'rna-XM_018127237.1'), crv3.2]))

hma_hd[[1]][, gfam:=as.character(gfam)] 
hma_hd[[1]][geneid %in% c('Tc02cons_g000500', 'Tc03cons_g026810'), gfam:='HMA5'] 

hma_hd[[2]][, gfam:=as.character(gfam)] 
hma_hd[[2]][geneid %in% c('Tc02cons_g000500', 'Tc03cons_g026810'), gfam:='HMA5'] 

hma_hd[[3]][, gfam:=as.character(gfam)] 
hma_hd[[3]][geneid %in% c('Tc01cons_g005630'),  gfam:='HMA1'] 

hma_hd[[4]][, gfam:=as.character(gfam)] 
hma_hd[[4]][geneid %in% c('Tc01cons_g005630'), gfam:='HMA1'] 

hma_hd[[3]][, gfam:=as.character(gfam)] 
hma_hd[[3]][geneid %in% c('Tc09cons_g016910'),  gfam:='HMA3'] 

hma_hd[[4]][, gfam:=as.character(gfam)] 
hma_hd[[4]][geneid %in% c('Tc09cons_g016910'), gfam:='HMA3'] 


keep <- !is.na(apply(hma_hd$Rtpadj[, -1], 1, function(x) any(x < 0.05)))
keep <- !is.na(apply(hma_hd$Lfpadj[, -1], 1, function(x) any(x < 0.05)))


# ZIP
zip_ids <- readLines('zachs_irt')
zip_ids <- zip_ids[2:length(zip_ids)]

zip_ids <- unlist(lapply(zip_ids,  substr, 1, 16))
zip_ids <- unique(zip_ids)

zip_hd <- hm_data2(gene_ids = zip_ids)
zip_hd <- lapply(zip_hd, add_labels_to_factor, column_name='gfam', labels='ZIP') 

zip_hd <- lapply(zip_hd, function(x, column_name, label) {
  x[[column_name]] <- label  # Assign the label to the specified column
  return(x)  # Return the modified data table
}, column_name = 'gfam', label = 'ZIP')  # Specify the column_name and label arguments

zip_hd

lis <- list(hma_hd, zip_hd, vit_hd)

hm_l2(lis=test_hm_lst, )

#first only HMA and ZIP



# VIT
vit_ids <- df_res$root_48h_0h[grep('vacuolar iron', tolower(annot)),.(geneid, annot, log2FoldChange, padj)][,geneid]
df_res$root_pa121_48h_0h[grep('vacuolar iron', tolower(annot)),.(geneid, annot, log2FoldChange, padj)][,geneid]
vit_hd <- hm_data2(gene_ids = vit_ids)
vit_hd <- lapply(vit_hd, function(x) x[, gfam:=as.character(gfam)])
vit_hd <- lapply(vit_hd, function(x) x[gfam %in% 662, gfam:='VIT4'])
vit_hd <- lapply(vit_hd, function(x) x[gfam %in% 3820, gfam:='VIT1'])


# PCR
pcr_ids <- df_res$root_48h_0h[grep('plant cadmium', tolower(annot)),.(geneid, annot, log2FoldChange, padj)][,geneid]
pcr_hd <- hm_data2(gene_ids = pcr_ids)
pcr_hd[c(1,2)] <- lapply(pcr_hd[c(1,2)], add_labels_to_factor, column_name='gfam', labels=c('PCR', 'PCR2'))
pcr_hd[c(3,4)] <- lapply(pcr_hd[c(3,4)], add_labels_to_factor, column_name='gfam', labels=c('PCR', 'PCR2'))

#NRAMP
nramp_ids <- df_res$root_48h_0h[grep('nramp', tolower(annot)),.(geneid, annot, log2FoldChange, padj)][,geneid]
nramp_hd <- hm_data2(gene_ids = nramp_ids)
nramp_hd[c(1,2)] <- lapply(nramp_hd[c(1,2)], add_labels_to_factor, column_name='gfam', labels='NRAMP')


# MATE
mate_ids <- df_res$root_48h_0h[grep('protein detoxification', tolower(annot)),.(geneid, annot, araport, annot2, log2FoldChange, padj)][padj<0.05,geneid]
mate_ids <- append(mate_ids, df_res$root_pa121_48h_0h[grep('protein detoxification', tolower(annot)),.(geneid, annot, araport, annot2, log2FoldChange, padj)][padj<0.05,geneid])
mate_ids <- append(mate_ids, df_res$root_tsh660_48h_0h[grep('protein detoxification', tolower(annot)),.(geneid, annot, araport, annot2, log2FoldChange, padj)][padj<0.05,geneid])
mate_ids <- append(mate_ids, df_res$leaf_48h_0h[grep('protein detoxification', tolower(annot)),.(geneid, annot, araport, annot2, log2FoldChange, padj)][padj<0.05,geneid])
mate_ids <- append(mate_ids, df_res$leaf_pa121_48h_0h[grep('protein detoxification', tolower(annot)),.(geneid, annot, araport, annot2, log2FoldChange, padj)][padj<0.05,geneid])
mate_ids <- append(mate_ids, df_res$leaf_tsh660_48h_0h[grep('protein detoxification', tolower(annot)),.(geneid, annot, araport, annot2, log2FoldChange, padj)][padj<0.05,geneid])
mate_ids <- unique(mate_ids)

mate_hd <- hm_data2(gene_ids = mate_ids)
mate_hd[c(1,2)] <- lapply(mate_hd[c(1,2)], add_labels_to_factor, column_name='gfam', labels=c('MATE1', 'MATE2', 'MATE3', 'MATE4', 'MATE5', 'MATE6', 'MATE7', 'MATE8'))

mate_hd[c(3:4)] <- lapply(mate_hd[c(3:4)], add_labels_to_factor, column_name='gfam', labels=c('MATE9', 'MATE10'))

# Magnesium 
mrs2_ids <- df_res$root_48h_0h[grep('mrs', tolower(annot)),.(geneid, annot, araport, annot2, log2FoldChange, padj)][,geneid]
mrs2_hd <- hm_data2(gene_ids = mrs2_ids)
mrs2_hd[c(1,2)] <- lapply(mrs2_hd[c(1,2)], add_labels_to_factor, column_name='gfam', labels=c('MRS2a', 'MRS2b'))

# CorA Magnesium transporter
cora_ids <- df_res$root_48h_0h[grep('cora', tolower(annot)),.(geneid, annot, araport, annot2, log2FoldChange, padj)][,geneid]
cora_hd <- hm_data2(gene_ids = cora_ids)
cora_hd[c(1,2)] <- lapply(cora_hd[c(1,2)], add_labels_to_factor, column_name='gfam', labels=c('CORA1', 'CORA2'))

# Metallochaperones
# ATX1
atx1_ids <- df_res$root_48h_0h[grep('atx', tolower(annot)),.(geneid, annot, araport, annot2, log2FoldChange, padj)][,geneid]
atx1_hd <- hm_data2(gene_ids = atx1_ids)
atx1_hd[c(1,2)] <- lapply(atx1_hd[c(1,2)], add_labels_to_factor, column_name='gfam', labels=c('ATXa', 'ATXb'))
# Heavy metal associated
hmp_ids <- df_res$root_48h_0h[grep('heavy metal transport', tolower(annot)),.(geneid, annot, araport, annot2, log2FoldChange, padj)][,geneid]
hmp_hd <- hm_data2(gene_ids = hmp_ids)
hmp_hd[c(1,2)] <- lapply(hmp_hd[c(1,2)], add_labels_to_factor, column_name='gfam', labels=c('HMP1', 'HMP2'))

metallochap_og <- unique(df_res$root_48h_0h[geneid %in% c(hmp_ids, atx1_ids)][padj<0.05, orthogroup])
metch_id <- df_res$root_48h_0h[orthogroup %in% metallochap_og,.(geneid, annot, araport, annot2, log2FoldChange, padj)][,geneid]
metch_hd <- hm_data2(gene_ids = metch_id)
metch_hd[1:2] <- lapply(metch_hd[1:2], add_labels_to_factor, column_name='gfam', labels=c('HIPP', 'CCH', 'HMP'))

# Copt Copper transporter
copt2_ids <- df_res$root_48h_0h[grep('copper transporter ', tolower(annot)),.(geneid, annot, araport, annot2, log2FoldChange, padj)][,geneid]
 df_res$leaf_48h_0h[grep('copper transporter ', tolower(annot)),.(geneid, annot, araport, annot2, log2FoldChange, padj)][,geneid]
copt2_hd <- hm_data2(gene_ids = copt2_ids)
copt2_hd[c(1,2)] <- lapply(copt2_hd[c(1,2)], add_labels_to_factor, column_name='gfam', labels='COPT')

# YSL
ysl_ids <- df_res$root_48h_0h[grep('metal-nicotianamine', tolower(annot)),.(geneid, annot, araport, annot2, log2FoldChange, padj)][,geneid]
ysl_hd <- hm_data2(gene_ids = ysl_ids)
ysl_hd <- lapply(ysl_hd, add_labels_to_factor, column_name='gfam', labels='YSL')

#CAX
cax_ids <- df_res$root_48h_0h[grep('vacuolar cation', tolower(annot)),.(geneid, annot, araport, annot2, log2FoldChange, padj)][,geneid]
cax_hd <- hm_data2(gene_ids = cax_ids)
cax_hd <- lapply(cax_hd, function(x) x[,gfam := as.character(gfam)])
cax_hd <- lapply(cax_hd, function(x) x[gfam %in% 1980, gfam:='CAX1'])
cax_hd <- lapply(cax_hd, function(x) x[gfam %in% 2335, gfam:='CAX2'])


#ABC
abc_ids <- df_res$root_48h_0h[grep('abc transporter', tolower(annot)),.(sig, geneid, annot, araport, araport2,  annot2, log2FoldChange, padj)][sig=='FDR<0.05',geneid]
abc_hd <- hm_data2(gene_ids = abc_ids)

abc_hd$Rtl2f[,gfam:=as.character(gfam)]

abc_hd$Rtl2f[gfam %in% '50', gfam:='PDR']
abc_hd$Rtl2f[gfam %in% '71', gfam:='ABCC']
abc_hd$Rtl2f[gfam %in% '92', gfam:='ABCB']
abc_hd$Rtl2f[gfam %in% '177', gfam:='ABCp']
abc_hd$Rtl2f[gfam %in% '374', gfam:='ABCG2']
abc_hd$Rtl2f[gfam %in% '748', gfam:='ABCG3']
abc_hd$Rtl2f[gfam %in% '838', gfam:='ABCG4']
abc_hd$Rtl2f[gfam %in% '1285', gfam:='ABCA7']
abc_hd$Rtl2f[gfam %in% '1314', gfam:='ABCG5']
abc_hd$Rtl2f[gfam %in% '1614', gfam:='ABCG24']
abc_hd$Rtl2f[gfam %in% '3409', gfam:='ABCA2']
abc_hd$Rtl2f[gfam %in% '3487', gfam:='ABCD']
abc_hd$Rtl2f[gfam %in% '3625', gfam:='ABCF']
abc_hd$Rtl2f[gfam %in% '3851', gfam:='ABCB25']
abc_hd$Rtl2f[gfam %in% '6778', gfam:='ABCG3']
abc_hd$Rtl2f[gfam %in% '7835', gfam:='ABCA1']
abc_hd$Rtl2f[gfam %in% '8093', gfam:='ABCG7']

abc_hd$Rtpadj[,gfam:=as.character(gfam)]
abc_hd$Rtpadj[gfam %in% '50', gfam:='PDR']
abc_hd$Rtpadj[gfam %in% '71', gfam:='ABCC']
abc_hd$Rtpadj[gfam %in% '92', gfam:='ABCB']
abc_hd$Rtpadj[gfam %in% '177', gfam:='ABCp']
abc_hd$Rtpadj[gfam %in% '374', gfam:='ABCG2']
abc_hd$Rtpadj[gfam %in% '748', gfam:='ABCG3']
abc_hd$Rtpadj[gfam %in% '838', gfam:='ABCG4']
abc_hd$Rtpadj[gfam %in% '1285', gfam:='ABCA7']
abc_hd$Rtpadj[gfam %in% '1314', gfam:='ABCG5']
abc_hd$Rtpadj[gfam %in% '1614', gfam:='ABCG24']
abc_hd$Rtpadj[gfam %in% '3409', gfam:='ABCA2']
abc_hd$Rtpadj[gfam %in% '3487', gfam:='ABCD']
abc_hd$Rtpadj[gfam %in% '3625', gfam:='ABCF']
abc_hd$Rtpadj[gfam %in% '3851', gfam:='ABCB25']
abc_hd$Rtpadj[gfam %in% '6778', gfam:='ABCG3']
abc_hd$Rtpadj[gfam %in% '7835', gfam:='ABCA1']
abc_hd$Rtpadj[gfam %in% '8093', gfam:='ABCG7']


abc_hd$Lfpadj[,gfam:=as.character(gfam)]
abc_hd$Lfpadj[gfam %in% '748', gfam:='ABCG3']
abc_hd$Lfpadj[gfam %in% '3487', gfam:='ABCD']


abc_hd$Lfl2f[,gfam:=as.character(gfam)]
abc_hd$Lfl2f[gfam %in% '748', gfam:='ABCG3']
abc_hd$Lfl2f[gfam %in% '3487', gfam:='ABCD']



#print fasta peptides for lcalization prediction

# read fasta CDS of the b97
library(ape)
b97cds <- read.FASTA('TC_B97_consensusV2_rev2021.cds.fna')

        lis=list(hma_hd,zip_hd,vit_hd,pcr_hd,nramp_hd,hmp_hd, atx1_hd,mate_hd,mrs2_hd,cora_hd,copt2_hd,ysl_hd,cax_hd,abc_hd)
names(lis) <- c('hma',  'zip','vit',  'pcr','nramp', 'hmp',  'atx1',  'mate','mrs2', 'cora', 'copt2', 'ysl',  'cax','abc')

#write.FASTA(trans(b97cds[paste(gsub('g', 't', x = l2f$geneid), '.1', sep='')]), file = 'metal_transporters2.fa')
write.FASTA(trans(b97cds[paste(gsub('g', 't', x = l2f$geneid), '.1', sep='')]), file = 'metal_transporters3.fa')

#read data from TMHHM2 

#tmh <- read.csv('probability_summaries_metal_transporters.csv')
tmh <- read.csv('probability_summaries_metal_transporters2.csv')

tmh <- data.table(tmh)
tmh[, geneid:= substr(gsub(pattern = 't', replacement = 'g', x= Protein_ID), 1, 16)]
tmh[, annot:= df_res$root_48h_0h[geneid %in% tmh$geneid,annot]]
setkey(tmh, geneid)


library(stringr)
# Define the mapping of localizations to abbreviations
shorthand_mapping <- c(
  "Cytoplasm" = "C",
  "Nucleus" = "N",
  "Cell membrane" = "CM",
  "Lysosome/Vacuole" = "LV",
  "Golgi apparatus" = "GA",
  "Plastid" = "P"
)


#Function to apply shorthand mapping
get_abbreviations <- function(localizations, mapping) {
  sapply(str_split(localizations, "\\|"), function(loc_list) {
    paste(sapply(loc_list, function(loc) mapping[loc]), collapse = "|")
  }, USE.NAMES = FALSE)
}

# Apply the mapping to create the 'abbreviations' column
tmh[, abbreviations := get_abbreviations(Localizations, shorthand_mapping)] #these have the same order as rowns (it was printed after...)




row_gfam_ha <- rowAnnotation(family=transporters_l2f$gfam, col = list(family=c("ZIP" = "#8DD3C7",   "HMA" = "#FFFFB3",   "MTP" = "#BEBADA",
                                                                               "VIT" = "#FB8072", "PCR" = "#80B1D3", "NRAMP" = "#FDB462", 
                                                                               "HMP" = "#B3DE69", "MATE" = "#FCCDE5", "MRS2" = "#D9D9D9", 
                                                                               "CorA" = "#BC80BD", "ATX" = "#CCEBC5", "COPT" = "#FFED6F", 
                                                                               "YSL" = "#E41A1C", "CAX" = "#377EB8", "ZIP" = "#4DAF4A", 'ABC'="#0EAD9A")))



#lis=list(hma_hd,zip_hd,vit_hd,pcr_hd,nramp_hd,metch_hd,mate_hd,mrs2_hd,cora_hd,copt2_hd,ysl_hd,cax_hd)

metal_tr_hm <- hm_l2(lis = list(hma_hd,zip_hd,vit_hd,pcr_hd,nramp_hd,hmp_hd,mate_hd,mrs2_hd,cora_hd,atx1_hd,copt2_hd,ysl_hd,cax_hd,abc_hd), onlysig = T,
                     ncbi = T, tissue = 'root')



# just the metal transporter figure 

#colors 
#define some colors

c25 <- c(
  "dodgerblue2", "#E31A1C", # red
  "green4",
  "#6A3D9A", # purple
  "#FF7F00", # orange
  "black", "gold1",
  "skyblue2", "#FB9A99", # lt pink
  "palegreen2",
  "#CAB2D6", # lt purple
  "#FDBF6F", # lt orange
  "gray70", "khaki2",
  "maroon", "orchid1", "deeppink1", "blue1", "steelblue4",
  "darkturquoise", "green1", "yellow4", "yellow3",
  "darkorange4", "brown"
)
  
  c16 <- c("#8DD3C7",
           "#FFFFB3",
           "#BEBADA",
           "#FB8072", 
           "#80B1D3",
           "#FDB462",
           "#B3DE69",
           "#FCCDE5",
           "#D9D9D9",
           "#BC80BD", 
           "#CCEBC5", 
           "#FFED6F", 
           "#E41A1C",
           "#377EB8",
           "#4DAF4A", 
           "#0EAD9A")
  
  row_gfam_ha <- rowAnnotation(og=transporters_l2f$gfam,
                               col = list(family=c("ZIP" = "#8DD3C7",
                                                   "HMA" = "#FFFFB3",
                                                   "MTP" = "#BEBADA",
                                                   "VIT" = "#FB8072", 
                                                   "PCR" = "#80B1D3",
                                                   "NRAMP" = "#FDB462",
                                                   "HMP" = "#B3DE69",
                                                   "MATE" = "#FCCDE5",
                                                   "MRS2" = "#D9D9D9",
                                                   "CorA" = "#BC80BD", 
                                                   "ATX" = "#CCEBC5", 
                                                   "COPT" = "#FFED6F", 
                                                   "YSL" = "#E41A1C",
                                                   "CAX" = "#377EB8",
                                                   "ZIP" = "#4DAF4A", 
                                                   'ABC'="#0EAD9A")))
  
  # Load required packages
  require(data.table)
  require(RColorBrewer)
  require(circlize)
  
  #load cr data
  if (!exists("cr")) {
    print("The 'cr' data frame is not loaded in the environment.")
    load('cr.rda')
  }
  
  # Combine root and leaf log2 fold change data into separate data frames
  root_l2f <- do.call(rbind, lapply(lis, function(x) x['Rtl2f'][[1]]))
  setkey(root_l2f, geneid)
  leaf_l2f <- do.call(rbind, lapply(lis, function(x) x['Lfl2f'][[1]]))
  setkey(leaf_l2f, geneid)
  # Combine root and leaf padj data into separate data frames
  root_padj <- do.call(rbind, lapply(lis, function(x) x['Rtpadj'][[1]]))
  setkey(root_padj, geneid)
  leaf_padj <- do.call(rbind, lapply(lis, function(x) x['Lfpadj'][[1]]))
  setkey(leaf_padj, geneid)
  
  
  # Rearrange columns for root and leaf l2f
  root_l2f <- root_l2f[,c(3:5, 1:2)]
  leaf_l2f <- leaf_l2f[,c(3:5, 1:2)]
  
  # Initialize l2f and padj 
  l2f <- root_l2f[leaf_l2f]
  padj <- root_padj[leaf_padj]
  keep <- apply(padj[, c(3:5, 7:9)], 1, function(x) any(x < 0.05))
  l2f <- l2f[keep, ]
  padj <- padj[keep, ]
  l2f[, gfam := fcoalesce(as.character(gfam), as.character(i.gfam) )]
  padj[, 'gfam'] <- l2f$gfam
  
  
  colsplt <- factor(rep(c("Root", "Leaf"), each=3), ordered = T, levels = c('Root', 'Leaf'))
  idx <- c(1:3, 6:8)
  
  lis_order <- unlist(lapply(lis, function(x) x[[1]][, geneid]))
  
  l2f <- l2f[match(lis_order, geneid, nomatch = 0)]
  padj <- padj[match(lis_order, geneid, nomatch = 0)]
  
  
  steps=c("HMA", 'ZIP', 'VIT', 'PCR', 'NRAMP', "MetChap" , 'MATE', 'MRS2', 'CorA', "COPT",'YSL', 'CAX', 'ABC')
  #steps <- steps[order(steps)]
  num_repeats=c( 4, 3, 3, 8, 2, 8, 14, 3, 2, 2, 2, 4, 27)
  
  STEP <- character()
  
  for (i in seq_along(steps)) {
    #num_repeats <- length(grep(tolower(steps[[i]]), tolower(l2f$gfam)))
    print(num_repeats[i])
    STEP <- c(STEP, rep(steps[i], times = num_repeats[i]))
    
  }
  
  # Define colors using RColorBrewer
  #num_categories <- length(unique(STEP))
  
  # Convert STEP into a named vector with colors assigned to each unique cate
  
  row_gfam_ha <- rowAnnotation(family=STEP, 
                               col = list(Family=c("ZIP" = "#8DD3C7",   "HMA" = "#FFFFB3",   "MTP" = "#BEBADA",
                                                   "VIT" = "#FB8072", "PCR" = "#80B1D3", "NRAMP" = "#FDB462", 
                                                   "MetChap" = "#B3DE69", "MATE" = "#FCCDE5", "MRS2" = "#D9D9D9", 
                                                   "CorA" = "#BC80BD", "COPT" = "#FFED6F", 
                                                   "YSL" = "#E41A1C", "CAX" = "#377EB8", "ZIP" = "#4DAF4A")))#, 'ABC'="#0EAD9A")))
  
  row_gfam_ha <- rowAnnotation(family=STEP, 
                               col = list(family=c("ZIP" = "#8DD3C7",   "HMA" = "#FFFFB3",   "MTP" = "#BEBADA",
                                                   "VIT" = "#FB8072", "PCR" = "#80B1D3", "NRAMP" = "#FDB462", 
                                                   "MetChap" = "#B3DE69", "MATE" = "#FCCDE5", "MRS2" = "#D9D9D9", 
                                                   "CorA" = "#BC80BD", "COPT" = "#FFED6F", 
                                                   "YSL" = "#E41A1C", "CAX" = "#377EB8", "ZIP" = "#4DAF4A", 'ABC'="#0EAD9A")))
  

  #rownames
  rowns <- df_res$root_48h_0h[match(l2f[, geneid] , geneid), ncbi]
  rowns <- paste(rowns, tmh[match(l2f$geneid, geneid), abbreviations], sep = '-')
  
  #legend for *
  lgd_sig = Legend(pch = c("*", '**', '***'), type = "points", 
                   labels = c('< 0.05', '< 0.01', '< 0.001'), 
                   title = 'Significance')
  
  rowsplit=factor(as.character(l2f$gfam), levels=unique(l2f$gfam))
  
  # heatmap
  set.seed(1300)
  hm <- Heatmap(matrix = as.matrix(l2f[,..idx]), name = "Log2FoldChange",cluster_columns = F, cluster_rows = F, row_labels = rowns,
                col=colorRamp2(c(min(l2f[,..idx], na.rm = T), 0, max(l2f[,..idx], na.rm = T)), colors = c("red", "white", "green")),
                #column_labels = c('RtLA48h/0h', 'RtHA48Hh/0h'),
                row_names_gp = gpar(fontsize=10),
                left_annotation = row_gfam_ha,
                #heatmap_legend_param = ,
                border = T,
                column_gap = unit(0, "mm"),
                column_split = colsplt,
                row_split = rowsplit,
                cell_fun = function(j, i, x, y, width, height, fill) {
                  if(is.na(as.matrix(padj[, c(-1, -2, -6,-10)])[i, j]) ){
                    grid.text(sprintf("%s", ''), x, y, gp = gpar(fontsize = 12))
                  }else if(as.matrix(padj[, c(-1, -2, -6,-10)])[i, j] <0.001 ){
                    grid.text(sprintf("%s", '***'), x, y, gp = gpar(fontsize = 12))
                  }else if(as.matrix(padj[, c(-1, -2, -6,-10)])[i, j] <0.01) {
                    grid.text(sprintf("%s", '**'), x, y, gp = gpar(fontsize = 12))
                  }else if(as.matrix(padj[, c(-1, -2, -6,-10)])[i, j] <0.05) {
                    grid.text(sprintf("%s", '*'), x, y, gp = gpar(fontsize = 12 ))
                  }
                }
  )
  
  pdf(file = 'metal_transporters5.pdf', width = 5.5, height = 10)
  hm
  draw(lgd_sig, x = unit(.88, "npc"), y = unit(.4, "npc"), just='center')
  dev.off()
  
  pdf(file = 'metal_transporters6.pdf', width = 5.5, height = 13)
  hm
  draw(lgd_sig, x = unit(.88, "npc"), y = unit(.4, "npc"), just='center')
  dev.off()
  
  # Define row labels including the Arabidopsis RBBH
  # Create second Arabidopsis RBBH labels
  
  araNames <- cr[crv3.2 %in% unique(padj$geneid), .(crv3.2, araport2)]
  
  araNames[is.na(araport2), 'araport2'] <- ''
  
  araNames <- unique(araNames)
  
  
  if(length(rowns)<nrow(padj)){
    
    # Create a new rowns_transport vector
    rowns <- character(nrow(padj))  # Initialize with empty character vector
    
    # Identify the indices where summary_stats$gene is in rowns_transport$crv3.2
    match_indices <- rowns %in% cr$crv3.2
    
    # Assign the original labels where there's a match
    
    rowns[!match_indices] <- padj[!match_indices,geneid]
    rowns[match_indices] <- substring((cr[match(padj[match_indices,geneid], crv3.2),][!is.na(cr), ncbi]), first = 5, last = 18)
    
  
  #define cell_fun
  cellfun = function(j, i, x, y, width, height, fill) {
    if(is.na(as.matrix(padj[, c(-1, -2, -6)])[i, j]) ){
      grid.text(sprintf("%s", ''), x, y, gp = gpar(fontsize = 12))
    }else if(as.matrix(padj[, c(-1, -2, -6)])[i, j] <0.001 ){
      grid.text(sprintf("%s", '***'), x, y, gp = gpar(fontsize = 12))
    }else if(as.matrix(padj[, c(-1, -2, -6)])[i, j] <0.01) {
      grid.text(sprintf("%s", '**'), x, y, gp = gpar(fontsize = 12))
    }else if(as.matrix(padj[, c(-1, -2, -6)])[i, j] <0.05) {
      grid.text(sprintf("%s", '*'), x, y, gp = gpar(fontsize = 12 ))
    }
  }
  # Make Heatmap left annotation
  
  
  



padj[padj[, .I[apply(.SD, 1, function(x) any(x < 0.05))], .SDcols = c("BOTH.x",   "TSH660.x", "PA121.x",  "BOTH.y",   "TSH660.y", "PA121.y")]]

apply(padj[, ..idx], 1, function(x) any(x < 0.05))

keep <- !is.na(apply(padj[, ..idx], 1, function(x) any(x < 0.05)))


l2f <- l2f[keep, ]
padj <- padj[keep, ]


rowns <- df_res$root_48h_0h[ match(l2f[, geneid], geneid), ncbi]
          


hm <- Heatmap(matrix = as.matrix(l2f[,..idx]), name = "Log2FoldChange",cluster_columns = F, cluster_rows = F, row_labels = rowns,
              col=colorRamp2(c(min(l2f[,..idx], na.rm = T), 0, max(l2f[,..idx], na.rm = T)), colors = c("red", "white", "green")),
              #column_labels = c('RtLA48h/0h', 'RtHA48Hh/0h'),
              row_names_gp = gpar(fontsize=10),
              left_annotation = row_gfam_ha,
              #heatmap_legend_param = ,
              border = T,
              column_gap = unit(0, "mm"),
              column_split = colsplt,
              row_split = l2f[[gfam_col]],
              cell_fun = cellfun
)


hm <- Heatmap(matrix = as.matrix(l2f[,..idx]), name = "Log2FoldChange",cluster_columns = F, cluster_rows = F, row_labels = rowns_transport,
              col=colorRamp2(c(min(l2f[,..idx], na.rm = T), 0, max(l2f[,..idx], na.rm = T)), colors = c("red", "white", "green")),
              #column_labels = c('RtLA48h/0h', 'RtHA48Hh/0h'),
              row_names_gp = gpar(fontsize=10),
              #left_annotation = row_gfam_ha,
              #heatmap_legend_param = ,
              border = T,
              column_gap = unit(0, "mm"),
              column_split = colsplt,
              row_split = l2f[[gfam_col]],
              cell_fun = function(j, i, x, y, width, height, fill) {
                if(is.na(as.matrix(padj)[i, j]) ){
                  grid.text(sprintf("%s", ''), x, y, gp = gpar(fontsize = 12))
                }else if(as.matrix(padj)[i, j] <0.001 ){
                  grid.text(sprintf("%s", '***'), x, y, gp = gpar(fontsize = 12))
                }else if(as.matrix(padj)[i, j] <0.01) {
                  grid.text(sprintf("%s", '**'), x, y, gp = gpar(fontsize = 12))
                }else if(as.matrix(padj)[i, j] <0.05) {
                  grid.text(sprintf("%s", '*'), x, y, gp = gpar(fontsize = 12 ))
                }
              }
)
