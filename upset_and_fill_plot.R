# make figure 2, DEG results
library(dplyr)
library(ggplot2)
library(magrittr)
library(data.table)
library(patchwork)
library(ComplexUpset)
load('df_leaf_geno.rda')
load('df_root_geno.rda')
load('cr.rda')
load('go_terms_cacao.rda')
load('df_res.rda')

cr <- data.table(cr)

source('../Agroavia_PSU_collab/ORA/Intersect_complement_functions.R') #Source very important functions for later


annot <- read.table('TC_B97_consensusV2_rev2021_annotated.gff3', skip = 2, fill = T, stringsAsFactors = F)
annot <- annot[grep('gene', annot$V3), ]
annot$id <- sapply(strsplit(annot$V9, split = ';') , function(x) x[1])
annot$id <- sapply(strsplit(annot$id, split = '=') , function(x) x[2])
annot$annot <- apply(annot, 1, function(x) paste(x[9:length(x)], collapse = ' ')) #annotations
annot <- annot[!is.na(annot$id), ]

#Create the matrix of upregulated and downregulated genes


deg_tf_matrix <- data.frame(geneid=annot$id, #create data frame for the gene id matrix
                            leaf_pa121_24h_0h=NA,
                            leaf_tsh660_24h_0h=NA,
                            leaf_24h_0h=NA,
                            leaf_pa121_48h_0h=NA,
                            leaf_tsh660_48h_0h=NA,
                            leaf_48h_0h=NA,
                            leaf_cd_no_cd=NA,
                            leaf_pa121_tsh660_0h_0h=NA,
                            leaf_pa121_tsh660_24h_24h=NA,
                            leaf_pa121_tsh660_48h_48h=NA,
                            
                            root_pa121_48h_0h=NA,
                            root_tsh660_48h_0h=NA,
                            root_48h_0h=NA,
                            
                            root_pa121_tsh660_0h_0h=NA,
                            root_pa121_tsh660_48h_48h=NA)

#root_geno_comparison_names <- c('root_pa121_48h_0h', 'root_tsh660_48h_0h', 'root_pa121_tsh660_0h_0h', 'root_pa121_tsh660_48h_48h', 'root_48h_0h')


deg_tf_matrix <- data.table::data.table(deg_tf_matrix)

upreg_tf <- deg_tf_matrix #one for the upregulated genes 

downreg_tf <- deg_tf_matrix # one for the down regulated genes

df_res <- df_leaf_geno[c('leaf_pa121_24h_0h', 'leaf_tsh660_24h_0h','leaf_pa121_48h_0h','leaf_tsh660_48h_0h','leaf_48h_0h','leaf_24h_0h', 'leaf_pa121_tsh660_0h_0h', 'leaf_pa121_tsh660_24h_24h','leaf_pa121_tsh660_48h_48h')]
df_res <- append(df_res, df_root_geno[c('root_pa121_48h_0h','root_tsh660_48h_0h','root_pa121_tsh660_0h_0h','root_pa121_tsh660_48h_48h','root_48h_0h')]) #list of results data frames for the comparisons of interest

df_res <- lapply(df_res, function(x) return(data.table(x)))

deg_id_list <- lapply(df_res, function(x) x[x$padj<0.05, 'geneid']) #the ids of all deg in leaves
upreg_deg_id <- lapply(df_res, function(x) x[x$padj<0.05 & x$log2FoldChange>0, 'geneid']) #good
downreg_deg_id <- lapply(df_res, function(x) x[x$padj<0.05 & x$log2FoldChange<0, 'geneid']) #good

deg_id_list <- lapply(deg_id_list, function(x) return(x[!is.na(geneid)]))
upreg_deg_id <- lapply(upreg_deg_id, function(x) x[!is.na(geneid)])
downreg_deg_id <- lapply(downreg_deg_id, function(x)x[!is.na(geneid)])

sum(deg_tf_matrix$geneid %in% deg_id_list[['root_48h_0h']][,geneid])


invisible(lapply(names(deg_tf_matrix[, 2:ncol(deg_tf_matrix)]), function(x) {
  deg_tf_matrix[, x] <<- deg_tf_matrix$geneid %in% deg_id_list[[x]]$geneid
}))

invisible(lapply(names(upreg_tf[, 2:ncol(upreg_tf)]), function(x) {
  upreg_tf[, x] <<- upreg_tf$geneid %in% upreg_deg_id[[x]]$geneid
}))

invisible(lapply(names(downreg_tf[, 2:ncol(downreg_tf)]), function(x) {
  downreg_tf[, x] <<- downreg_tf$geneid %in% downreg_deg_id[[x]]$geneid
}))

deg_tf_matrix <- data.table(deg_tf_matrix)
upreg_tf <- data.table(upreg_tf)
downreg_tf <- data.table(downreg_tf)

#upset



in2way <- function(name1, name2, new_col='intersection'){
  fill_deg[(upreg_tf[, name1] == 1 & upreg_tf[, name2]==1), paste(new_col)] <<- 'up'
  fill_deg[(downreg_tf[, name1] == 1 & downreg_tf[, name2]==1), paste(new_col)] <<- 'down'
  fill_deg[(upreg_tf[, name1] == 1 & downreg_tf[, name2]==1), paste(new_col)] <<- 'opposite'
  fill_deg[(downreg_tf[, name1] == 1 & upreg_tf[, name2]==1), paste(new_col)] <<- 'opposite'
}


in3way <- function(name1, name2, name3, new_col='intersection'){
  fill_deg[upreg_tf[, name1]==1 & upreg_tf[,name2]==1 & upreg_tf[, name3]==1, paste(new_col)] <<- 'up'
  fill_deg[downreg_tf[, name1]==1 & downreg_tf[,name2]==1 & downreg_tf[, name3]==1, paste(new_col)] <<- 'down'
  fill_deg[upreg_tf[, name1]==1 & downreg_tf[,name2]==1 & downreg_tf[, name3]==1, paste(new_col)] <<- 'opposite'
  fill_deg[downreg_tf[, name1]==1 & upreg_tf[,name2]==1 & downreg_tf[, name3]==1, paste(new_col)] <<- 'opposite'
  fill_deg[downreg_tf[, name1]==1 & downreg_tf[,name2]==1 & upreg_tf[, name3]==1, paste(new_col)] <<- 'opposite'
  fill_deg[upreg_tf[, name1]==1 & downreg_tf[,name2]==1 & upreg_tf[, name3]==1, paste(new_col)] <<- 'opposite'
  fill_deg[upreg_tf[, name1]==1 & upreg_tf[,name2]==1 & downreg_tf[, name3]==1, paste(new_col)] <<- 'opposite'
}



fill_deg <- cbind(deg_tf_matrix[, .(geneid)], data.table(root_48h_0h = NA,
                                                         root_tsh660_48h_0h = NA,
                                                         root_pa121_48h_0h = NA,
                                                         root_pa121_tsh660_0h_0h = NA,
                                                         root_pa121_tsh660_48h_48h = NA,
                                                         leaf_pa121_24h_0h=NA,
                                                         leaf_tsh660_24h_0h=NA,
                                                         leaf_pa121_48h_0h=NA,
                                                         leaf_tsh660_48h_0h=NA,
                                                         leaf_48h_0h=NA,
                                                         leaf_24h_0h=NA,
                                                         leaf_pa121_tsh660_0h_0h = NA,
                                                         leaf_pa121_tsh660_24h_24h = NA,
                                                         leaf_pa121_tsh660_48h_48h = NA))


fill_deg <- as.data.frame(fill_deg)
upreg_tf <- as.data.frame(upreg_tf)
downreg_tf <- as.data.frame(downreg_tf)

lapply(colnames(fill_deg)[-1], function(x){
  print(x)
  fill_deg[upreg_tf[,x] ==  1, paste(x)] <<- 'up'
  fill_deg[downreg_tf[,x] == 1, paste(x)] <<- 'down'
})



in2way(name1 = 'root_tsh660_48h_0h', name2 = 'root_pa121_48h_0h', new_col = 'intersection_1')
in3way(name1 = 'root_48h_0h', name2 = 'root_pa121_48h_0h', name3 = 'root_tsh660_48h_0h', new_col = 'intersection_2')
in2way(name1 = 'leaf_tsh660_24h_0h', name2 = 'leaf_pa121_24h_0h', new_col = 'intersection_3')
in3way(name1 = 'leaf_24h_0h', name2 = 'leaf_pa121_24h_0h', name3 = 'leaf_tsh660_24h_0h', new_col = 'intersection_4')
in2way(name1 = 'leaf_tsh660_48h_0h', name2 = 'leaf_pa121_48h_0h', new_col = 'intersection_5')
in3way(name1 = 'leaf_48h_0h', name2 = 'leaf_pa121_48h_0h', name3 = 'leaf_tsh660_48h_0h', new_col = 'intersection_6')
in2way(name1 = 'leaf_tsh660_48h_0h', name2 = 'leaf_tsh660_24h_0h', new_col = 'intersection_7')
in2way(name1 = 'leaf_pa121_48h_0h', name2 = 'leaf_pa121_24h_0h', new_col = 'intersection_8')
in2way(name1 = 'root_pa121_tsh660_0h_0h', name2 = 'root_pa121_tsh660_48h_48h', new_col = 'intersection_9')




fill_deg <- reshape2::melt(data = fill_deg[,], id.vars = 'geneid')
fill_deg <- na.omit(fill_deg)



fill_deg$variable <- factor( fill_deg$variable, levels = unique(fill_deg$variable)[c(1:3, 11, 7, 6, 10, 9, 8, 4:5, 12:13, 14:23)], labels = c('Root Both genotypes 48h ACT (1)',
                                                                                                                                              'Root TSH660 48h vs 0h ACT (2)',
                                                                                                                                              'Root PA121 48h vs 0h ACT (3)',
                                                                                                                                              'Leaf Both genotypes 24h ACT (4)',
                                                                                                                                              'Leaf TSH660 24h vs 0h ACT (5)',
                                                                                                                                              'Leaf PA121 24h vs 0h ACT (6)',
                                                                                                                                              'Leaf Both genotypes 48h ACT (7)',
                                                                                                                                              'Leaf TSH660 48h vs 0h ACT (8)',
                                                                                                                                              'Leaf PA121 48h vs 0h ACT (9)',
                                                                                                                                              'Root PA121 0h vs TSH660 0h ACT (10)',
                                                                                                                                              'Root PA121 48h vs TSH660 48h ACT (11)',
                                                                                                                                              'Leaf PA121 0h vs TSH660 0h ACT (12)',
                                                                                                                                              'Leaf PA121 24h vs TSH660 24h ACT (13)',
                                                                                                                                              'Leaf PA121 48h vs TSH660 48h ACT (14)',
                                                                                                                                              '(2) in (3)', #root 48 two gneotypes INT1
                                                                                                                                              '(1) in (2) in (3)', #including the both genotype comparison INT2
                                                                                                                                              '(5) in (6)', #leaf 24h INT3
                                                                                                                                              '(4) in (5) in (6)', #leaf 24h plus both INT4
                                                                                                                                              '(8) in (9)', #leaf 48h INT5
                                                                                                                                              '(7) in (8) in (9)', #leaf 48h plus bth INT6
                                                                                                                                              '(5) in (8)', #leaf 24 and 48 INT7
                                                                                                                                              '(6) in (9)', #leaf 24 and 48 INT8
                                                                                                                                              '(13) in (14)')) #leaf 24 diff and 48h diff INT9



fill_deg$variable2 <- factor( fill_deg$variable, levels = levels(fill_deg$variable), labels = c('Rt48h/Rt0h',
                                                                                                'RtLA48h/RtLA0h',
                                                                                                'RtHA48h/RtHA0h',
                                                                                                'Lf24h/Lf0h',
                                                                                                'LfLA24h/LfLA0h',
                                                                                                'LfHA24h/LfHA0h',
                                                                                                'Lf48h/Lf0h',
                                                                                                'LfLA48h/LfLA0h',
                                                                                                'LfHA48h/LfHA0h',
                                                                                                'RtHA0h/RtLA0h',
                                                                                                'RtHA48h/RtLA48h',
                                                                                                'LfHA0h/LfLA0h',
                                                                                                'LfHA24h/LfLA24h',
                                                                                                'LfHA48h/LfLA48h',
                                                                                                '(2) in (3)', #root 48 two gneotypes INT1
                                                                                                '(1) in (2) in (3)', #including the both genotype comparison INT2
                                                                                                '(5) in (6)', #leaf 24h INT3
                                                                                                '(4) in (5) in (6)', #leaf 24h plus both INT4
                                                                                                '(8) in (9)', #leaf 48h INT5
                                                                                                '(7) in (8) in (9)', #leaf 48h plus bth INT6
                                                                                                '(5) in (8)', #leaf 24 and 48 INT7
                                                                                                '(6) in (9)', #leaf 24 and 48 INT8
                                                                                                '(13) in (14)')) #leaf 24 diff and 48h diff INT9



intersections_upset <- list(c("leaf_pa121_24h_0h"), 
  c("leaf_tsh660_24h_0h"),
  c("leaf_pa121_48h_0h"),
  c("leaf_tsh660_48h_0h"),
  c("leaf_24h_0h"),
  c("leaf_48h_0h"),
  c("root_pa121_48h_0h"),
  c("root_tsh660_48h_0h"),
  c("root_48h_0h"),
  c('root_pa121_tsh660_0h_0h'),
  c('root_pa121_tsh660_48h_48h'),
  c('leaf_pa121_tsh660_0h_0h'),
  c('leaf_pa121_tsh660_24h_24h'),
  c('leaf_pa121_tsh660_48h_48h'), 
  c('root_tsh660_48h_0h', 'root_pa121_48h_0h'),
  c('root_48h_0h', 'root_tsh660_48h_0h', 'root_pa121_48h_0h'),
  c('leaf_tsh660_24h_0h', 'leaf_pa121_24h_0h' ),
  c('leaf_24h_0h', 'leaf_tsh660_24h_0h', 'leaf_pa121_24h_0h'),
  c('leaf_tsh660_48h_0h', 'leaf_pa121_48h_0h'),
  c('leaf_48h_0h', 'leaf_tsh660_48h_0h', 'leaf_pa121_48h_0h'),
  c('leaf_tsh660_24h_0h', 'leaf_tsh660_48h_0h'),
  c('leaf_pa121_24h_0h', 'leaf_pa121_48h_0h'),
  c('leaf_pa121_tsh660_24h_24h', 'leaf_pa121_tsh660_48h_48h'))

intersections_upset2 <- list(c('Root Both genotypes 48h ACT (1)'),
                            c('Root TSH660 48h vs 0h ACT (2)'),
                            c('Root PA121 48h vs 0h ACT (3)'),
                            c('Leaf Both genotypes 24h ACT (4)'),
                            c('Leaf TSH660 24h vs 0h ACT (5)'),
                            c('Leaf PA121 24h vs 0h ACT (6)'),
                            c('Leaf Both genotypes 48h ACT (7)'),
                            c('Leaf TSH660 48h vs 0h ACT (8)'),
                            c('Leaf PA121 48h vs 0h ACT (9)'),
                            c('Root PA121 0h vs TSH660 0h ACT (10)'),
                            c('Root PA121 48h vs TSH660 48h ACT (11)'),
                            c('Leaf PA121 0h vs TSH660 0h ACT (12)'),
                            c('Leaf PA121 24h vs TSH660 24h ACT (13)'),
                            c('Leaf PA121 48h vs TSH660 48h ACT (14)'),
                            c('Root TSH660 48h vs 0h ACT (2)', 'Root PA121 48h vs 0h ACT (3)'),
                            c('Root Both genotypes 48h ACT (1)', 'Root TSH660 48h vs 0h ACT (2)', 'Root PA121 48h vs 0h ACT (3)'),
                            c('Leaf TSH660 24h vs 0h ACT (5)', 'Leaf PA121 24h vs 0h ACT (6)' ),
                            c('Leaf Both genotypes 24h ACT (4)', 'Leaf TSH660 24h vs 0h ACT (5)', 'Leaf PA121 24h vs 0h ACT (6)'),
                            c('Leaf TSH660 48h vs 0h ACT (8)', 'Leaf PA121 48h vs 0h ACT (9)'),
                            c('Leaf Both genotypes 48h ACT (7)', 'Leaf TSH660 48h vs 0h ACT (8)', 'Leaf PA121 48h vs 0h ACT (9)'),
                            c('Leaf TSH660 24h vs 0h ACT (5)', 'Leaf TSH660 48h vs 0h ACT (8)'),
                            c('Leaf PA121 24h vs 0h ACT (6)', 'Leaf PA121 48h vs 0h ACT (9)'),
                            c('Leaf PA121 24h vs TSH660 24h ACT (13)', 'Leaf PA121 48h vs TSH660 48h ACT (14)')
)


intersections_upset3 <- list(c('Rt48h/Rt0h'),
                             c('RtLA48h/RtLA0h'),
                             c('RtHA48h/RtHA0h'),
                             c('Lf24h/Lf0h'),
                             c('LfLA24h/LfLA0h'),
                             c('LfHA24h/LfHA0h'),
                             c('Lf48h/Lf0h'),
                             c('LfLA48h/LfLA0h'),
                             c('LfHA48h/LfHA0h'),
                             c('RtHA0h/RtLA0h'),
                             c('RtHA48h/RtLA48h'),
                             c('LfHA0h/LfLA0h'),
                             c('LfHA24h/LfLA24h'),
                             c('LfHA48h/LfLA48h'),
                             c('RtLA48h/RtLA0h', 'RtHA48h/RtHA0h'),
                             c('Rt48h/Rt0h', 'RtLA48h/RtLA0h', 'RtHA48h/RtHA0h'),
                             c('LfLA24h/LfLA0h', 'LfHA24h/LfHA0h'),
                             c('Lf24h/Lf0h', 'LfLA24h/LfLA0h', 'LfHA24h/LfHA0h'),
                            c('LfLA48h/LfLA0h', 'LfHA48h/LfHA0h'),
                            c('Lf24h/Lf0h', 'LfLA48h/LfLA0h', 'LfHA48h/LfHA0h'),
                            c('LfLA24h/LfLA0h', 'LfLA48h/LfLA0h'),
                            c('LfHA24h/LfHA0h', 'LfHA48h/LfHA0h'),
                            c('LfHA24h/LfLA24h', 'LfHA48h/LfLA48h')
)



labs_upset <- c('Root Both genotypes 48h ACT (1)',
                'Root TSH660 48h vs 0h ACT (2)',
                'Root PA121 48h vs 0h ACT (3)',
                'Leaf Both genotypes 24h ACT (4)',
                'Leaf TSH660 24h vs 0h ACT (5)',
                'Leaf PA121 24h vs 0h ACT (6)',
                'Leaf Both genotypes 48h ACT (7)',
                'Leaf TSH660 48h vs 0h ACT (8)',
                'Leaf PA121 48h vs 0h ACT (9)',
                'Root PA121 0h vs TSH660 0h ACT (10)',
                'Root PA121 48h vs TSH660 48h ACT (11)',
                'Leaf PA121 0h vs TSH660 0h ACT (12)',
                'Leaf PA121 24h vs TSH660 24h ACT (13)',
                'Leaf PA121 48h vs TSH660 48h ACT (14)',
                '(2) in (3)', #root 48 two gneotypes INT1
                '(1) in (2) in (3)', #including the both genotype comparison INT2
                '(5) in (6)', #leaf 24h INT3
                '(4) in (5) in (6)', #leaf 24h plus both INT4
                '(8) in (9)', #leaf 48h INT5
                '(7) in (8) in (9)', #leaf 48h plus bth INT6
                '(5) in (8)', #leaf 24 and 48 INT7
                '(6) in (9)', #leaf 24 and 48 INT8
                '(13) in (14)') #leaf 24 diff and 48h diff INT9

fill_plot <- fill_deg %>% 
  count(variable, value) %>%
  group_by(variable) %>%
  mutate(pct=prop.table(n) * 100) %>%
  ggplot() + aes(variable, pct, fill=value)+
  geom_bar(stat='identity')+
  theme_minimal()+
  theme(axis.text.x = element_text(angle = 45, size = 5, vjust = 0.865, hjust = 0.85), axis.text.y = element_text(size = 7), axis.title.y = element_text(size = 9, hjust = .5, vjust=.5), axis.title.x = element_text(size = 7), legend.text = element_text(size=5), legend.title = element_text(size=7))+
  ylab(label='Proportion of upreglated,\ndownregulated and opposite\nregulated DEG in set')+
  xlab(label = 'Gene expression comparison and intersection (in) set')+
  labs(fill = 'Regulation\npolarity')+
  geom_text(aes(label=paste0(sprintf("%1.1f", pct),"%")), size=2.5,
            position=position_stack(vjust=0.5))

ggsave(filename = 'fill_plot_PAG.pdf', plot = fill_plot, device = 'pdf', width = 9, height = 4.5, dpi = 300)


upset_1 <- ComplexUpset::upset(data = deg_tf_matrix,
                    wrap = T,
                    width_ratio = 0.1, 
                    intersect = 
                      c("leaf_pa121_24h_0h", 
                        "leaf_tsh660_24h_0h",
                        "leaf_pa121_48h_0h",
                        "leaf_tsh660_48h_0h",
                        "leaf_24h_0h",
                        "leaf_48h_0h",
                        "root_pa121_48h_0h",
                        "root_tsh660_48h_0h",
                        "root_48h_0h",
                        'root_pa121_tsh660_0h_0h',
                        'root_pa121_tsh660_48h_48h',
                        'leaf_pa121_tsh660_0h_0h',
                        'leaf_pa121_tsh660_24h_24h',
                        'leaf_pa121_tsh660_48h_48h'), 
                    mode = 'inclusive_intersection',
                    name = 'Gene Contrast Intersection', 
                    keep_empty_groups=T,
                    sort_sets=F,
                    set_sizes=F,
                    sort_intersections=F, 
                    intersections = intersections_upset,
                    base_annotations=list('Intersection size'=intersection_size(mode = 'inclusive_intersection', text = list(size=2))),
                    themes = upset_modify_themes(list('intersections_matrix' = theme(axis.text.y=element_text(size=5), axis.title.x = element_text(size=7)),
                                                      'Intersection size' = list(theme(axis.title=element_text(size=7)))))
                    #labeller=ggplot2::as_labeller(c(
                    #  "root_48h_0h" = "Both Genotype 48h vs 0h ACT",
                    #  "root_tsh660_48h_0h" = "TSH660 48h vs 0h ACT",
                    #  "root_pa121_48h_0h" = 'PA121 48h vs 0h ACT')),
                    
                    ) 
)


upset_2 <- ComplexUpset::upset(data = deg_tf_matrix,
                               wrap = T,
                               width_ratio = 0.1, 
                               intersect = 
                                 c("leaf_pa121_24h_0h", 
                                   "leaf_tsh660_24h_0h",
                                   "leaf_pa121_48h_0h",
                                   "leaf_tsh660_48h_0h",
                                   "leaf_24h_0h",
                                   "leaf_48h_0h",
                                   "root_pa121_48h_0h",
                                   "root_tsh660_48h_0h",
                                   "root_48h_0h",
                                   'root_pa121_tsh660_0h_0h',
                                   'root_pa121_tsh660_48h_48h',
                                   'leaf_pa121_tsh660_0h_0h',
                                   'leaf_pa121_tsh660_24h_24h',
                                   'leaf_pa121_tsh660_48h_48h'), 
                               mode = 'inclusive_intersection',
                               name = 'Gene Contrast Intersection', 
                               keep_empty_groups=T,
                               sort_sets=F,
                               set_sizes=F,
                               sort_intersections=F, 
                               intersections = intersections_upset,
                               base_annotations=list('Intersection size'=intersection_size(mode = 'inclusive_intersection', text = list(size=2))),
                               themes = upset_modify_themes(list('intersections_matrix' = theme(axis.text.y=element_text(size=5), axis.title.x = element_text(size=7)),
                                                                 'Intersection size' = list(theme(axis.title=element_text(size=7))))),
                               #annotations = list('Polarity'=(fill_plot)),
                               labeller=ggplot2::as_labeller(c(
                                 "root_48h_0h" = "Root Both genotypes 48h ACT (1)",
                                 "root_tsh660_48h_0h" = "Root TSH660 48h vs 0h ACT (2)",
                                 "root_pa121_48h_0h" = "Root PA121 48h vs 0h ACT (3)",
                                 "leaf_24h_0h" = "Leaf Both genotypes 24h ACT (4)",
                                 "leaf_tsh660_24h_0h" = "Leaf TSH660 24h vs 0h ACT (5)",
                                 "leaf_pa121_24h_0h" = "Leaf PA121 24h vs 0h ACT (6)",
                                 "leaf_48h_0h" = "Leaf Both genotypes 48h ACT (7)",
                                 "leaf_tsh660_48h_0h" = "Leaf TSH660 48h vs 0h ACT (8)",
                                 "leaf_pa121_48h_0h" = "Leaf PA121 48h vs 0h ACT (9)",
                                 "root_pa121_tsh660_0h_0h" = "Root PA121 0h vs TSH660 0h ACT (10)",
                                 "root_pa121_tsh660_48h_48h" = "Root PA121 48h vs TSH660 48h ACT (11)",
                                 "leaf_pa121_tsh660_0h_0h" = "Leaf PA121 0h vs TSH660 0h ACT (12)",
                                 "leaf_pa121_tsh660_24h_24h" = "Leaf PA121 24h vs TSH660 24h ACT (13)",
                                 "leaf_pa121_tsh660_48h_48h" = "Leaf PA121 48h vs TSH660 48h ACT (14)")),
                               
                               
) 





upset_3 <- ComplexUpset::upset(data = deg_tf_matrix,
                               wrap = T,
                               width_ratio = 0.1, 
                               intersect = 
                                 c("leaf_pa121_24h_0h", 
                                   "leaf_tsh660_24h_0h",
                                   "leaf_pa121_48h_0h",
                                   "leaf_tsh660_48h_0h",
                                   "leaf_24h_0h",
                                   "leaf_48h_0h",
                                   "root_pa121_48h_0h",
                                   "root_tsh660_48h_0h",
                                   "root_48h_0h",
                                   'root_pa121_tsh660_0h_0h',
                                   'root_pa121_tsh660_48h_48h',
                                   'leaf_pa121_tsh660_0h_0h',
                                   'leaf_pa121_tsh660_24h_24h',
                                   'leaf_pa121_tsh660_48h_48h'), 
                               mode = 'inclusive_intersection',
                               name = 'Gene Contrast Intersection', 
                               keep_empty_groups=T,
                               sort_sets=F,
                               set_sizes=F,
                               sort_intersections=F, 
                               intersections = intersections_upset,
                               base_annotations=list('Intersection size'=intersection_size(mode = 'inclusive_intersection', text = list(size=2))),
                               themes = upset_modify_themes(list('intersections_matrix' = theme(axis.text.y=element_text(size=5), axis.title.x = element_text(size=7)),
                                                                 'Intersection size' = list(theme(axis.title=element_text(size=7))))),
                               #annotations = list('Polarity'=(fill_plot)),
                               labeller=ggplot2::as_labeller(c(
                                 "root_48h_0h" = "Rt48h/Rt0h",
                                 "root_tsh660_48h_0h" = "RtLA48h/RtLA0h",
                                 "root_pa121_48h_0h" = "RtHA48h/RtHA0h",
                                 "leaf_24h_0h" = "Lf24h/Lf0h",
                                 "leaf_tsh660_24h_0h" = "LfLA24h/LfLA0h",
                                 "leaf_pa121_24h_0h" = "LfHA24h/LfHA0h",
                                 "leaf_48h_0h" = "Lf48h/Lf0h",
                                 "leaf_tsh660_48h_0h" = "LfLA48h/LfLA0h",
                                 "leaf_pa121_48h_0h" = "LfHA48h/LfHA0h",
                                 "root_pa121_tsh660_0h_0h" = "RtHA0h/RtLA0h",
                                 "root_pa121_tsh660_48h_48h" = "RtHA48h/RtLA48h",
                                 "leaf_pa121_tsh660_0h_0h" = "LfHA0h/LfLA0h",
                                 "leaf_pa121_tsh660_24h_24h" = "LfHA24h/LfLA24h",
                                 "leaf_pa121_tsh660_48h_48h" = "LfHA48h/LfLA48h"),
                               
                               
) 
)



bar_inter <- c(4185, 2396, 2297, 1162, 57, 604, 1999, 180, 1538, 544, 530, 491, 884, 1183, 1378, 1377, 32, 32, 114, 90, 30, 414, 530)
bar_inter <- data.frame(x=levels(fill_deg$variable), y= bar_inter)
bar_inter$x <- factor(bar_inter$x, levels = bar_inter$x, ordered = T)


(ggplot(data = bar_inter, aes(x=x, y=y))+geom_bar(stat = 'identity'))/fill_plot

#(ggplot(data = bar_inter, aes(x=x, y=log10(y)))+geom_bar(stat = 'identity')) #log transformed 

dot_mat <- data.frame(y=factor(levels(fill_deg$variable)[1:14], levels=levels(fill_deg$variable)[14:1], ordered = T), x=factor(rep(levels(fill_deg$variable), each = 14), levels = levels(fill_deg$variable), ordered = T))
dot_mat <- data.table(dot_mat)

dot_mat_plot <- ggplot(data=dot_mat, aes(x=x, y=y))+geom_point(colour='grey')+
  theme_minimal()+
  theme(axis.text.x = element_blank(), axis.title.y = element_blank(), axis.title.x = element_blank(), panel.grid.major.x = element_blank(), panel.grid.major.y = element_blank())+
  geom_point(data=dot_mat[dot_mat$x==intersections_upset2[[1]] & dot_mat$y==intersections_upset2[[1]],], aes(x=x, y=y))+
  geom_point(data=dot_mat[dot_mat$x==intersections_upset2[[2]] & dot_mat$y==intersections_upset2[[2]],], aes(x=x, y=y))+
  geom_point(data=dot_mat[dot_mat$x==intersections_upset2[[3]] & dot_mat$y==intersections_upset2[[3]],], aes(x=x, y=y))+
  geom_point(data=dot_mat[dot_mat$x==intersections_upset2[[4]] & dot_mat$y==intersections_upset2[[4]],], aes(x=x, y=y))+
  geom_point(data=dot_mat[dot_mat$x==intersections_upset2[[5]] & dot_mat$y==intersections_upset2[[5]],], aes(x=x, y=y))+
  geom_point(data=dot_mat[dot_mat$x==intersections_upset2[[6]] & dot_mat$y==intersections_upset2[[6]],], aes(x=x, y=y))+
  geom_point(data=dot_mat[dot_mat$x==intersections_upset2[[7]] & dot_mat$y==intersections_upset2[[7]],], aes(x=x, y=y))+
  geom_point(data=dot_mat[dot_mat$x==intersections_upset2[[8]] & dot_mat$y==intersections_upset2[[8]],], aes(x=x, y=y))+
  geom_point(data=dot_mat[dot_mat$x==intersections_upset2[[9]] & dot_mat$y==intersections_upset2[[9]],], aes(x=x, y=y))+
  geom_point(data=dot_mat[dot_mat$x==intersections_upset2[[10]] & dot_mat$y==intersections_upset2[[10]],], aes(x=x, y=y))+
  geom_point(data=dot_mat[dot_mat$x==intersections_upset2[[11]] & dot_mat$y==intersections_upset2[[11]],], aes(x=x, y=y))+
  geom_point(data=dot_mat[dot_mat$x==intersections_upset2[[12]] & dot_mat$y==intersections_upset2[[12]],], aes(x=x, y=y))+
  geom_point(data=dot_mat[dot_mat$x==intersections_upset2[[13]] & dot_mat$y==intersections_upset2[[13]],], aes(x=x, y=y))+
  geom_point(data=dot_mat[dot_mat$x==intersections_upset2[[14]] & dot_mat$y==intersections_upset2[[14]],], aes(x=x, y=y))+
  
  geom_point(data=dot_mat[x %in% "(2) in (3)" & y %in% intersections_upset2[[15]]], aes(x=x, y=y))+
  geom_point(data=dot_mat[x %in% "(1) in (2) in (3)" & y %in% intersections_upset2[[16]]], aes(x=x, y=y))+
  geom_point(data=dot_mat[x %in% "(5) in (6)" & y %in% intersections_upset2[[17]]], aes(x=x, y=y))+
  geom_point(data=dot_mat[x %in% "(4) in (5) in (6)" & y %in% intersections_upset2[[18]]], aes(x=x, y=y))+
  geom_point(data=dot_mat[x %in% "(8) in (9)" & y %in% intersections_upset2[[19]]], aes(x=x, y=y))+
  geom_point(data=dot_mat[x %in% "(7) in (8) in (9)" & y %in% intersections_upset2[[20]]], aes(x=x, y=y))+
  geom_point(data=dot_mat[x %in% "(5) in (8)" & y %in% intersections_upset2[[21]]], aes(x=x, y=y))+
  geom_point(data=dot_mat[x %in% "(6) in (9)" & y %in% intersections_upset2[[22]]], aes(x=x, y=y))+
  geom_point(data=dot_mat[x %in% "(13) in (14)" & y %in% intersections_upset2[[23]]], aes(x=x, y=y))+
  
  geom_line(data=dot_mat[x %in% "(2) in (3)" & y %in% intersections_upset2[[15]]], aes(group=1))+
  geom_line(data=dot_mat[x %in% "(1) in (2) in (3)" & y %in% intersections_upset2[[16]]], aes(group=1))+
  geom_line(data=dot_mat[x %in% "(5) in (6)" & y %in% intersections_upset2[[17]]], aes(group=1))+
  geom_line(data=dot_mat[x %in% "(4) in (5) in (6)" & y %in% intersections_upset2[[18]]], aes(group=1))+
  geom_line(data=dot_mat[x %in% "(8) in (9)" & y %in% intersections_upset2[[19]]], aes(group=1))+
  geom_line(data=dot_mat[x %in% "(7) in (8) in (9)" & y %in% intersections_upset2[[20]]], aes(group=1))+
  geom_line(data=dot_mat[x %in% "(5) in (8)" & y %in% intersections_upset2[[21]]], aes(group=1))+
  geom_line(data=dot_mat[x %in% "(6) in (9)" & y %in% intersections_upset2[[22]]], aes(group=1))+
  geom_line(data=dot_mat[x %in% "(13) in (14)" & y %in% intersections_upset2[[23]]], aes(group=1))

dot_mat2 <- data.frame(y=factor(levels(fill_deg$variable2)[1:14], levels=levels(fill_deg$variable2)[14:1], ordered = T), x=factor(rep(levels(fill_deg$variable2), each = 14), levels = levels(fill_deg$variable2), ordered = T))
dot_mat2 <- data.table(dot_mat2)

dot_mat_plot2 <- ggplot(data=dot_mat2, aes(x=x, y=y))+geom_point(colour='grey')+
  theme_minimal()+
  theme(axis.text.x = element_blank(), axis.text.y = element_text(size=12), axis.title.y = element_blank(), axis.title.x = element_blank(), panel.grid.major.x = element_blank(), panel.grid.major.y = element_blank())+
  geom_point(data=dot_mat2[dot_mat2$x==intersections_upset3[[1]] & dot_mat2$y==intersections_upset3[[1]],], aes(x=x, y=y))+
  geom_point(data=dot_mat2[dot_mat2$x==intersections_upset3[[2]] & dot_mat2$y==intersections_upset3[[2]],], aes(x=x, y=y))+
  geom_point(data=dot_mat2[dot_mat2$x==intersections_upset3[[3]] & dot_mat2$y==intersections_upset3[[3]],], aes(x=x, y=y))+
  geom_point(data=dot_mat2[dot_mat2$x==intersections_upset3[[4]] & dot_mat2$y==intersections_upset3[[4]],], aes(x=x, y=y))+
  geom_point(data=dot_mat2[dot_mat2$x==intersections_upset3[[5]] & dot_mat2$y==intersections_upset3[[5]],], aes(x=x, y=y))+
  geom_point(data=dot_mat2[dot_mat2$x==intersections_upset3[[6]] & dot_mat2$y==intersections_upset3[[6]],], aes(x=x, y=y))+
  geom_point(data=dot_mat2[dot_mat2$x==intersections_upset3[[7]] & dot_mat2$y==intersections_upset3[[7]],], aes(x=x, y=y))+
  geom_point(data=dot_mat2[dot_mat2$x==intersections_upset3[[8]] & dot_mat2$y==intersections_upset3[[8]],], aes(x=x, y=y))+
  geom_point(data=dot_mat2[dot_mat2$x==intersections_upset3[[9]] & dot_mat2$y==intersections_upset3[[9]],], aes(x=x, y=y))+
  geom_point(data=dot_mat2[dot_mat2$x==intersections_upset3[[10]] & dot_mat2$y==intersections_upset3[[10]],], aes(x=x, y=y))+
  geom_point(data=dot_mat2[dot_mat2$x==intersections_upset3[[11]] & dot_mat2$y==intersections_upset3[[11]],], aes(x=x, y=y))+
  geom_point(data=dot_mat2[dot_mat2$x==intersections_upset3[[12]] & dot_mat2$y==intersections_upset3[[12]],], aes(x=x, y=y))+
  geom_point(data=dot_mat2[dot_mat2$x==intersections_upset3[[13]] & dot_mat2$y==intersections_upset3[[13]],], aes(x=x, y=y))+
  geom_point(data=dot_mat2[dot_mat2$x==intersections_upset3[[14]] & dot_mat2$y==intersections_upset3[[14]],], aes(x=x, y=y))+
  
  geom_point(data=dot_mat2[x %in% "(2) in (3)" & y %in% intersections_upset3[[15]]], aes(x=x, y=y))+
  geom_point(data=dot_mat2[x %in% "(1) in (2) in (3)" & y %in% intersections_upset3[[16]]], aes(x=x, y=y))+
  geom_point(data=dot_mat2[x %in% "(5) in (6)" & y %in% intersections_upset3[[17]]], aes(x=x, y=y))+
  geom_point(data=dot_mat2[x %in% "(4) in (5) in (6)" & y %in% intersections_upset3[[18]]], aes(x=x, y=y))+
  geom_point(data=dot_mat2[x %in% "(8) in (9)" & y %in% intersections_upset3[[19]]], aes(x=x, y=y))+
  geom_point(data=dot_mat2[x %in% "(7) in (8) in (9)" & y %in% intersections_upset3[[20]]], aes(x=x, y=y))+
  geom_point(data=dot_mat2[x %in% "(5) in (8)" & y %in% intersections_upset3[[21]]], aes(x=x, y=y))+
  geom_point(data=dot_mat2[x %in% "(6) in (9)" & y %in% intersections_upset3[[22]]], aes(x=x, y=y))+
  geom_point(data=dot_mat2[x %in% "(13) in (14)" & y %in% intersections_upset3[[23]]], aes(x=x, y=y))+
  
  geom_line(data=dot_mat2[x %in% "(2) in (3)" & y %in% intersections_upset3[[15]]], aes(group=1))+
  geom_line(data=dot_mat2[x %in% "(1) in (2) in (3)" & y %in% intersections_upset3[[16]]], aes(group=1))+
  geom_line(data=dot_mat2[x %in% "(5) in (6)" & y %in% intersections_upset3[[17]]], aes(group=1))+
  geom_line(data=dot_mat2[x %in% "(4) in (5) in (6)" & y %in% intersections_upset3[[18]]], aes(group=1))+
  geom_line(data=dot_mat2[x %in% "(8) in (9)" & y %in% intersections_upset3[[19]]], aes(group=1))+
  geom_line(data=dot_mat2[x %in% "(7) in (8) in (9)" & y %in% intersections_upset3[[20]]], aes(group=1))+
  geom_line(data=dot_mat2[x %in% "(5) in (8)" & y %in% intersections_upset3[[21]]], aes(group=1))+
  geom_line(data=dot_mat2[x %in% "(6) in (9)" & y %in% intersections_upset3[[22]]], aes(group=1))+
  geom_line(data=dot_mat2[x %in% "(13) in (14)" & y %in% intersections_upset3[[23]]], aes(group=1))


fill_plot <- fill_deg %>% 
  count(variable, value) %>%
  group_by(variable) %>%
  mutate(pct=prop.table(n) * 100) %>%
  ggplot() + aes(variable, pct, fill=value)+
  geom_bar(stat='identity')+
  theme_minimal()+
  theme(axis.text.x = element_text(angle = 45, size = 5, vjust = 0.865, hjust = 0.85), axis.text.y = element_text(size = 12), axis.title.y = element_text(size = 12, hjust = .5, vjust=.5), axis.title.x = element_text(size = 12), legend.text = element_text(size=12), legend.title = element_text(size=12))+
  ylab(label='Proportion of upreglated,\ndownregulated and opposite\nregulated DEG in set')+
  xlab(label = 'Gene expression comparison and intersection (in) set')+
  labs(fill = 'Genetic\nnorm')+
  geom_text(aes(label=paste0(sprintf("%1.1f", pct),"%")), size=3,
            position=position_stack(vjust=0.5))

observed_both_geno24 <- matrix(c(56, 44, 43.2, 56.8), nrow = 2, byrow = F)
# Row names and column names for better understanding of the output
colnames(observed_both_geno24) <- c("root", "leaf")
rownames(observed_both_geno24) <- c("downregulation", "upregulation")
chisq.test(observed_both_geno24)

observed_both_geno48 <- matrix(c(56, 44, 43.6, 56.4), nrow = 2, byrow = F)
# Row names and column names for better understanding of the output
colnames(observed_both_geno48) <- c("root", "leaf")
rownames(observed_both_geno48) <- c("downregulation", "upregulation")
chisq.test(observed_both_geno48) #not different

observed_between_geno_root48 <- matrix(c(51.4, 48.6, 62, 38), nrow = 2, byrow = F)
# Row names and column names for better understanding of the output
colnames(observed_between_geno_root48) <- c("TSH660", "PA121")
rownames(observed_between_geno_root48) <- c("downregulation", "upregulation")
chisq.test(observed_between_geno_root48) #not differnet

observed_between_geno_leaves48 <- matrix(c(37.8, 62.2, 44.6, 55.4), nrow = 2, byrow = F)
# Row names and column names for better understanding of the output
colnames(observed_between_geno_leaves48) <- c("TSH660", "PA211")
rownames(observed_between_geno_leaves48) <- c("downregulation", "upregulation")
chisq.test(observed_between_geno_leaves48)

# whole leaves 
observed_between_geno_leaves48_whole <- matrix(c(111.92, 68.04, 852.032, 685.948), nrow = 2, byrow = F)
# Row names and column names for better understanding of the output
colnames(observed_between_geno_leaves48_whole) <- c("leaf48_TSH660", "leaf48_PA211")
rownames(observed_between_geno_leaves48_whole) <- c("upregulation", "downregulation")
chisq.test(observed_between_geno_leaves48_whole)

observed_between_geno_leaves48_test <- matrix(c(956.636, 581.364, 852.032, 685.948), nrow = 2, byrow = F) #test with same # of DEG 4 both
# Row names and column names for better understanding of the output
colnames(observed_between_geno_leaves48_test) <- c("leaf48_TSH660", "leaf48_PA211")
rownames(observed_between_geno_leaves48_test) <- c("upregulation", "downregulation")
chisq.test(observed_between_geno_leaves48_test)

observed_between_geno_roots48 <- matrix(c(51.4, 48.6, 62, 38), nrow = 2, byrow = F)
# Row names and column names for better understanding of the output
colnames(observed_between_geno_roots48) <- c("TSH660", "PA211")
rownames(observed_between_geno_roots48) <- c("downregulation", "upregulation")
chisq.test(observed_between_geno_roots48)

# with whole numbers proportion
observed_between_geno_roots48_whole <- matrix(c(1841.4, 1127.436, 2348.6, 871.564), nrow = 2, byrow = F)
colnames(observed_between_geno_roots48_whole) <- c("root48", "leaf48")
rownames(observed_between_geno_roots48_whole) <- c("upregulation", "downregulation")
chisq.test(observed_between_geno_roots48_whole)

observed_tsh48 <- matrix(c(51.4, 48.6, 37.8, 62.2), nrow = 2, byrow = F)
colnames(observed_tsh48) <- c("root48", "leaf48")
rownames(observed_tsh48) <- c("downregulation", "upregulation")
chisq.test(observed_tsh48) #not different

#whole 
observed_tsh48_whole <- matrix(c(1164.456, 1231.544, 111.92, 68.04), nrow = 2, byrow = F)
colnames(observed_tsh48_whole) <- c("root48_tsh", "leaf48_tsh")
rownames(observed_tsh48_whole) <- c("upregulation", "downregulation")
chisq.test(observed_tsh48_whole) # different

observed_pa121 <- matrix(c(62, 38, 41.6, 58.4), nrow = 2, byrow = F)
colnames(observed_pa121) <- c("root", "leaf")
rownames(observed_pa121) <- c("downregulation", "upregulation")
chisq.test(observed_pa121)

observed_pa121_48_whole <- matrix(c(872.86, 1424.14, 852.032, 685.948), nrow = 2, byrow = F)
colnames(observed_pa121_48_whole) <- c("root48_pa121", "leaf48_pa121")
rownames(observed_pa121_48_whole) <- c("upregulation", "downregulation")
chisq.test(observed_pa121_48_whole)

#McNemar's
observed_tissues_48_whole <- matrix(c(1841.4, 2343.6, 1127.436, 871.564), nrow = 2, byrow = F)
colnames(observed_tissues_48_whole) <- c("root48_pa121", "leaf48_pa121")
rownames(observed_tissues_48_whole) <- c("upregulation", "downregulation")
mcnemar.test(observed_tissues_48_whole)

observed_geno_rt <- matrix(c(872.86, 1424.14, 1231.544, 1164.456), nrow=2, byrow = T)
observed_geno_rt <- matrix(c(873, 1424, 1232, 1164), nrow=2, byrow = T)
colnames(observed_geno_rt) <- c("rt48h_up", "rt48h_dn")
rownames(observed_geno_rt) <- c("pa121", "tsh660")
mcnemar.test(observed_geno_rt)

observed_geno_lf <- matrix(c(852, 686, 112, 68), nrow=2, byrow = T)
colnames(observed_geno_lf) <- c("lf48h_up", "lf48h_dn")
rownames(observed_geno_lf) <- c("pa121", "tsh660")
mcnemar.test(observed_geno_lf)


bar_inter <- bar_inter %>% 
  mutate(textvjust = ifelse(y<200, -.2, 1.5),
         colortext = case_when(y>200 ~ 'white', y<200 ~ 'black')
  )
bar_inter_plot <- bar_inter %>%
  ggplot(aes(x=x, y=y, label=y))+
  theme_minimal()+
  geom_bar(stat='identity')+
  ylab('Number of DEG in set')+
  theme(axis.text.x = element_blank(), axis.text.y = element_text(size = 12), axis.title.y = element_text(size=12), axis.title.x = element_blank())+
  geom_text(aes(vjust=textvjust), color=bar_inter$colortext, size=3,
            position=position_dodge(width = 1))

bar_inter_plot_PAG <- bar_inter %>%
  ggplot(aes(x=x, y=y, label=y))+
  theme_minimal()+
  geom_bar(stat='identity')+
  ylab('Number of DEG in set')+
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1), axis.text.y = element_text(size = 12), axis.title.y = element_text(size=12))+
  geom_text(aes(vjust=textvjust), color=bar_inter$colortext, size=3,
            position=position_dodge(width = 1))
ggsave('bar_inter_plot_PAG.pdf', plot = bar_inter_plot_PAG, device = 'pdf', width = 10, height = 5, dpi = 300)

#bar_inter_plot+scale_y_reverse()

pdf(file = 'upset_plot_prop_manual.pdf', width = 14, height = 10)
(fill_plot+theme(axis.text.x = element_blank(), axis.title.x = element_blank(), axis.title.y = element_text(size = 8, hjust = .5, vjust=-1)))/bar_inter_plot/dot_mat_plot+plot_layout(heights = c(6, 6, 4))
dev.off()


pdf(file = 'upset_plot_prop_manual2.pdf', width = 14, height = 10)
(dot_mat_plot+theme(axis.text.x = element_blank()))/(bar_inter_plot+scale_y_reverse())/(fill_plot+theme(axis.text.x = element_blank(), axis.title.x = element_blank(), axis.title.y = element_text(size = 8, hjust = .5, vjust=-1)))+plot_layout(heights = c(6, 6, 4))
dev.off()

pdf(file = 'upset_plot_prop_manual3.pdf', width = 14, height = 10)
(dot_mat_plot2+theme(axis.text.x = element_blank()))/(bar_inter_plot)/(fill_plot+theme(axis.text.x = element_blank(), axis.title.x = element_blank(), axis.title.y = element_text(size = 12, hjust = .5, vjust=-1)))+plot_layout(heights = c(6, 6, 4))
dev.off()

pdf(file = 'upset_plot_prop_manual4.pdf', width = 14, height = 10)
(bar_inter_plot)/(dot_mat_plot2+theme(axis.text.x = element_blank()))/(fill_plot+theme(axis.text.x = element_blank(), axis.title.x = element_blank(), axis.title.y = element_text(size = 12, hjust = .5, vjust=-1)))+plot_layout(heights = c(6, 6, 4))
dev.off()
