# hormone charts
library(patchwork)
# preliminary data for the hormones

ho <- read.csv('../Cd_treatment_exp/Phytohormones_Burns_082024.csv', na.strings = 'N/A')

ho <- data.table(ho)
ho[, aba_m:=mean(ABA), .(tissue, trt, time)]
ho[, cinnamic_m:=mean(cinnamic.acid), .(tissue, trt, time)]
ho[, benzoic_m:=mean(benzoic.acid), .(tissue, trt, time)]
ho[, iaa_m:=mean(IAA, na.rm=T), .(tissue, trt, time)]
ho[, ja_m:=mean(JA, na.rm=T), .(tissue, trt, time)]
ho[, oxd_m:=mean(oxd.acid), .(tissue, trt, time)]
ho[, sa_m:=mean(SA), .(tissue, trt, time)]
ho[, aba_sd:=sd(ABA), .(tissue, trt, time)]
ho[, aba_sem:=sd(ABA)/sqrt(.N), .(tissue, trt, time)]
ho[, ja_sd:=sd(JA), .(tissue, trt, time)]
ho[, ja_sem:=sd(JA)/sqrt(.N), .(tissue, trt, time)]
ho[, iaa_sd:=sd(IAA), .(tissue, trt, time)]
ho[, iaa_sem:=sd(IAA)/sqrt(.N), .(tissue, trt, time)]
ho[, sa_sd:=sd(SA), .(tissue, trt, time)]
ho[, sa_sem:=sd(SA)/sqrt(.N), .(tissue, trt, time)]
ho[, cinnamic_sd:=sd(cinnamic_m), .(tissue, trt, time)]
ho[, cinnamic_sem:=sd(cinnamic_m)/sqrt(.N), .(tissue, trt, time)]
ho[, benzoic_sd:=sd(benzoic.acid), .(tissue, trt, time)]
ho[, benzoic_sem:=sd(benzoic.acid)/sqrt(.N), .(tissue, trt, time)]

hol <- melt(data = ho, id.vars = c('metabo_core', 'pn', 'tissue', 'trt', 'time'), measure.vars = c('aba_m', 'cinnamic_m', 'benzoic_m', 'iaa_m', 'ja_m', 'oxd_m', 'sa_m', 'aba_sem'))

ho_rt <- ggplot(data = hol[tissue=='root'], aes(x = factor(time), y = value, fill=trt))+
  geom_bar(stat='identity', position = position_dodge(width = 0.9))+
  xlab('time after Cd 12PPM')+
  ylab('nMol')+
  ggtitle('Roots hormones')+
  guides(fill='none')+
  facet_wrap(~variable, scales = 'free_y', ncol = 1)

ho_lf <- ggplot(data = hol[tissue=='e2'], aes(x = factor(time), y = value, fill=trt))+
  geom_bar(stat='identity', position = position_dodge(width = 0.9))+
  xlab('time after Cd 12PPM')+
  ylab(element_blank())+
  ggtitle('Leaves hormones')+
  facet_wrap(~variable, scales = 'free_y', ncol = 1)

rt_lf_ho <- ho_rt+ho_lf+plot_layout(ncol = 2)
ggsave(filename = 'root_leaves_hormones_cd_prelim.pdf', plot = rt_lf_ho, device = 'pdf', width = 5, height = 10, dpi = 300)

#just significant hormones

t.test(ho[tissue=='root' & time==36 & trt=='cd+', ABA], ho[tissue=='root' & time==36 & trt=='cd-', ABA], alternative = 'g')
t.test(ho[tissue=='e2' & time==36 & trt=='cd+', ABA], ho[tissue=='e2' & time==36 & trt=='cd-', ABA], alternative = 'g')

root_aba <- ggplot(data = ho[tissue=='root' & time==36], aes(x = factor(time), y = aba_m, fill=trt))+
  geom_bar(stat='identity', position = position_dodge(width = 0.9), color='black', linewidth=0.2)+
  geom_point(aes(y=ABA), position = position_jitterdodge(jitter.width = 0.1, dodge.width = 0.9), alpha=0.5)+
  geom_errorbar(aes(ymin = aba_m - aba_sem, ymax = aba_m + aba_sem), width = 0.2,
                position = position_dodge(width = 0.9)) +
  annotate('text', x = '36', y=3.1, label='1.95 greater than 0.21 nMol, p=0.026', fontface='italic', size=2.5)+
  geom_segment(aes(x = 0.75, xend = 1.25), y = 3, yend = 3, color = "black", linetype = "solid", linewidth = 0.5) +
  geom_segment(aes(x = 0.75, xend = 0.75), y = 2.95, yend = 3.05, color = "black", linetype = "solid", linewidth = 0.5) +
  geom_segment(aes(x = 1.25, xend = 1.25), y = 2.95, yend = 3.05, color = "black", linetype = "solid", linewidth = 0.5) +
  xlab(element_blank())+
  ylab('nMol')+
  ggtitle('Root ABA at 36h')+
  theme_light()+
  theme(legend.position = 'inside', legend.position.inside = c(0.9, 0.8),legend.key.size = unit(0.4, 'cm'),legend.title = element_blank(),  # Adjust position as needed (x, y)
        legend.background = element_rect(fill = alpha('white', 0.5))) +
  scale_x_discrete(labels = c("36" = "36h"))+
  scale_fill_discrete(name='Cd treatments', type =  c('white', 'grey'))


leaf_aba <- ggplot(data = ho[tissue=='e2' & time==36], aes(x = factor(time), y = aba_m, fill=trt))+
  geom_bar(stat='identity', position = position_dodge(width = 0.9), color='black', linewidth=0.2)+
  geom_point(aes(y=ABA), position = position_jitterdodge(jitter.width = 0.1, dodge.width = 0.9), alpha=0.5)+
  geom_errorbar(aes(ymin = aba_m - aba_sem, ymax = aba_m + aba_sem), width = 0.2,
                position = position_dodge(width = 0.9)) +
  annotate('text', x = '36', y=47, label='28.6 no greater than 25.1 nMol, p = 0.4', fontface='italic', size=2.5)+
  geom_segment(aes(x = 0.75, xend = 1.25), y = 45, yend = 45, color = "black", linetype = "solid", linewidth = 0.5) +
  geom_segment(aes(x = 0.75, xend = 0.75), y = 44, yend = 46, color = "black", linetype = "solid", linewidth = 0.5) +
  geom_segment(aes(x = 1.25, xend = 1.25), y = 44, yend = 46, color = "black", linetype = "solid", linewidth = 0.5) +
  xlab('Time after Cd 12 PPM')+
  ylab('nMol')+
  ggtitle('Leaf ABA at 36h')+
  theme_light()+
  theme(legend.position = 'inside', legend.position.inside = c(0.9, 0.8),legend.key.size = unit(0.4, 'cm'),legend.title = element_blank(),  # Adjust position as needed (x, y)
        legend.background = element_rect(fill = alpha('white', 0.5))) +  # Add semi-transparent background
  scale_x_discrete(labels = c("36" = "36h"))+
  scale_fill_discrete(name='Cd treatments', type =  c('white', 'grey'))

root_aba+leaf_aba+plot_layout(ncol=1)

iaa_root <- ggplot(data = ho[tissue=='root' & time==36], aes(x = factor(time), y = iaa_m, fill=trt))+
  geom_bar(stat='identity', position = position_dodge(width = 0.9), color='black', linewidth=0.2)+
  geom_point(aes(y=IAA), position = position_jitterdodge(jitter.width = 0.1, dodge.width = 0.9), alpha=0.5)+
  geom_errorbar(aes(ymin = iaa_m - iaa_sem, ymax = iaa_m + iaa_sem), width = 0.2,
                position = position_dodge(width = 0.9)) +
  xlab('time after Cd 12PPM')+
  ylab('nMol')+
  ggtitle('Root IAA')+
  theme_light()+
  scale_fill_discrete(name='Cd treatments', type =  c('white', 'grey'))+
  guides(fill='none')

iaa_leaf <- ggplot(data = ho[tissue=='e2' & time==36], aes(x = factor(time), y = iaa_m, fill=trt))+
  geom_bar(stat='identity', position = position_dodge(width = 0.9), color='black', linewidth=0.2)+
  geom_point(aes(y=IAA), position = position_jitterdodge(jitter.width = 0.1, dodge.width = 0.9), alpha=0.5)+
  geom_errorbar(aes(ymin = iaa_m - iaa_sem, ymax = iaa_m + iaa_sem), width = 0.2,
                position = position_dodge(width = 0.9)) +
  xlab('time after Cd 12PPM')+
  ylab('nMol')+
  ggtitle('Root IAA')+
  theme_light()+
  scale_fill_discrete(name='Cd treatments', type =  c('white', 'grey'))+
  guides(fill='none')

iaa_root+iaa_leaf

ja_root <- ggplot(data = ho[tissue=='root' & time==36], aes(x = factor(time), y = ja_m, fill=trt))+
  geom_bar(stat='identity', position = position_dodge(width = 0.9), color='black', linewidth=0.2)+
  geom_point(aes(y=JA), position = position_jitterdodge(jitter.width = 0.1, dodge.width = 0.9), alpha=0.5)+
  geom_errorbar(aes(ymin = ja_m - ja_sem, ymax = ja_m + ja_sem), width = 0.2,
                position = position_dodge(width = 0.9)) +
  xlab('time after Cd 12PPM')+
  ylab('nMol')+
  ggtitle('Root JA')+
  theme_light()+
  scale_fill_discrete(name='Cd treatments', type =  c('white', 'grey'))+
  guides(fill='none')

ja_leaf <- ggplot(data = ho[tissue=='e2' & time==36], aes(x = factor(time), y = ja_m, fill=trt))+
  geom_bar(stat='identity', position = position_dodge(width = 0.9), color='black', linewidth=0.2)+
  geom_point(aes(y=JA), position = position_jitterdodge(jitter.width = 0.1, dodge.width = 0.9), alpha=0.5)+
  geom_errorbar(aes(ymin = ja_m - ja_sem, ymax = ja_m + ja_sem), width = 0.2,
                position = position_dodge(width = 0.9)) +
  xlab('time after Cd 12PPM')+
  ylab('nMol')+
  ggtitle('Leaf JA')+
  theme_light()+
  scale_fill_discrete(name='Cd treatments', type =  c('white', 'grey'))+
  guides(fill='none')

ja_root+ja_leaf

sa_root <- ggplot(data = ho[tissue=='root' & time==36], aes(x = factor(time), y = sa_m, fill=trt))+
  geom_bar(stat='identity', position = position_dodge(width = 0.9), color='black', linewidth=0.2)+
  geom_point(aes(y=SA), position = position_jitterdodge(jitter.width = 0.1, dodge.width = 0.9), alpha=0.5)+
  geom_errorbar(aes(ymin = sa_m - sa_sem, ymax = sa_m + sa_sem), width = 0.2,
                position = position_dodge(width = 0.9)) +
  xlab('time after Cd 12PPM')+
  ylab('nMol')+
  ggtitle('Root SA')+
  theme_light()+
  scale_fill_discrete(name='Cd treatments', type =  c('white', 'grey'))+
  guides(fill='none')

sa_leaf <- ggplot(data = ho[tissue=='e2' & time==36], aes(x = factor(time), y = sa_m, fill=trt))+
  geom_bar(stat='identity', position = position_dodge(width = 0.9), color='black', linewidth=0.2)+
  geom_point(aes(y=SA), position = position_jitterdodge(jitter.width = 0.1, dodge.width = 0.9), alpha=0.5)+
  geom_errorbar(aes(ymin = sa_m - sa_sem, ymax = sa_m + sa_sem), width = 0.2,
                position = position_dodge(width = 0.9)) +
  xlab('time after Cd 12PPM')+
  ylab('nMol')+
  ggtitle('Leaf SA')+
  theme_light()+
  scale_fill_discrete(name='Cd treatments', type =  c('white', 'grey'))+
  guides(fill='none')

sa_root+sa_leaf

benzoic_root <- ggplot(data = ho[tissue=='root' & time==36], aes(x = factor(time), y = benzoic_m, fill=trt))+
  geom_bar(stat='identity', position = position_dodge(width = 0.9), color='black', linewidth=0.2)+
  geom_point(aes(y=benzoic.acid), position = position_jitterdodge(jitter.width = 0.1, dodge.width = 0.9), alpha=0.5)+
  geom_errorbar(aes(ymin = benzoic_m - benzoic_sem, ymax = benzoic_m + benzoic_sem), width = 0.2,
                position = position_dodge(width = 0.9)) +
  xlab('time after Cd 12PPM')+
  ylab('nMol')+
  ggtitle('Root Benzoic acid')+
  theme_light()+
  scale_fill_discrete(name='Cd treatments', type =  c('white', 'grey'))+
  guides(fill='none')

benzoic_leaf <- ggplot(data = ho[tissue=='e2' & time==36], aes(x = factor(time), y = benzoic_m, fill=trt))+
  geom_bar(stat='identity', position = position_dodge(width = 0.9), color='black', linewidth=0.2)+
  geom_point(aes(y=benzoic.acid), position = position_jitterdodge(jitter.width = 0.1, dodge.width = 0.9), alpha=0.5)+
  geom_errorbar(aes(ymin = benzoic_m - benzoic_sem, ymax =benzoic_m + benzoic_sem), width = 0.2,
                position = position_dodge(width = 0.9)) +
  xlab('time after Cd 12PPM')+
  ylab('nMol')+
  ggtitle('Leaf Benzoic Acid')+
  theme_light()+
  scale_fill_discrete(name='Cd treatments', type =  c('white', 'grey'))+
  guides(fill='none')
benzoic_root+benzoic_leaf

load('../Cd_treatment_exp/dec1415_first_Cdtrt_run.rda')
load('../Cd_treatment_exp/plots_licor_dec15.rda')


cond_dec15_plot <- cond_dec15_plot+xlab(element_blank())+ggtitle('Gs at 36h')

cond_date_plot20 <- cond_date_plot20+ylab(element_blank())+xlab(element_blank())+ggtitle('Gs across 2 days')

photo_date_plot20 <- photo_date_plot20+ylab(element_blank())+ggtitle('A across 2 days')

photo_dec15_plot <- photo_dec15_plot+ggtitle('A at 36h')


figure7_cond_hormones <- cond_dec15_plot+cond_date_plot20+root_aba+photo_dec15_plot+photo_date_plot20+leaf_aba+patchwork::plot_layout(ncol = 3, axes='keep')+patchwork::plot_annotation(tag_levels = list(c('A', 'C', 'E', 'B', 'D', 'F')) )

ggsave(filename = 'figure7_cond_hormones.pdf', plot = figure7_cond_hormones, device = 'pdf', width = 10.5, height = 6,dpi = 300)


