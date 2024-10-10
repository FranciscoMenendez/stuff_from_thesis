# Leaf Photosynthetic rate at developmental Stage
library(ggplot2)
library(magrittr)
library(dplyr)
library(data.table)

library(gridExtra)
library(broom)
library(pwr)
setwd("/Users/francisco/Documents/PSU/Lab/Cd/1project/")
lfstsu <- read.csv('LicorMeasurementDataTransfer/leaf_stage_survey_12823.csv', skip = 11, header = T ) #leaf stage survey
lfstsu <- lfstsu[2:nrow(lfstsu), ]
lfstsu <- data.table(lfstsu)
lfstsu[, (5:ncol(lfstsu)) := lapply(.SD, as.numeric), .SDcols = 5:ncol(lfstsu)]
lfstsu[, c("stage") := lapply(.SD, function(x) ifelse(x == "D", "E", x)), .SDcols = c("stage")]
lfstsu[, c("stage") := lapply(.SD, function(x) ifelse(x == "C", "D", x)), .SDcols = c("stage")]

lfstsu[, .(meanPhoto=mean(Photo), meanCond=mean(Cond),sdPhoto=sd(Photo), sdcond=sd(Cond)), .(stage) ]
t.test(lfstsu[stage=='D', Photo], lfstsu[ stage=='E', Photo])
wilcox.test(lfstsu[stage=='D', Photo], lfstsu[ stage=='E', Photo])


photo_stage_D_E <- ggplot(data = lfstsu, aes(x = stage, y = Photo))+geom_boxplot()+
  labs(x='Leaf Stage', y=expression(atop('Photosynthesis', '(umol CO2 m'^-2~'s'^-1~')')))+
  annotate('text', x = 1.5, y=2.8, label='p=1.667e-07', fontface='italic')+
  geom_segment(aes(x = 'D', xend = "E"), y=2.6, yend=2.6, color = "black", linetype = "solid")+
  geom_segment(aes(x = 'D', xend = "D"), y=2.5, yend=2.6, color = "black", linetype = "solid")+
  geom_segment(aes(x = 'E', xend = "E"), y=2.5, yend=2.6, color = "black", linetype = "solid")


t.test(lfstsu[stage=='D', Cond], lfstsu[ stage=='E', Cond])

cond_stage_D_E <- ggplot(data = lfstsu, aes(x = stage, y = Cond))+geom_boxplot()+
  labs(x='Leaf Stage', y=expression(atop('Stomatal Conductivity', '(mol H2O m'^-2~'s'^-1~')')))+
  annotate('text', x = 1.5, y=.062, label='p = 0.07608', fontface='italic')+
  geom_segment(aes(x = 'D', xend = "E"), y=0.06, yend=0.06, color = "black", linetype = "solid")+
  geom_segment(aes(x = 'D', xend = "D"), y=0.059, yend=0.06, color = "black", linetype = "solid")+
  geom_segment(aes(x = 'E', xend = "E"), y=0.059, yend=0.06, color = "black", linetype = "solid")

pdf(file = 'stage_e_d_photo_cond_boxplot.pdf', width = 7, height = 4)
photo_stage_D_E+cond_stage_D_E
dev.off()  


density_photo <- ggplot(data = lfstsu, aes(x=Photo, stat='identity', fill=stage))+
  geom_density()+
  labs(x=expression('Photosynthesis (umol CO2 m'^-2~'s'^-1~')'), y='Frequency')+
  annotate("text", x = -2, y = 0.90, 
           label = bquote(atop(bar(x) == .(round(-0.5145831, 3)), italic(s) == .(round(0.5409939, 5)))), 
           vjust = "inward", hjust = "inward")+
  annotate("text", x = 2, y = 0.90, 
           label = bquote(atop(bar(x) == .(round(1.2413025, 3)), italic(s) == .(round(0.6162215, 5)))), 
           vjust = "inward", hjust = "inward")+
  geom_segment(aes(x = -1, xend = -0.9), y=.64, yend=.625, color = "black", linetype = "solid")+
  geom_segment(aes(x = 1.2, xend = 1.1), y=.64, yend=.6, color = "black", linetype = "solid")+
  theme(legend.position = 'none')

density_cond <- ggplot(data = lfstsu, aes(x=Cond, stat='identity', fill=stage))+
  geom_density()+
  labs(x=expression('Stomatal Conductivity (mol H2O m'^-2~'s'^-1~')'), y='Frequency')+
  annotate("text", x = 0.007, y = 65, 
           label = bquote(atop(bar(x) == .(round(0.01846242, 3)), italic(s) == .(round(0.007229708, 5)))), 
           vjust = "inward", hjust = "inward")+
  annotate("text", x = 0.028, y = 65, 
           label = bquote(atop(bar(x) == .(round(0.02483390, 3)), italic(s) == .(round(0.0102472, 5)))), 
           vjust = "inward", hjust = "inward")+
  geom_segment(aes(x = 0.0125, xend = 0.0135), y=48, yend=43, color = "black", linetype = "solid")+
  geom_segment(aes(x = 0.03, xend = 0.0275), y=48, yend=43, color = "black", linetype = "solid")


pdf(file = 'stage_e_d_photo_cond_boxplot_density.pdf', width = 10, height = 8)
photo_stage_D_E+cond_stage_D_E+density_photo+density_cond+patchwork::plot_layout(ncol = 2, heights = c(2,1))+patchwork::plot_annotation(tag_levels = 'A')
dev.off()  


# CO2 Response Curve

co2flw500 <- read.csv('LicorMeasurementDataTransfer/co2 response morning auto const flow 500_.csv', skip = 11, header = T) # contatn flow 500 measured at Dec8th
co2flw500 <- co2flw500[3:nrow(co2flw500), ]
co2flw500 <- data.table(co2flw500)
co2flw500[, (5:ncol(co2flw500)) := lapply(.SD, as.numeric), .SDcols = 5:ncol(co2flw500)]
co2flw500$time <- '0950'

ggplot(data = co2flw500, aes(Ci, Photo))+
  geom_point()+
  geom_smooth()


#CO2 Response Curve Solar Noon
co2sn <- read.csv('LicorMeasurementDataTransfer/fmb-co2 response 1pm const flow_.csv', skip = 10, header = T) # contant at 1 flow 500 measured at Dec8th
co2sn <- co2sn[10:nrow(co2sn), ]
co2sn <- data.table(co2sn)
co2sn[, 5:ncol(co2sn) := lapply(.SD, as.numeric), .SDcols = 5:ncol(co2sn)]  
co2sn$time <- '1300'

ggplot(co2sn, aes(Ci, Photo))+
  geom_point()+
  geom_smooth()
  

#combining both morning and noon

cicf3 <- rbind(co2flw500, co2sn)
cicf3 <- data.table(cicf3)

co2response_photo_lm <- ggplot(cicf3, aes(Ci, Photo, color=time))+
  geom_point()+
  geom_smooth(method = 'lm')+
  theme(legend.position = 'none')+
  labs(x='Ci (uMol CO2)', y=expression('Photosynthesis (umol CO2 m'^-2~'s'^-1~')'))

co2response_cond_lm <- ggplot(cicf3, aes(Ci, Cond, color=time))+
  geom_point()+
  geom_smooth(method = 'lm')+
  theme(legend.position = 'bottom')+
  labs(x='Ci (uMol CO2)', y=expression('Stomatal Conductivity (mol H2O m'^-2~'s'^-1~')'))

co2LM <- lm(Photo~Ci, data = co2flw500)
co2LM2 <- lm(Photo~Ci+time, data = cicf3) # add the 
co2LM3 <- lm(Photo~Ci, data=cicf3)

co2_condLM <- lm(Cond ~ Ci+time, data = cicf3)
co2_condLM2 <- lm(Cond ~ Ci, data = cicf3)

co2response_photo_lm+co2response_cond_lm


# Convert lm summary to a table
co2_condLM_summary <- tidy(summary(co2_condLM))
co2_photoLM_summary <- tidy(summary(co2LM2))

# Convert table to a grob
table_grob_co2cond <- tableGrob(co2_condLM_summary)
table_grob_co2photo <- tableGrob(co2_photoLM_summary)

co2_response_lm_plot <- grid.arrange(co2response_photo_lm, table_grob_co2photo, co2response_cond_lm, table_grob_co2cond, ncol=2)


ggsave(filename =  'co2_response_lm.pdf', plot = co2_response_lm_plot, device = 'pdf', width = 14, height = 10, dpi = 300)

cicf3$wue <- cicf3$Photo/cicf3$Cond

ggplot(cicf3[Ci<400,], aes( Ci, wue)) +
  geom_point()+
  geom_smooth(method = 'lm', formula = y~log(x))

ggplot(cicf3[Ci<400,], aes( Ci, wue)) +
  geom_point()+
  geom_smooth(method = 'lm', formula = y~exp(x))


#comparing the R2 of both A~ci noon vs morning

ggplot(cicf3[Ci<400,], aes(time, Photo))+geom_violin()+geom_point()
ggplot(cicf3[Ci<400,], aes(time, Cond))+geom_violin()+geom_point()

# Fit the models for each time point
library(lme4)

model_time <- lm(Photo ~ Ci+time, data =cicf3)
model_notime <- lm(Photo ~ Ci, data = cicf3)
anova(model_time, model_notime)

model_timegs <- lm(Cond ~ Ci+time, data =cicf3)
model_notimegs <- lm(Cond ~ Ci, data = cicf3)
anova(model_timegs, model_notimegs)

# Fit a combined model including the interaction term
combined_model <- lm(Photo ~ Ci * time, data = cicf3)




# data from barroso

barr <- read.csv('../LicorMeasurementDataTransfer/barroso_sup_data.csv') #check the data source 
barr <- data.table(barr)
ggplot(data = barr[mn==0], aes(x =  cd, y=photo, color=factor(day)))+
  geom_point()+
  geom_line()
  

barr$wue2 <- NA

barr[, wue2 := as.numeric(wue2)]
barr[!is.na(gs), wue2 := photo/gs]

barr[is.na(gs), gs2 := photo/wue]

barr[, gs3 := ifelse(cd == 28, gs, gs2)]
barr[, gs3 := fifelse(cd == 28, gs, gs2)]
barr[day == 28, gs3 := gs] 


ggplot(data = barr[mn==0], aes(x =  cd, y=photo, color=factor(day)))+
  geom_point()+
  geom_line()



photosynthesis_barr2023 <- ggplot(data = barr[mn==0, ], aes(x = factor(cd), y = photo, fill = factor(day))) +
  geom_bar(stat = "identity", position = "dodge") +
  labs(x = "Cd Treatment Level", y = "Photosynthesis Value", fill = "Days After Treatment") +
  theme_minimal()

gs_barr2023 <- ggplot(data = barr[mn==0 & !is.na(gs2), ], aes(x = factor(cd), y = gs2, fill = factor(day))) +
  geom_bar(stat = "identity", position = "dodge") +
  labs(x = "Cd Treatment Level", y = "Gs (calculated)", fill = "Days After Treatment") +
  theme_minimal()

gsboth_barr2023 <- ggplot(data = barr[mn==0 & !is.na(gs3), ], aes(x = factor(cd), y = gs3, fill = factor(day))) +
  geom_bar(stat = "identity", position = "dodge") +
  labs(x = "Cd Treatment Level", y = "Gs (calculated day 1-14,\n day 28 reported)", fill = "Days After Treatment") +
  theme_minimal()

pdf(file = 'barr_data.pdf', width = 7, height = 7 )
photosynthesis_barr2023/gs_barr2023/gsboth_barr2023
dev.off()

#wide photo
wide_photo <- dcast(barr[mn==0, ], day ~ cd, value.var = "photo")
wide_photo[, rowMean := rowMeans(.SD, na.rm = TRUE), .SDcols = c("0", "0.2", "0.4", "0.6", "0.8")]
photo_t <- mean(wide_photo$rowMean) #overall mean
wide_photo[, phmt_mi2 := (photo_t-rowMean)^2]

#wide photo
wide_photo <- dcast(barr[mn==0, ], cd ~ day, value.var = "photo")
wide_photo[, rowMean := rowMeans(.SD, na.rm = TRUE), .SDcols = c("1", "4", "7", "14","21", "28")]
photo_t <- mean(wide_photo$rowMean) #overall mean
wide_photo[, phmt_mi2 := (photo_t-rowMean)^2]

sd_m_photo <- sqrt(sum(wide_photo$phmt_mi2)/6)
f_photo <- sd_m_photo/sd(as.matrix(wide_photo[, 2:7]))

#power test with f_photo
pwr.anova.test(k = 3, f = f_photo, sig.level = 0.05, power = 0.8) #effect size calculated using the gs sd from leaf stage e of samples

#n=11

#wide gs
wide_gs <- dcast(barr[mn==0 & !is.na(gs3), ], cd ~ day, value.var = "gs3")
wide_gs[, rowMean := rowMeans(.SD, na.rm = TRUE), .SDcols = c("1", "4", "7", "14", "28")]
gs_t <- mean(wide_gs$rowMean) #overall mean
wide_gs[, phmt_mi2 := (gs_t-rowMean)^2]

sd_m_gs <- sqrt(sum(wide_gs$phmt_mi2)/6)
f_gs <- sd_m_gs/sd(as.matrix(wide_gs[, 2:7]))


#power test with f_gs
pwr.anova.test(k = 3, f = f_gs, sig.level = 0.05, power = 0.8) #effect size calculated using the gs sd from leaf stage e of samples

#n=16


#Plot power curves

sample_sizes <- seq(from = 3, to = 40, by = 2)

powers_gs <- sapply(sample_sizes, function(n) {
  pwr.anova.test(n = n, f = f_gs, k = 3, sig.level = 0.05)$power
})
powers_photo <- sapply(sample_sizes, function(n) {
  pwr.anova.test(n = n, f = f_photo, k = 3, sig.level = 0.05)$power
})


pdf("power_curves_photo_gs.pdf", width = 7, height = 5)
par(mfrow = c(1,2))
plot(sample_sizes, powers_gs, type = "l", xlab = "Group Size", ylab = "Power", main = 'Gs Power Curve')
abline(h = .80, v = 16.26815)
text(x = 10, y=0.9,  'n~16' )
segments(x0 = 10, y0 = .88, x1 = 16.26815, y1 = 0.8)

plot(sample_sizes, powers_photo, type = "l", xlab = "Group Size", ylab = "Power", main = "A Power Curve")
abline(h = .80, v = 10.8963)
text(x = 22, y=0.9,  'n~11' )
segments(x0 = 20, y0 = .88, x1 = 10.8963, y1 = 0.8)

dev.off()

write.table(wide_gs, file = 'barroso_data_cd_only.txt', quote = F, sep = '\t', row.names = F, col.names = T)
write.table(wide_photo, file = 'barroso_data_cd_only.txt', quote = F, sep = '\t', row.names = F, col.names = T, append = T)


#Since the n=16 might be too high we will discard the three level. 
#Test prediction of gs reduction based on a lm of the gs of Barroso and predict
#it to level 10ppm. 
#Then use that to estimate the effect and propose a 2 group t.test experiment. 

lmgsbarre <- lm(formula = gs3~exp(cd), data = barr[mn==0 &!is.na(gs3) &day!=28, ])
lmgsbarr <- lm(formula = gs3~cd, data = barr[mn==0 &!is.na(gs3) &day!=28, ])
summary(lmgsbarr)

predicted <- predict(lmgsbarr, newdata = data.frame(cd=seq(from=0.8, to=0.01, by=-0.01)))
seq(from=0.8, to=0.01, by=-0.01)
abline(lmgsbarr)
plot(y=barr[mn==0 &day!=28, gs3], x=barr[mn==0 & day!=28, cd], xlab = 'Cd level', ylab = 'gs3')
lines(seq(from=0.8, to=0.01, by=-0.01), predicted)
pdf()

lmphotobarr <- lm(formula = photo~cd, data = barr[mn==0 &!is.na(gs3) &day!=28, ])
summary(lmphotobarr)

predictedphoto <- predict(lmphotobarr, newdata = data.frame(cd=seq(from=0.8, to=0.01, by=-0.01)))
plot(y=barr[mn==0 &day!=28, photo], x=barr[mn==0 & day!=28, cd], xlab = 'Cd level', ylab = 'Photsynthesis')
lines(seq(from=0.8, to=0.01, by=-0.01), predictedphoto)

predictedgs <- data.frame(cd=seq(from=0.8, to=0.01, by=-0.01), predicted)

average_table <- barr[, .(average_photo = mean(photo, na.rm = TRUE),
                        average_gs = mean(gs, na.rm = TRUE)), by = cd]
average_table[, phmt_mi2 := (mean(average_photo)-average_photo)^2] #total mean - mean squared
average_table[, gsmt_mi2 := (mean(average_gs)-average_gs)^2] 

sd_m_photo <- sqrt(sum(average_table$phmt_mi2)/nrow(average_table))
f_photo <- sd_m_photo/sd(average_table$average_photo)

sd_m_gs <- sqrt(sum(average_table$gsmt_mi2)/nrow(average_table))
f_gs <- sd_m_gs/sd(average_table$average_gs)




# power estimation for sample size using t test

pwr.t.test(d = 1.653016, sig.level = 0.05, power = 0.80, type = 'one', alternative = 'greater')
pwr.t.test(d = 1.653016, sig.level = 0.05, power = 0.80, type = 'two', alternative = 'greater')
pwr.t.test(d = 1.653016, sig.level = 0.05, power = 0.80, type = 'pair', alternative = 'greater')

pwr.t.test(d = 1.653016, sig.level = 0.05, power = 0.80, type = 'two', alternative = 'two.sided')

# A powr
pwr.t.test(d = 1.503749, sig.level = 0.05, power = 0.80, type = 'two', alternative = 'two.sided')

# using CICf3 data

pwr.t.test(d = 0.4775756, sig.level = 0.05, power = 0.80, type = 'two', alternative = 'two.sided')
pwr.t.test(d = 0.5237367, sig.level = 0.05, power = 0.80, type = 'two', alternative = 'two.sided')

# leaf stage

pwr.t.test(d = 1.954281, sig.level = 0.05, power = 0.80, type = 'two', alternative = 'two.sided')
pwr.t.test(d = 1.387696, sig.level = 0.05, power = 0.80, type = 'two', alternative = 'two.sided')



# power anova
# first calculate the F table using the leaf stage data
stagelm <- lm(Cond~stage, data = lfstsu)
anova(stagelm)

sqrt(0.00024982/0.00204550) # sum of srqs over error sum of sqrs (residuals)
#[1] 0.3494732
pwr.anova.test(k = 3, f = 0.3494732, sig.level = 0.05, power = 0.8) #effect size calculated using the gs sd from leaf stage e of samples





stageplm <- lm(Photo~stage, data = lfstsu)
anova(stageplm)
sqrt(18.973/8.330)
#[1] 1.509196
pwr.anova.test(k = 3, f = 1.509196, sig.level = 0.05, power = 0.8) #effect size calculated using the gs sd from leaf stage e of samples

#using the stage e cohen d divided by 2

pwr.anova.test(k = 3, f = 1.653016, sig.level = 0.05, power = 0.8)

pwr.anova.test(k = 3, f = 0.5237367, sig.level = 0.05, power = 0.8) #effet size calculated using the gs sd from cicf3

pwr.anova.test(k = 3, f = 1.387696/2, sig.level = 0.05, power = 0.8) #effect size calculated using the gs sd from leaf stage e of samples

sqrt(0.00024982/0.00204550)
#[1] 0.3494732

pwr.anova.test(k = 3, f = 0.3494732, sig.level = 0.05, power = 0.8) #effect size calculated using the gs sd from leaf stage e of samples

pwr.anova.test(k = 3, f = f_photo, sig.level = 0.05, power = 0.8) #effect size calculated using the gs sd from leaf stage e of samples

pwr.anova.test(k = 3, f = f_gs, sig.level = 0.05, power = 0.8) #effect size calculated using the gs sd from leaf stage e of samples






#CO2 Response Curve in quarentine room (No intracellular CO2??)
co2qua <- read.csv('LicorMeasurementDataTransfer/fmb-co2 response 1pm const flow  quarent_2.csv', skip = 11, header = T) # contant at 1 flow 500 measured at Dec8th
co2qua <- co2qua[3:nrow(co2qua), ]
co2qua <- data.table(co2qua)
co2qua[, 5:ncol(co2qua) := lapply(.SD, as.numeric), .SDcols = 5:ncol(co2qua)]

ggplot(co2qua, aes( Photo, Ci))+
  geom_point()+
  labs(y='Ci (uMol CO2)', x=expression('Photosynthesis (umol CO2 m'^-2~'s'^-1~')'))
ggplot(co2qua, aes( Ci, Photo))+
  geom_point()+
  labs(x='Ci (uMol CO2)', y=expression('Photosynthesis (umol CO2 m'^-2~'s'^-1~')'))



plantnumtb <- read.table('LicorMeasurementDataTransfer/PLANT_TRATMENT_NUMBER_TREVB.txt', skip = 2, fill = T)
pntb <- plantnumtb[, 1:3]
pnta <- plantnumtb[, 4:6]
names(pnta) <- pnta[1, ]
names(pntb) <- pnta[1, ]
pnta <- pnta[2:22, ]
pntb  <-  pntb[2:28, ]

pnt <- rbind(pntb, pnta)

write.table(pnt[nrow(pnt):1, ], file = 'plant_number_tray_trevb.txt', append = F, quote = F, sep = '\t', row.names = F, col.names = T)

read.csv('../Cd_treatment_exp/plant_number_lvl.csv')
              


#scratch
# Extract the sum of squared residuals and degrees of freedom
# first model made, later replaced with the likelihood test ratio
SSR_early <- sum(residuals(model_early)^2)
SSR_late <- sum(residuals(model_late)^2)
df_early <- df.residual(model_early)
df_late <- df.residual(model_late)

# Calculate the F-statistic
F_statistic <- ((SSR_early - SSR_late) / (df_early - df_late)) / (SSR_late / df_late)

# Calculate the p-value
p_value <- pf(F_statistic, df_early - df_late, df_late, lower.tail = FALSE)

# Output the results
cat("F-statistic:", F_statistic, "\n")
cat("p-value:", p_value, "\n")
