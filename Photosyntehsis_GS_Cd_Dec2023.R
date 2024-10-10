# Photosynthesis and GS after Cd treatment

# read data

# plant treatment level will be assinged based on plant number, which was recorded at the moment of measurment. 

# libraries 
library(ggplot2)
library(data.table)
library(lme4)

# set working directory

setwd('../Cd_treatment_exp/')

trt <- read.csv('plant_number_lvl.csv')
trt <- data.table(trt)
trt[,bioreps := paste(Tray, lvl, sep='')]
trt[, .(plants = paste(Plant, collapse = ', '), .N), .(bioreps)]


# function

read_licor_mean <- function(file_name, date, alle=F){
  dt <- read.csv(file_name)
  dt <- dt[2:nrow(dt), ]
  dt <- data.table(dt)
  dt[, code := paste(plant.number, positon, Cd.level, sep='')]
  dt[, bioreps := paste(positon, Cd.level, sep='')]
  
  
  dt[, .(plants = paste(plant.number, collapse = ", "),
             count = .N), 
         by = bioreps]
  
  trt[match(dt[, plant.number], Plant), .(Plant, lvl)]
  
  dt[match(trt[, Plant], plant.number),]
  
  print(data.frame(correct_lvl=trt[match(dt[, plant.number], Plant), lvl], 
             correct_plant=trt[match(dt[, plant.number], Plant), Plant],
             written_lvl=dt$Cd.level, 
             written_plant=dt$plant.number))
  
  #corrected treatment level
  
  dt$Cd.level <- trt[match(dt[, plant.number], Plant), lvl]
  dt$positon <- trt[match(dt[, plant.number], Plant), Tray]
  
  dt[, code := paste(plant.number, positon, Cd.level, sep='')]
  dt[, bioreps := paste(positon, Cd.level, sep='')]
  dt[, .(plants = paste(plant.number, collapse = ", "),
             count = .N), 
         by = bioreps]
  
  # AVERAGE BY PLANT CODE 
  #Each plant should have one observation, per data policy we will keep the 
  #average between repeat measurements
  
  dt[, .N, .(code)]
  
  dt[, 9:(ncol(dt)-2) := lapply(.SD, as.numeric), .SDcols = 9:(ncol(dt)-2)]
  
  dt[, photo_m := mean(Photo), .(code)]
  dt[, cond_m := mean(Cond), .(code)]
  dt[, ci_m := mean(Ci), .(code)]
  dt[, wue := photo_m/cond_m]
  dt[, vpdl_m := mean(VpdL), .(code)]
  dt[, rhs_m := mean(RH_S), .(code)]
  dt[, tleaf_m := mean(Tleaf), .(code)]
  dt[, tair_m := mean(Tair), .(code)]
  dt[, rh_r_m := mean(RH_R), .(code)]
  dt[, pari_m :=mean(PARi), .(code)]
  dt[, ci_photo := photo_m/ci_m, .(code)]
  dt[, ci_cond := cond_m/ci_m, .(code)]
  
  
  #dt[, Cd.level := factor(Cd.level, levels = c('0', '10', '20'), ordered = T)]
  #dt[, Tray := factor(positon, levels = c('A', 'B', 'C', 'D', 'E', 'F', 'G', 'H'), ordered = T)]
  
  # Split the time into hours, minutes, and seconds
  dt[, c("hours", "minutes", "seconds") := tstrsplit(HHMMSS, ":", type.convert = TRUE)]
  
  # Convert hours and minutes to seconds and add them up
  dt[, secondsOD := hours * 3600 + minutes * 60 + seconds]
  
  dt[secondsOD<20000, 'secondsOD'] <- dt[ secondsOD<20000, 'secondsOD']+(12*3600)
  
  dt[, hours := as.numeric(hours)]
  dt[, minutes := as.numeric(minutes)]
  dt[, secondsOD := as.numeric(secondsOD)]
  
  dt[, hours := mean(hours), .(code)]
  dt[, minutes := mean(minutes), .(code)]
  dt[, secondsOD := mean(secondsOD), .(code)]
  
  #convert to factors Cd.level
  
  dt$Cd.level <- factor(dt$Cd.level)
  
  if (alle==T){
    dt <- unique(dt[, .(code, hours, minutes, secondsOD, positon, plant.number, Cd.level, photo_m, cond_m, wue, ci_m, ci_photo, ci_cond, vpdl_m, rhs_m, tleaf_m, tair_m, rh_r_m, pari_m)]) # 48 observations
  }else{
    dt <- unique(dt[, .(code, hours, minutes, secondsOD, positon, plant.number, Cd.level, photo_m, cond_m, wue, ci_m)]) # 48 observations    
  }
  
  
  dt$date <- date
  
  return(dt)
  
}




#proceed to update the function to include a more general approach replacing also the positon

dec14am <- read_licor_mean(file_name = 'leaf cd dec14_.csv', date = 'dec14am', alle = T)

dec14pm <- read_licor_mean(file_name = 'leaf Cd dec14 1pm_.csv', date = 'dec14pm', alle = T) # this now produced only 48 observations

dec15am <- read_licor_mean(file_name = 'leaf cd survey dec15 morning_.csv', date = 'dec15am', alle = T)

dec15pm <- read_licor_mean('leaf Cd survey dec15 1pm_.csv', date = 'dec15pm', alle = T)

dec1415 <- do.call(rbind, list(dec14am, dec14pm, dec15am, dec15pm))

save(dec1415, file = 'dec1415_first_Cdtrt_run.rda')

dec1415[ci_m<0, ci_m := NA]

# plot the dist of each parameter by date and treatment

ggplot(data = dec1415[!is.na(ci_m),], aes(ci_m, stat='identity', fill=Cd.level))+
  geom_density()+
  facet_wrap(~date)

ggplot(data = dec1415[,], aes(vpdl_m, stat='identity', fill=Cd.level))+
  geom_density()+
  facet_wrap(~date)

ggplot(data = dec1415[,], aes(pari_m, stat='identity', fill=Cd.level))+
  geom_density()+
  facet_wrap(~date)

ggplot(data = dec1415[,], aes(tleaf_m, stat='identity', fill=Cd.level))+
  geom_density()+
  facet_wrap(~date)

ggplot(data = dec1415[,], aes(rhs_m, stat='identity', fill=Cd.level))+
  geom_density()+
  facet_wrap(~date)

ggplot(data = dec1415[,], aes(ci_cond, stat='identity', fill=Cd.level))+
  geom_density()+
  facet_wrap(~date)

ggplot(data = dec1415[,], aes(ci_photo, stat='identity', fill=Cd.level))+
  geom_density()+
  facet_wrap(~date)

# all seem to be independent of the treatment

# now we can see if any correlation occurs between the predictors and the cond_m or photo_m

#Ci and photo cond
ggplot(data = dec1415[!is.na(ci_m),], aes(x=ci_m, y=photo_m, color=Cd.level))+
  geom_point()+
  facet_wrap(~date) #there seems to be some variaton in Ci that might be confounding the Photo
ggplot(data = dec1415[!is.na(ci_m),], aes(x=ci_m, y=cond_m, color=Cd.level))+
  geom_point()+
  facet_wrap(~date)

#Ci and ci normalized photo and cond
ggplot(data = dec1415[!is.na(ci_m),], aes(x=ci_m, y=ci_cond, color=Cd.level))+
  geom_point()+
  facet_wrap(~date)
ggplot(data = dec1415[!is.na(ci_m),], aes(x=ci_m, y=ci_photo, color=Cd.level))+
  geom_point()+
  facet_wrap(~date)

#temperature leaf and photo and cond
ggplot(data = dec1415[!is.na(tleaf_m),], aes(x=tleaf_m, y=cond_m, color=Cd.level))+
  geom_point()+
  facet_wrap(~date)
ggplot(data = dec1415[!is.na(tleaf_m),], aes(x=tleaf_m, y=photo_m, color=Cd.level))+
  geom_point()+
  facet_wrap(~date)

#vapor pressure deficit and cond and photo
ggplot(data = dec1415[!is.na(vpdl_m),], aes(x=vpdl_m, y=photo_m, color=Cd.level))+
  geom_point()+
  facet_wrap(~date)
ggplot(data = dec1415[!is.na(vpdl_m),], aes(x=vpdl_m, y=cond_m, color=Cd.level))+
  geom_point()+
  facet_wrap(~date)

#par in and photo and cond
ggplot(data = dec1415[!is.na(pari_m),], aes(x=pari_m, y=photo_m, color=Cd.level))+
  geom_point()+
  facet_wrap(~date)
ggplot(data = dec1415[!is.na(pari_m),], aes(x=pari_m, y=cond_m, color=Cd.level))+
  geom_point()+
  facet_wrap(~date)

#relative humidity in sample and photo and cond
ggplot(data = dec1415[!is.na(rhs_m),], aes(x=rhs_m, y=cond_m, color=Cd.level))+
  geom_point()+
  facet_wrap(~date)
ggplot(data = dec1415[!is.na(rhs_m),], aes(x=rhs_m, y=photo_m, color=Cd.level))+
  geom_point()+
  facet_wrap(~date)

# cond and photo
ggplot(data = dec1415[,], aes(x=cond_m, y=photo_m, color=Cd.level))+
  geom_point()+
  facet_wrap(~date)

#WUE
ggplot(data = dec1415[wue<200,], aes(x=Cd.level, y=wue, color=Cd.level))+
  geom_point()+
  geom_boxplot()+
  facet_wrap(~date)


# Linear models by day for photo_m

dec14am_photolm <- lm(photo_m ~ Cd.level + positon, data = dec14am)
dec14pm_photolm <- lm(photo_m ~ Cd.level + positon, data = dec14pm)

dec15am_photolm <- lm(photo_m ~ Cd.level + positon, data = dec15am)
dec15pm_photolm <- lm(photo_m ~ Cd.level + positon, data = dec15pm)

dec_photolm <- lm(photo_m ~ date, data = dec1415)
dec_photo0lm <- lm(photo_m ~ date, data = dec1415[Cd.level=='0'])
dec_photo10lm <- lm(photo_m ~ date, data = dec1415[Cd.level=='10'])
dec_photo20lm <- lm(photo_m ~ date, data = dec1415[Cd.level=='20'])

# linear model for ci normalized photo
dec14am_ciphotolm <- lm(ci_photo ~ Cd.level + positon, data = dec14am)
dec14pm_ciphotolm <- lm(ci_photo ~ Cd.level + positon, data = dec14pm)

dec15am_ciphotolm <- lm(ci_photo ~ Cd.level + positon, data = dec15am)
dec15pm_ciphotolm <- lm(ci_photo ~ Cd.level + positon, data = dec15pm)

# Linear models by day for cond_m

dec14am_condlm <- lm(cond_m ~ Cd.level + positon, data = dec14am)
dec14pm_condlm <- lm(cond_m ~ Cd.level + positon, data = dec14pm)

dec15am_condlm <- lm(cond_m ~ Cd.level + positon, data = dec15am)
dec15pm_condlm <- lm(cond_m ~ Cd.level + positon, data = dec15pm)

summary(dec15pm_condlm) # this shows significant effect of Cd.level ***

dec_cond20lm <- lm(cond_m ~ date, data = dec1415[Cd.level=='20'])
dec_cond10lm <- lm(cond_m ~ date, data = dec1415[Cd.level=='10'])
dec_cond0lm <- lm(cond_m ~ date, data = dec1415[Cd.level=='0'])
dec_condlm <- lm(cond_m ~ date, data = dec1415)

dec15pm_condlmr <- lmer(cond_m ~ positon + (1|Cd.level), data = dec1415)
red_dec15pm_lmr <- lmer(cond_m ~ (1|Cd.level), data = dec1415) #reduced model

anova(red_dec15pm_lmr, dec15pm_condlmr)
anova(dec15pm_condlmr, red_dec15pm_lmr)

# Linear models by day for ci normalized cond_m

dec14am_cicondlm <- lm(ci_cond ~ Cd.level + positon, data = dec14am)
dec14pm_cicondlm <- lm(ci_cond ~ Cd.level + positon, data = dec14pm)

dec15am_cicondlm <- lm(ci_cond ~ Cd.level + positon, data = dec15am)
dec15pm_cicondlm <- lm(ci_cond ~ Cd.level + positon, data = dec15pm)

# Linear models by day for wue

dec14am_wuelm <- lm(wue ~ Cd.level + positon, data = dec14am)
dec14pm_wuelm <- lm(wue ~ Cd.level + positon, data = dec14pm)

dec15am_wuelm <- lm(wue ~ Cd.level + positon, data = dec15am)
dec15pm_wuelm <- lm(wue ~ Cd.level + positon, data = dec15pm)



# ANOVA of allllll
photo_lm <- list(dec14am_photolm, dec14pm_photolm, dec15am_photolm, dec15pm_photolm, dec_photo0lm, dec_photo10lm, dec_photo20lm, dec_photolm)
cond_lm <- list(dec14am_condlm, dec14pm_condlm, dec15am_condlm, dec15pm_condlm, dec_cond0lm, dec_cond10lm, dec_cond20lm, dec_condlm)
cicond_lm <- list(dec14am_cicondlm, dec14pm_cicondlm, dec15am_cicondlm, dec15pm_cicondlm)
ciphoto_lm <- list(dec14am_ciphotolm, dec14pm_ciphotolm, dec15am_ciphotolm, dec15pm_ciphotolm)
wue_lm <- list(dec14am_wuelm, dec14pm_wuelm, dec15am_wuelm, dec15pm_wuelm)

photo_anova <- lapply(photo_lm, anova)
names(photo_anova) <- c('dec14am_photolm', 'dec14pm_photolm', 'dec15am_photolm', 'dec15pm_photolm',  'dec_photo0lm', 'dec_photo10lm', 'dec_photo20lm', 'dec_photolm')

cond_anova <- lapply(cond_lm, anova)
names(cond_anova) <- c('dec14am_condlm', 'dec14pm_condlm', 'dec15am_condlm', 'dec15pm_condlm', 'dec_cond0lm', 'dec_cond10lm', 'dec_cond20lm', 'dec_condlm')

ciphoto_anova <- lapply(ciphoto_lm, anova)
names(ciphoto_anova) <- c('dec14am_ciphotolm', 'dec14pm_ciphotolm', 'dec15am_ciphotolm', 'dec15pm_ciphotolm')

cicond_anova <- lapply(cicond_lm, anova)
names(cicond_anova) <- c('dec14am_cicondlm', 'dec14pm_cicondlm', 'dec15am_cicondlm', 'dec15pm_cicondlm')

wue_anova <- lapply(wue_lm, anova) 
names(wue_anova) <- c('dec14am_wuelm', 'dec14pm_wuelm', 'dec15am_wuelm', 'dec15pm_wuelm')


# Convert lm summary to a table
cond_anova_table <- lapply(cond_anova, tidy)
cond_anova_grob <- lapply(cond_anova_table, tableGrob)
cond_anova_grob[[1]] <- arrangeGrob(cond_anova_grob[[1]], top = "Dec 14 morning, Conductivity")
cond_anova_grob[[2]] <- arrangeGrob(cond_anova_grob[[2]], top = "Dec 14 afternoon, Conductivity")
cond_anova_grob[[3]] <- arrangeGrob(cond_anova_grob[[3]], top = "Dec 15 morning, Conductivity")
cond_anova_grob[[4]] <- arrangeGrob(cond_anova_grob[[4]], top = "Dec 15 afternoon, Conductivity")

grid.arrange(cond_anova_grob[[1]], cond_anova_grob[[2]], cond_anova_grob[[3]], cond_anova_grob[[4]], ncol=1)

photo_anova_table <- lapply(photo_anova, tidy)
photo_anova_grob <- lapply(photo_anova_table, tableGrob)
photo_anova_grob[[1]] <- arrangeGrob(photo_anova_grob[[1]], top = 'Dec 14 morning, Photosynthesis')
photo_anova_grob[[2]] <- arrangeGrob(photo_anova_grob[[2]], top = 'Dec 14 afternoon, Photosynthesis')
photo_anova_grob[[3]] <- arrangeGrob(photo_anova_grob[[3]], top = 'Dec 15 morning, Photosynthesis')
photo_anova_grob[[4]] <- arrangeGrob(photo_anova_grob[[4]], top = 'Dec 15 afternoon, Photosynthesis')

wue_anova_table <- lapply(wue_anova, tidy)
wue_anova_grob <- lapply(wue_anova_table, tableGrob)
wue_anova_grob[[1]] <- arrangeGrob(wue_anova_grob[[1]], top = 'Dec 14 morning WUE')
wue_anova_grob[[2]] <- arrangeGrob(wue_anova_grob[[2]], top = 'Dec 14 afternoon WUE')
wue_anova_grob[[3]] <- arrangeGrob(wue_anova_grob[[3]], top = 'Dec 15 morning WUE')
wue_anova_grob[[4]] <- arrangeGrob(wue_anova_grob[[4]], top = 'Dec 15 afternoon WUE')

grid.arrange(photo_anova_grob[[1]], photo_anova_grob[[2]], photo_anova_grob[[3]], photo_anova_grob[[4]],
             cond_anova_grob[[1]], cond_anova_grob[[2]], cond_anova_grob[[3]], cond_anova_grob[[4]],
             ncol=2)

photo_tables <- grid.arrange(photo_anova_grob[[1]], photo_anova_grob[[2]], photo_anova_grob[[3]], photo_anova_grob[[4]])
cond_tables <- grid.arrange(cond_anova_grob[[1]], cond_anova_grob[[2]], cond_anova_grob[[3]], cond_anova_grob[[4]])
wue_tables <- grid.arrange(wue_anova_grob[[1]],wue_anova_grob[[2]], wue_anova_grob[[3]], wue_anova_grob[[4]])



co2_response_lm_plot <- grid.arrange(co2response_photo_lm, table_grob_co2photo, co2response_cond_lm, table_grob_co2cond, ncol=2)


# Convert table to a grob
table_grob_co2cond <- tableGrob(co2_condLM_summary)
table_grob_co2photo <- tableGrob(co2_photoLM_summary)

#linear model with interaction term
#  for cond_m

dec14amInt_lm <- lm(cond_m ~ Cd.level + positon + Cd.level*positon, data = dec14am)
dec14pmInt_lm <- lm(cond_m ~ Cd.level + positon + Cd.level*positon, data = dec14pm)

dec15amInt_lm <- lm(cond_m ~ Cd.level + positon + Cd.level*positon, data = dec15am)
dec15pmInt_lm <- lm(cond_m ~ Cd.level + positon + Cd.level*positon, data = dec15pm)

anova(dec15pmInt_lm)
summary(dec15pmInt_lm)

#imputation for tray f

imputation_dec15pm <- dec15pm[, .(mean_photo_m = mean(photo_m, na.rm = TRUE),
                                  mean_cond_m = mean(cond_m, na.rm = TRUE),
                                  mean_wue = mean(wue, na.rm = TRUE), 
                                  mean_ci = mean(ci_m, na.rm = TRUE)),
                              by = Cd.level]


dec15pmIMP <- rbind(dec15pm, data.table(code=c('31F10', '35F20', '36F0', '40G10'), #impute plants 31F10, 35F20, 36F0
                                     hours=c(0.6044792, 0.6055671, 0.6069560, 0.6156713), 
                                     minutes=c(0.6044792, 0.6055671, 0.6069560, 0.6156713), 
                                     secondsOD=c(45413.00, 45416.98, 45422.07, 45453.97), 
                                     positon=c('F', 'F', 'F', "G"), 
                                     plant.number=c(31, 35, 36, 40), 
                                     Cd.level=c(10, 20, 0, 10), 
                                     photo_m=c(1.462827, 1.750753, 1.681702, 1.462827), 
                                     cond_m=c(0.03466625, 0.02793819, 0.03980343, 0.03466625), 
                                     wue=c(49.29941, 130.66710, 45.96496, 49.29941), 
                                     ci_m=c(406.3976, 271.7678, 411.6007, 406.3976), 
                                     date=c('dec15pm', 'dec15pm', 'dec15pm', 'dec15pm')))

dec15pmIMP_condlm <- lm(cond_m ~ Cd.level + positon, data = dec15pmIMP)

# using a liner model to predit the values (omit data to avoid overfitting)

lm_cond_predMV <- lm(cond_m~Cd.level+positon, data = dec1415) #linear model
                                                              #used to predict 
                                                              #the missing values, 
                                                              #model leaves the date out
                                                              #to avoid overfitting

#using a linear model with all environmental var to predict

lm_cond_palle <- lm(cond_m) #predicted all environment


predicted_f_g_MV <- predict(object = lm_cond_predMV, newdata = dec15pmIMP)

dec15pmIMP$cond_pred <- predicted_f_g_MV # prediction based on model using all dates

residuals_dec15imp <- data.table(code=dec15pmIMP$code, fitted= dec15pmIMP$cond_pred, residuals=dec15pmIMP$cond_m-dec15pmIMP$cond_pred) # residuals 

ggplot(residuals_dec15imp, aes(x = fitted, y = residuals)) +
  geom_point() +  # Adds the points
  geom_hline(yintercept = 0, linetype = "dashed") +  # Adds a horizontal line at y = 0
  xlab("Fitted Values") + 
  ylab("Residuals") +
  ggtitle("Residuals vs Fitted Values")


ggplot(residuals_dec15imp, aes(residuals))+
  geom_point()

# Photosynthesis quick plotting 
ggplot(dec14am, aes(x = factor(Cd.level), y = photo_m, group=factor(Cd.level)))+
  geom_boxplot()

ggplot(dec14pm, aes(x = Cd.level, y = photo_m, group=Cd.level))+
  geom_boxplot()

ggplot(dec15am, aes(x = Cd.level, y = photo_m, group=Cd.level))+
  geom_boxplot()

ggplot(dec15pm, aes(x = Cd.level, y = photo_m, group=Cd.level))+
  geom_boxplot()

photo_dec15_plot <- ggplot(dec1415[date %in% c('dec15pm'),], aes(x = Cd.level, y = photo_m, group=factor(Cd.level)))+
  geom_boxplot(size=0.25, outlier.shape = NA)+
  geom_point(position = position_jitter(width = 0.2), alpha = 0.2) +
  scale_x_discrete(labels = c("0"="0 PPM", "10" = "6 PPM", "20" = "12 PPM"))+
  labs(x="Cd Treatment Level (PPM)", y=expression(atop('Photosynthesis', paste( '(A) (', mu, "mol CO"[2]* " m"^-2~"s"^-1, ")"))))+
  theme_light()

# Cond
ggplot(dec14am, aes(x = Cd.level, y = cond_m, group=Cd.level))+
  geom_boxplot()

ggplot(dec14pm, aes(x = Cd.level, y = cond_m, group=Cd.level))+
  geom_boxplot()

ggplot(dec15am, aes(x = Cd.level, y = cond_m, group=Cd.level))+
  geom_boxplot()

# Post Hoc Test Tukey's 

dec15pm$Cd.level <- factor(dec15pm$Cd.level)
dec15pm$positon <- factor(dec15pm$positon)

dec15pm_condaov <- aov(formula = cond_m ~ Cd.level + positon, data = dec15pm)
dec15pm_condposth <- TukeyHSD(dec15pm_condaov)
pdf('post_hoc_dec15pm.pdf', width = 7, height = 20)
plot(dec15pm_condposth)
dev.off()

dec_condaov <- aov(formula = cond_m ~ date, data = dec1415)
dec_condaov <- TukeyHSD(dec_condaov)

dec_condaov20 <- aov(formula = cond_m ~ date, data = dec1415[Cd.level=='20'])
dec_condaov20 <- TukeyHSD(dec_condaov20)

dec_condaov0 <- aov(formula = cond_m ~ date, data = dec1415[Cd.level=='0'])
dec_condaov0 <- TukeyHSD(dec_condaov0)

## This is in the text
dec15pm_20ppm <- dec15pm[Cd.level==20, ]
dec14am_20ppm <- dec14am[Cd.level==20, ]
t.test(dec15pm_20ppm[, cond_m], dec14am_20ppm[, cond_m])
t.test(dec15pm_20ppm[, cond_m], dec14am_20ppm[, cond_m], alternative = 'less')
t.test(dec15pm_20ppm[, cond_m], dec14am_20ppm[, cond_m], alternative = 'g')

dec15pm_0ppm <- dec15pm[Cd.level==0, ]
dec14am_0ppm <- dec14am[Cd.level==0, ]
t.test(dec15pm_0ppm[, cond_m], dec14am_0ppm[, cond_m])
t.test(dec15pm_0ppm[, cond_m], dec14am_0ppm[, cond_m], alternative = 'less')

dec15pm_10ppm <- dec15pm[Cd.level==10, ]
dec14am_10ppm <- dec14am[Cd.level==10, ]
t.test(dec15pm_10ppm[, cond_m], dec14am_10ppm[, cond_m])
t.test(dec15pm_10ppm[, cond_m], dec14am_10ppm[, cond_m], alternative = 'less')



#pdf('post_hoc_dec15pm.pdf', width = 7, height = 20)
#plot(dec15pm_condposth)
#dev.off()
dec_photoaov <- aov(formula = photo_m ~ date, data = dec1415)
dec_photoaov <- TukeyHSD(dec_photoaov)

dec_photoaov20 <- aov(formula = photo_m ~ date, data = dec1415[Cd.level=='20'])
dec_photoaov20 <- TukeyHSD(dec_photoaov20)

dec_photoaov0 <- aov(formula = photo_m ~ date, data = dec1415[Cd.level=='0'])
dec_photoaov0 <- TukeyHSD(dec_photoaov0)

#plot by date
scale_x_discrete(labels = c("dec14am" = "Day 1 AM", "dec14pm" = "Day 1 PM", 
                            "dec15am" = "Day 2 AM", "dec15pm" = "Day 2 PM"))


photo_date_plot20 <- ggplot(dec1415[Cd.level=='20',], aes(x = date, y = photo_m, group=factor(date)))+
  geom_boxplot(size=0.25, outlier.shape = NA)+
  geom_point(position = position_jitter(width = 0.2), alpha = 0.2) +
  scale_x_discrete(labels = c("dec14am" = "Day 1\nAM", "dec14pm" = "Day 1\nPM", 
                              "dec15am" = "Day 2\nAM", "dec15pm" = "Day 2\nPM"))+
  labs(x="Time of measurement at 12 PPM Cd", y=expression(atop('Photosynthesis', paste( '(A) (', mu, "mol CO"[2]* " m"^-2~"s"^-1, ")"))))+
  theme_light()

photo_date_plot <- ggplot(dec1415[,], aes(x = date, y = photo_m, group=factor(date)))+
  geom_boxplot(size=0.25, outlier.shape = NA)+
  geom_point(position = position_jitter(width = 0.2), alpha = 0.2) +
  scale_x_discrete(labels = c("dec14am" = "Day 1\nAM", "dec14pm" = "Day 1\nPM", 
                              "dec15am" = "Day 2\nAM", "dec15pm" = "Day 2\nPM"))+
  labs(x="Time of measurement at 12 PPM Cd", y=expression(atop('Photosynthesis', paste( '(', mu, "mol CO"[2]* " m"^-2~"s"^-1, ")"))))+
  facet_wrap(~Cd.level, nrow=3, scale = 'free')+
  theme_light()
  
  annotate('text', x = 'dec14am', y = 5, label='a')+
  annotate('text', x = 'dec14pm', y = 5, label='b')+
  annotate('text', x = 'dec15am', y = 5, label='abc')+
  annotate('text', x = 'dec15pm', y = 5, label='bcd')+
  
  
  annotate('text', x = 2.0, y=4.9, label='diff=-0.56662320, padj=0.0252164', fontface='italic', size=2.5)+
  geom_segment(aes(x = 'dec14am', xend = "dec14pm"), y=4.6, yend=4.6, color = "black", linetype = "solid", linewidth=.25)+
  geom_segment(aes(x = 'dec14am', xend = "dec14am"), y=4.5, yend=4.7, color = "black", linetype = "solid", linewidth=.25, alpha=0.5)+
  geom_segment(aes(x = 'dec14pm', xend = "dec14pm"), y=4.5, yend=4.7, color = "black", linetype = "solid", linewidth=.25, alpha=0.5)+
  annotate('text', x = 2.5, y=5.5, label='diff=-0.64832314, padj=0.0091594', fontface='italic', size=2.5)+
  geom_segment(aes(x = 'dec14am', xend = "dec15pm"), y=5.2, yend=5.2, color = "black", linetype = "solid", linewidth=.25, alpha=0.5)+
  geom_segment(aes(x = 'dec14am', xend = "dec14am"), y=5.1, yend=5.3, color = "black", linetype = "solid", linewidth=.25, alpha=0.5)+
  geom_segment(aes(x = 'dec15pm', xend = "dec15pm"), y=5.1, yend=5.3, color = "black", linetype = "solid", linewidth=.25, alpha=0.5)+
  theme_light()

cond_date_plot20 <- ggplot(dec1415[Cd.level=="20",], aes(x = date, y = cond_m, group=factor(date)))+
  geom_boxplot(size=0.25, outlier.shape = NA)+
  geom_point(position = position_jitter(width = 0.2), alpha = 0.2) +
  scale_x_discrete(labels = c("dec14am" = "Day 1\nAM", "dec14pm" = "Day 1\nPM", 
                              "dec15am" = "Day 2\nAM", "dec15pm" = "Day 2\nPM"))+
  labs(x="Time of measurement at 12 PPM Cd", y=expression(atop('Stomatal Conductance', '(Gs) (mol H'[2]*'O m'^-2~'s'^-1~')')))+
  annotate('text', x = 'dec14am', y = 0.125, label='a')+
  annotate('text', x = 'dec14pm', y = 0.125, label='ab')+
  annotate('text', x = 'dec15am', y = 0.125, label='ab')+
  annotate('text', x = 'dec15pm', y = 0.125, label='b')+
  theme_light()

cond_date_plot <- ggplot(dec1415[,], aes(x = date, y = cond_m, group=factor(date)))+
  geom_boxplot(size=0.25, outlier.shape = NA)+
  geom_point(position = position_jitter(width = 0.2), alpha = 0.2) +
  scale_x_discrete(labels = c("dec14am" = "Day 1\nAM", "dec14pm" = "Day 1\nPM", 
                              "dec15am" = "Day 2\nAM", "dec15pm" = "Day 2\nPM"))+
  labs(x="Time of measurement at 12 PPM Cd", y=expression(atop('Stomatal Conductance', '(Gs) (mol H'[2]*'O m'^-2~'s'^-1~')')))+
  facet_wrap(~Cd.level, nrow = 3, scales = 'free')+
  theme_light()
##

cond_dec15_plot <- ggplot(dec1415[date %in% c('dec15pm'),], aes(x = Cd.level, y = cond_m, group=factor(Cd.level)))+
  geom_boxplot(size=0.25, outlier.shape = NA)+
  geom_point(position = position_jitter(width = 0.2), alpha = 0.2) +
  scale_x_discrete(labels = c("0"="0 PPM", "10" = "6 PPM", "20" = "12 PPM"))+
  labs(x="Cd Treatment Level (PPM)",  y=expression(atop('Stomatal Conductance', '(Gs) (mol H'[2]*'O m'^-2~'s'^-1~')')))+
  annotate('text', x = '20', y = 0.075, label='b')+
  annotate('text', x = '0', y = 0.075, label='a')+
  annotate('text', x = '10', y = 0.075, label='ab')+
  theme_light()



ggsave('cond_dec15_plot_summaryRes', plot = cond_dec15_plot, device = 'pdf', width = 7, height = 5, units = 'cm', dpi = 'screen')

res_summary_plots <- cond_dec15_plot+photo_dec15_plot



ggsave('graphical_plot_summaryRes', plot = res_summary_plots, device = 'pdf', width = 15, height = 5, units = 'cm', dpi = 'screen')

res_summary_plots2 <- cond_dec15_plot+cond_date_plot+photo_dec15_plot+photo_date_plot+patchwork::plot_layout(ncol = 2, axes = 'collect')+patchwork::plot_annotation(tag_levels = 'A')
ggsave('graphical_plot_summaryRes2', plot = res_summary_plots2, device = 'pdf', width = 15, height = 5, units = 'cm', dpi = 'screen')

res_summary_plots3 <- cond_dec15_plot+cond_date_plot20+photo_dec15_plot+photo_date_plot20+patchwork::plot_layout(ncol = 2, axes='collect')+patchwork::plot_annotation(tag_levels = list(c('A', 'C', 'B', 'D')) )
save(file = 'plots_licor_dec15.rda', list = c('cond_dec15_plot','cond_date_plot20', 'photo_dec15_plot', 'photo_date_plot20'))

ggsave('graphical_plot_summaryRes3.pdf', plot = res_summary_plots3, device = 'pdf', width = 20, height = 12, units = 'cm', dpi = 'screen')

# WUE
ggplot(dec14am, aes(x = Cd.level, y = wue, group=Cd.level))+
  geom_boxplot()

ggplot(dec14pm, aes(x = Cd.level, y = wue, group=Cd.level))+
  geom_boxplot()

ggplot(dec15am, aes(x = Cd.level, y = wue, group=Cd.level))+
  geom_boxplot()

ggplot(dec15pm, aes(x = Cd.level, y = wue, group=Cd.level))+
  geom_boxplot()


#all together


ggplot(dec1415, aes(x = Cd.level, y = cond_m, group=Cd.level))+
  geom_boxplot()
ggplot(dec1415, aes(x = Cd.level, y = photo_m, group=Cd.level))+
  geom_boxplot()
ggplot(dec1415, aes(x = Cd.level, y = wue, group=Cd.level))+
  geom_boxplot()

A_alldays <- ggplot(dec1415, aes(x = Cd.level, y = photo_m, group=Cd.level))+
  geom_boxplot()+
  labs(x="Cd Treatment Level (PPM)", y=expression('Photosynthesis (umol CO2 m'^-2~'s'^-1~')'))+
  theme(axis.title.x = element_blank())+
  facet_wrap(~date)

Gs_alldays <- ggplot(dec1415, aes(x = Cd.level, y = cond_m, group=Cd.level))+
  geom_boxplot()+
  labs(x="Cd Treatment Level (PPM)", y=expression('Stomatal Conductance (Gs) (mol H2O m'^-2~'s'^-1~')'))+
  theme(axis.title.x = element_blank())+
  facet_wrap(~date, scales = 'free_y')

WUE_alldays <- ggplot(dec1415, aes(x = Cd.level, y = wue, group=Cd.level))+
  geom_boxplot()+
  labs(x="Cd Treatment Level (PPM)", y=expression('WUE (umol CO2*mol H2O '^-1~')'))+
  facet_wrap(~date, scale = 'free_y')

boxplot_all_daysA_Gs_WUE <- A_alldays+/Gs_alldays/WUE_alldays
ggsave('boxplot_all_daysA_Gs_WUE', plot = boxplot_all_daysA_Gs_WUE, device = 'pdf', width = 10, height = 12, units = 'in', dpi = 300)

boxplot_all_daysA_Gs_WUE.ANOVA_tables <- (A_alldays+photo_tables)/(Gs_alldays+cond_tables)/(WUE_alldays+wue_tables)
ggsave('boxplot_all_daysA_Gs_WUE.ANOVA_tables.pdf', plot = boxplot_all_daysA_Gs_WUE.ANOVA_tables, device = 'pdf', width = 24, height = 21, units = 'in', dpi = 300)

# distribution of values

ggplot(dec1415, aes(x=photo_m, stat='identity'))+
  geom_density()
ggplot(dec14am, aes(x=photo_m, stat='identity'))+
  geom_density()
ggplot(dec14pm, aes(x=photo_m, stat='identity'))+
  geom_density()
ggplot(dec15am, aes(x=photo_m, stat='identity'))+
  geom_density()
ggplot(dec15pm, aes(x=photo_m, stat='identity'))+
  geom_density()

ggplot(dec1415, aes(x=cond_m, stat='identity'))+
  geom_density()
ggplot(dec14am, aes(x=cond_m, stat='identity'))+
  geom_density()
ggplot(dec14pm, aes(x=cond_m, stat='identity'))+
  geom_density()
ggplot(dec15am, aes(x=cond_m, stat='identity'))+
  geom_density()
ggplot(dec15pm, aes(x=cond_m, stat='identity'))+
  geom_density()

qqplot(x=rnorm(n = 48), y = dec1415$photo_m)
qqplot(x=rnorm(n = 48), y = dec1415$cond_m)

ggplot(dec1415[], aes(x = Cd.level, y=photo_m, color=factor(positon)))+
  geom_point()+
  geom_smooth(method = 'lm')

lm_photo_all <- lm(photo_m~Cd.level+positon+date, data = dec1415)  #evidence for a strong positional effect 
lm_cond_all <- lm(cond_m~Cd.level+positon+date, data = dec1415)  #evidence for a strong positional effect 



ggplot(data = dec1415[date=='dec15pm', ], aes(x = positon, y=photo_m))+
  geom_point()

ggplot(data = dec1415[date=='dec14pm', ], aes(x = positon, y=photo_m))+
  geom_point() #measurements made later in the day have lower Photosynthesis values 

ggplot(data = dec1415[date=='dec15pm', ], aes(x = positon, y=cond_m))+
  geom_point()

ggplot(data = dec1415[date=='dec14pm', ], aes(x = positon, y=cond_m))+
  geom_point()

lm_photo_posA <- lm(photo_m~Cd.level+date, data = dec1415[positon=='A'])  #evidence for a strong positional effect 
summary(lm_photo_posA)

lm_photo_posH <- lm(photo_m~Cd.level+date, data = dec1415[positon=='H'])  #evidence for a strong positional effect 
summary(lm_photo_posH)

# using time as seconds of the day

ggplot(data = dec1415[ grep('15am', date),], aes(x = secondsOD, y=photo_m, color=factor(Cd.level)))+
  geom_point()

ggplot(data = dec1415[ grep('15pm', date),], aes(x = secondsOD, y=photo_m, color=factor(Cd.level)))+
  geom_point()

ggplot(data = dec1415[ grep('14am', date),], aes(x = secondsOD, y=photo_m, color=factor(Cd.level)))+
  geom_point()

# Plot cond across days 
pdf('time_course_a_gs.pdf', width = 7, height = 10)
ggplot(data = dec1415, aes(date, cond_m))+geom_point()+geom_boxplot()+
ggplot(data = dec1415, aes(date, photo_m))+geom_point()+geom_boxplot()+patchwork::plot_layout(ncol = 2)+patchwork::plot_annotation('A')
dev.off()

ggplot(data = dec15pm, aes(positon, cond_m))+geom_point()+geom_boxplot()

dec1415[, cf1 := paste0(positon, Cd.level)] #not. really neccessary since ~wrap 

ggplot(data = dec1415[date=='dec15pm' & positon=='F'], aes(cf1, cond_m))+geom_point()+geom_boxplot()

ggplot(data = dec1415[date=='dec15pm' ], aes(Cd.level, cond_m))+geom_point()+geom_boxplot()+facet_wrap(~positon)
ggplot(data = dec1415[ ], aes(Cd.level, cond_m))+geom_point()+geom_boxplot()+facet_wrap(~positon)

pdf('confounding_cond.pdf', width = 14, height = 10)
ggplot(data = dec1415[positon=='A'], aes(date, cond_m))+geom_point()+geom_boxplot()+facet_wrap(~Cd.level)+
ggplot(data = dec1415[positon=='B'], aes(date, cond_m))+geom_point()+geom_boxplot()+facet_wrap(~Cd.level)+
ggplot(data = dec1415[positon=='C'], aes(date, cond_m))+geom_point()+geom_boxplot()+facet_wrap(~Cd.level)+
ggplot(data = dec1415[positon=='D'], aes(date, cond_m))+geom_point()+geom_boxplot()+facet_wrap(~Cd.level)+
ggplot(data = dec1415[positon=='E'], aes(date, cond_m))+geom_point()+geom_boxplot()+facet_wrap(~Cd.level)+
ggplot(data = dec1415[positon=='F'], aes(date, cond_m))+geom_point()+geom_boxplot()+facet_wrap(~Cd.level)+
ggplot(data = dec1415[positon=='H'], aes(date, cond_m))+geom_point()+geom_boxplot()+facet_wrap(~Cd.level)+
ggplot(data = dec1415[positon=='G'], aes(date, cond_m))+geom_point()+geom_boxplot()+facet_wrap(~Cd.level)+patchwork::plot_layout(ncol = 2)+patchwork::plot_annotation(tag_levels = 'A')
dev.off()

# since position has a significant effect on cond, 
# this means that will see that any cond values from positon G (***)
# will be significantly different form the rest, likewise positon F (*)

pdf('positon_effect_dec15pm.pdf', width = 7, height = 5)
ggplot(data = dec15pm, aes(x = positon, y=cond_m))+geom _point()+geom_boxplot()
dev.off()

#Using all dates
ggplot(data = dec1415, aes(x = positon, y=cond_m))+geom_point()+geom_boxplot()

dec15pm_condaov <- aov(formula = cond_m ~ Cd.level + positon, data = dec15pm)
dec15pm_condposth <- TukeyHSD(dec15pm_condaov)
plot(dec15pm_condposth)

# effect size 
fval <- summary(dec15pm_condaov)[[1]]$'F value'[[1]]
n <- nrow(dec15pm)
k <- 3
df_total <- n-k
cohens_f <- sqrt(fval / (fval + df_total))
library(pwr)
pwr.anova.test(k = 3, f = cohens_f, sig.level = 0.05, power = 0.8) #effect size calculated using the gs sd from leaf stage e of samples

wide_dt[, colMean := colMeans(.SD, na.rm = T), .SDcols = c()]

#wide cond
wide_cond <- t(wide_dt)
wide_cond <- wide_cond[2:nrow(wide_cond), ]
colnames(wide_cond) <- paste('p', 1:ncol(wide_cond), sep='')
wide_cond <- data.table(wide_cond)
wide_cond[, rowMean := rowMeans(.SD, na.rm = TRUE), .SDcols = paste('p', 1:ncol(wide_cond), sep='')]
cond_t <- mean(wide_cond$rowMean)

wide_cond[, condmt_mi2 := (cond_t-rowMean)^2]

sd_m_cond <- sqrt(sum(wide_cond$condmt_mi2)/16) #divided bt the number of reps
f_cond <- sd_m_cond/sd(as.matrix(wide_cond[, 1:15]))

#power test with f_photo
pwr.anova.test(k = 3, f = f_cond, sig.level = 0.05, power = 0.8) #effect size calculated using the gs sd from leaf stage e of samples

# this makes a ridiculous group size, ill change to only using cohens d between 0 and 20...

mean1 <- mean(wide_dt$level0, na.rm = TRUE)
mean2 <- mean(wide_dt$level20, na.rm = TRUE)

# Calculate pooled standard deviation
sd_pool <- sqrt(((var(wide_dt$level0, na.rm = TRUE) * (nrow(wide_dt[!is.na(level0)]) - 1)) +
                   (var(wide_dt$level10, na.rm = TRUE) * (nrow(wide_dt[!is.na(level10)]) - 1))) / 
                  (nrow(wide_dt[!is.na(level0)]) + nrow(wide_dt[!is.na(level20)]) - 2))

# Calculate Cohen's d
cohens_d <- (mean1 - mean2) / sd_pool
pwr.t.test(d = cohens_d, sig.level = 0.05, power = 0.80, type = 'one', alternative = 'greater')


write.table(x = wide_cond, file = 'wide_cond_dec15pm.txt', append = F, quote = F, sep = '\t', row.names = F, col.names = T)

sample_sizes <- seq(from = 3, to = 40, by = 2)

powers_gs <- sapply(sample_sizes, function(n) {
  pwr.t.test(d = cohens_d, sig.level = 0.05, n = n, type = 'one', alternative = 'greater')$power
})

pdf(file = 'gs_power_curve_second_exp.pdf', width = 5, height = 7)
plot(sample_sizes, powers_gs, type = "l", xlab = "Group Size", ylab = "Power", main = 'Gs Power Curve')
abline(h = .80, v = 9.5)
text(x = 5, y=0.9,  'n~10' )
segments(x0 = 9.5, y0 = .8, x1 = 5, y1 = 0.88)
dev.off()

#scratch


dec14m <- read.csv('leaf cd dec14_.csv')
dec14m <- dec14m[2:nrow(dec14m), ]
dec14m <- data.table(dec14m)
dec14m[, code := paste(plant.number, positon, Cd.level, sep='')]
dec14m[, bioreps := paste(positon, Cd.level, sep='')]


dec14m[, .(plants = paste(plant.number, collapse = ", "),
           count = .N), 
       by = bioreps]

trt[match(dec14m[, plant.number], Plant), .(Plant, lvl)]

dec14m[match(trt[, Plant], plant.number),]

data.frame(correct_lvl=trt[match(dec14m[, plant.number], Plant), lvl], 
           correct_plant=trt[match(dec14m[, plant.number], Plant), Plant],
           written_lvl=dec14m$Cd.level, 
           written_plant=dec14m$plant.number)

#corrected treatment level

dec14m$Cd.level <- trt[match(dec14m[, plant.number], Plant), lvl]


dec14m[, code := paste(plant.number, positon, Cd.level, sep='')]
dec14m[, bioreps := paste(positon, Cd.level, sep='')]
dec14m[, .(plants = paste(plant.number, collapse = ", "),
           count = .N), 
       by = bioreps]

# AVERAGE BY PLANT CODE 
#Each plant should have one observation, per data policy we will keep the 
#average between repeat measurements

dec14m[, .N, .(code)]

dec14m[, 9:(ncol(dec14m)-2) := lapply(.SD, as.numeric), .SDcols = 9:(ncol(dec14m)-2)]

dec14m[, photo_m := mean(Photo), .(code)]
dec14m[, cond_m := mean(Cond), .(code)]
dec14m[, ci_m := mean(Ci), .(code)]

dec14m <- unique(dec14m[, .(code, positon, plant.number, Cd.level, photo_m, cond_m, ci_m)]) # 48 observations


# manual for dec14 1pm

dt <- read.csv('leaf Cd dec14 1pm_.csv')
dt <- dt[2:nrow(dt), ]
dt <- data.table(dt)
dt[, code := paste(plant.number, positon, Cd.level, sep='')]
dt[, bioreps := paste(positon, Cd.level, sep='')]


dt[, .(plants = paste(plant.number, collapse = ", "),
       count = .N), 
   by = bioreps]

trt[match(dt[, plant.number], Plant), .(Plant, lvl)]


print(data.frame(correct_lvl=trt[match(dt[, plant.number], Plant), lvl], 
                 correct_plant=trt[match(dt[, plant.number], Plant), Plant],
                 correct_tray=trt[match(dt[, plant.number], Plant), Tray],
                 written_lvl=dt$Cd.level, 
                 written_plant=dt$plant.number, 
                 written_tray=dt$positon))

#corrected treatment level

dt$Cd.level <- trt[match(dt[, plant.number], Plant), lvl]
dt$positon <- trt[match(dt[, plant.number], Plant), Tray]




dt[, code := paste(plant.number, positon, Cd.level, sep='')]
dt[, bioreps := paste(positon, Cd.level, sep='')]
dt[, .(plants = paste(plant.number, collapse = ", "),
       count = .N), 
   by = bioreps]

# AVERAGE BY PLANT CODE 
#Each plant should have one observation, per data policy we will keep the 
#average between repeat measurements

dt[, .N, .(code)]

dt[, 9:(ncol(dt)-2) := lapply(.SD, as.numeric), .SDcols = 9:(ncol(dt)-2)]

dt[, photo_m := mean(Photo), .(code)]
dt[, cond_m := mean(Cond), .(code)]
dt[, ci_m := mean(Ci), .(code)]
library(data.table)

# Split the time into hours, minutes, and seconds
dt[, c("hours", "minutes", "seconds") := tstrsplit(HHMMSS, ":", type.convert = TRUE)]

# Convert hours and minutes to seconds and add them up
dt[, secondsOD := hours * 3600 + minutes * 60 + seconds]

dt[secondsOD<20000, 'secondsOD'] <- dt[ secondsOD<20000, 'secondsOD']+(12*3600)

dt <- unique(dt[, .(code, secondsOD, positon, plant.number, Cd.level, photo_m, cond_m, ci_m)]) # 48 observations


ggplot(data=data.frame(spp=c('Pelo Crespito', 'Lechuga', 'Cuerno de venado', "Pelo de caballo"), area=c(50,75,15,5)), aes(factor(spp), area, fill=spp))+
  geom_bar(stat='identity')+
  ylab(label = 'Area m2')+
  xlab(label = 'Especies')+
  theme(axis.text.x = element_text(angle = 45, hjust = 1))
