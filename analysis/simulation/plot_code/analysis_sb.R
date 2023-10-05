library(ggplot2)
library(reshape2)
library(ggpubr)
library(rstatix)
library(scales)
setwd('/data/gpfs/projects/punim1662/yidi_projects/stablemate')

############### Load results ###############
load('./simulation/r_object/sb_50.RData')
AUC_record_50 <- AUC_record
time_record_50 <- time_record
sensSB_record_50 <- sensSB_record 
specSB_record_50 <- specSB_record
sensNSB_record_50 <- sensNSB_record
specNSB_record_50 <- specNSB_record
sumSB_record_50 <- (sensSB_record_50 + specSB_record_50)/2 # Balanced accuracy
sumNSB_record_50 <- (sensNSB_record_50 + specNSB_record_50)/2
F1SB_record_50 <- F1SB_record ; F1SB_record_50[is.na(F1SB_record_50)] <- 0
F1NSB_record_50 <- F1NSB_record ; F1NSB_record_50[is.na(F1NSB_record_50)] <- 0
MSE_record_50 <- -log10(-MSE_record)

load('./simulation/r_object/sb_100.RData')
AUC_record_100 <- AUC_record
time_record_100 <- time_record
sensSB_record_100 <- sensSB_record 
specSB_record_100 <- specSB_record
sensNSB_record_100 <- sensNSB_record
specNSB_record_100 <- specNSB_record
sumSB_record_100 <- (sensSB_record_100 + specSB_record_100)/2
sumNSB_record_100 <- (sensNSB_record_100 + specNSB_record_100)/2
F1SB_record_100 <- F1SB_record ; F1SB_record_100[is.na(F1SB_record_100)] <- 0
F1NSB_record_100 <- F1NSB_record ; F1NSB_record_100[is.na(F1NSB_record_100)] <- 0
MSE_record_100 <- -log10(-MSE_record)

load('./simulation/r_object/sb_200.RData')
AUC_record_200<- AUC_record
time_record_200 <- time_record
sensSB_record_200 <- sensSB_record 
specSB_record_200 <- specSB_record
sensNSB_record_200 <- sensNSB_record
specNSB_record_200 <- specNSB_record
sumSB_record_200 <- (sensSB_record_200 + specSB_record_200)/2
sumNSB_record_200 <-( sensNSB_record_200 + specNSB_record_200)/2
F1SB_record_200 <- F1SB_record ; F1SB_record_200[is.na(F1SB_record_200)] <- 0
F1NSB_record_200 <- F1NSB_record ; F1NSB_record_200[is.na(F1NSB_record_200)] <- 0
MSE_record_200 <- -log10(-MSE_record)



########################### Combine and collapse ###############################
#Combine result of different simulation scenarios 
AUC_record <- rbind(AUC_record_50, AUC_record_100, AUC_record_200)
time_record <- rbind(time_record_50, time_record_100, time_record_200)
sensSB_record <- rbind(sensSB_record_50, sensSB_record_50, sensSB_record_50) 
specSB_record <- rbind(specSB_record_50, specSB_record_50, specSB_record_50)
sensNSB_record <- rbind(sensNSB_record_50, sensNSB_record_100, sensNSB_record_200)
specNSB_record <- rbind(sensNSB_record_50, sensNSB_record_100, sensNSB_record_200)
sumSB_record <-  rbind(sumSB_record_50, sumSB_record_100, sumSB_record_200)
sumNSB_record <-  rbind(sumNSB_record_50, sumNSB_record_100, sumNSB_record_200)
F1SB_record <- rbind(F1SB_record_50, F1SB_record_100, F1SB_record_200)
F1NSB_record <- rbind(F1NSB_record_50, F1NSB_record_100, F1NSB_record_200)
MSE_record <- rbind(MSE_record_50, MSE_record_100, MSE_record_200)

#Create simulation id and sizes of systems simulated
id <- 1:nrow(AUC_record)
size <- c(rep(50, nrow(AUC_record_50)), rep(100, nrow(AUC_record_100)), rep(200, nrow(AUC_record_200)))

#Collapse the data frame for ggplot
AUC_record <- melt(data.frame(AUC_record, id= id),
                   id.vars = 'id', variable.name =  'method', value.name = 'AUC')
sensSB_record <- melt(data.frame(sensSB_record, id= id),
                      id.vars = 'id', variable.name =  'method', value.name = 'SensitivitySB')
specSB_record <- melt(data.frame(specSB_record, id= id),
                      id.vars = 'id', variable.name =  'method', value.name = 'SpecificitySB')
sensNSB_record <- melt(data.frame(sensNSB_record, id= id),
                       id.vars = 'id', variable.name =  'method', value.name = 'SensitivityNSB')
specNSB_record <- melt(data.frame(specNSB_record, id= id),
                       id.vars = 'id', variable.name =  'method', value.name = 'SpecificityNSB')
sumSB_record <- melt(data.frame(sumSB_record, id= id),
                     id.vars = 'id', variable.name =  'method', value.name = 'SumSB')
sumNSB_record <- melt(data.frame(sumNSB_record, id= id),
                      id.vars = 'id', variable.name =  'method', value.name = 'SumNSB')
F1SB_record <- melt(data.frame(F1SB_record, id= id),
                    id.vars = 'id', variable.name =  'method', value.name = 'F1SB')
F1NSB_record <- melt(data.frame(F1NSB_record, id= id),
                     id.vars = 'id', variable.name =  'method', value.name = 'F1NSB')
MSE_record <- melt(data.frame(MSE_record, id= id),
                   id.vars = 'id', variable.name =  'method', value.name = 'NMSE')
time_record <- melt(data.frame(time_record, id= id),
                    id.vars = 'id', variable.name =  'method', value.name = 'Seconds')

#Add factor level to x variable (method)
AUC_record$method <- factor(AUC_record$method, levels = c('st2estar', 'SR'))
sensSB_record$method <- factor(sensSB_record$method, levels = c('st2estar', 'SR'))
specSB_record$method <- factor(specSB_record$method, levels = c('st2estar', 'SR'))
sensNSB_record$method <- factor(sensNSB_record$method, levels = c('st2estar', 'SR'))
specNSB_record$method <- factor(specNSB_record$method, levels = c('st2estar', 'SR'))
sumSB_record$method <- factor(sumSB_record$method, levels = c('st2estar', 'SR'))
sumNSB_record$method <- factor(sumNSB_record$method, levels = c('st2estar', 'SR'))
F1SB_record$method <- factor(F1SB_record$method, levels = c('st2estar', 'SR'))
F1NSB_record$method <- factor(F1NSB_record$method, levels = c('st2estar', 'SR'))
MSE_record$method <- factor(MSE_record$method, levels = c('st2estar', 'ols', 'SR', 'lasso', 'rf'))
time_record$method <- factor(time_record$method, levels = c('st2estar', 'SR'))

#Add sizes of systems simulated as a metadata
AUC_record$size <- size[AUC_record$id]
sensSB_record$size <- size[sensSB_record$id]
specSB_record$size <- size[specSB_record$id]
sensNSB_record$size <- size[sensNSB_record$id]
specNSB_record$size <- size[specNSB_record$id]
sumSB_record$size <- size[sumSB_record$id]
sumNSB_record$size <- size[sumNSB_record$id]
F1SB_record$size <- size[F1SB_record$id]
F1NSB_record$size <- size[F1NSB_record$id]
MSE_record$size <- size[MSE_record$id]
time_record$size <- size[time_record$id]


############################### Stat test ######################################
#Paired t-tests
AUC_test <- AUC_record %>%
  group_by(size) %>%  t_test(AUC ~ method, 
                             ref.group = 'st2estar', paired = TRUE, alternative = 'greater')
sensSB_test <- sensSB_record %>%
  group_by(size) %>%  t_test(SensitivitySB ~ method, 
                             ref.group = 'st2estar', paired = TRUE, alternative = 'greater')
specSB_test <- specSB_record %>%
  group_by(size) %>%  t_test(SpecificitySB ~ method, 
                             ref.group = 'st2estar', paired = TRUE, alternative = 'greater')
sensNSB_test <- sensNSB_record %>%
  group_by(size) %>%  t_test(SensitivityNSB ~ method, 
                             ref.group = 'st2estar', paired = TRUE, alternative = 'greater')
specNSB_test <- specNSB_record %>%
  group_by(size) %>%  t_test(SpecificityNSB ~ method, 
                             ref.group = 'st2estar', paired = TRUE, alternative = 'greater')
sumSB_test <- sumSB_record %>%
  group_by(size) %>%  t_test(SumSB ~ method,
                             ref.group = 'st2estar', paired = TRUE, alternative = 'greater')
sumNSB_test <- sumNSB_record %>%
  group_by(size) %>%  t_test(SumNSB ~ method,
                             ref.group = 'st2estar', paired = TRUE, alternative = 'greater')
F1SB_test <- F1SB_record %>%
  group_by(size) %>%  t_test(F1SB ~ method,
                             ref.group = 'st2estar', paired = TRUE, alternative = 'greater')
F1NSB_test <- F1NSB_record %>%
  group_by(size) %>%  t_test(F1NSB ~ method,
                             ref.group = 'st2estar', paired = TRUE, alternative = 'greater')

MSE_test <- MSE_record %>%
  group_by(size) %>%  t_test(NMSE ~ method,
                             ref.group = 'st2estar', paired = TRUE, alternative = 'greater')

time_test <- time_record %>%
  group_by(size) %>%  t_test(Seconds ~ method,
                             ref.group = 'st2estar', paired = TRUE)

#Prepare the p-value var for plot
AUC_test <- AUC_test %>% add_xy_position(x = "method") %>% add_significance()
AUC_test$y.position <- c(1.05)
sensSB_test <- sensSB_test %>% add_xy_position(x = "method")%>% add_significance()
specSB_test <- specSB_test %>% add_xy_position(x = "method")%>% add_significance()
sensNSB_test <- sensNSB_test %>% add_xy_position(x = "method")%>% add_significance()
specNSB_test <- specNSB_test %>% add_xy_position(x = "method")%>% add_significance()
sumSB_test <- sumSB_test %>% add_xy_position(x = "method")%>% add_significance()
sumSB_test$y.position <- c(1.05)
sumNSB_test <- sumNSB_test %>% add_xy_position(x = "method")%>% add_significance()
sumNSB_test$y.position <- c(1.05)
F1SB_test <- F1SB_test %>% add_xy_position(x = "method") %>% add_significance()
F1SB_test$y.position <- c(1.05)
F1NSB_test <- F1NSB_test %>% add_xy_position(x = "method")%>% add_significance()
F1SB_test$y.position <- c(1.05)
MSE_test <- MSE_test %>% add_xy_position(x = "method")%>% add_significance()
time_test <- time_test %>% add_xy_position(x = "method")%>% add_significance()


################################### Plot #######################################
# Scaling x-axis a
reverselog_trans <- function(base = exp(1)) {
  trans <- function(x) -log(x, base)
  inv <- function(x) base^(-x)
  trans_new(paste0("reverselog-", format(base)), trans, inv, 
            log_breaks(base = base), 
            domain = c(1e-100, Inf))
}

# Figure S2A
ggboxplot(
  AUC_record, x = "method", y = "AUC", fill = "method",
  facet = 'size', add = "jitter", color = 'black',  ylim = c(0,1.1)
) + stat_pvalue_manual(AUC_test) + theme(text = element_text(size = 15)) +
  theme_gray()+
  theme(text = element_text(size = 15), legend.position="none")+
  scale_fill_manual(values = c('white', 'lavenderblush4')) +
  scale_x_discrete(labels = c('StableMate',  'SR'))  +
  scale_y_continuous(breaks = c(0, 0.2,0.4,0.6,0.8,1))  + xlab('Methods')+
  ylab('AUC (SB)')

# Figure S2B
ggboxplot(
  sumSB_record, x = "method", y = "SumSB", fill = "method",
  facet = 'size', add = "jitter", color = 'black',  ylim = c(0,1.1)
) + stat_pvalue_manual(sumSB_test) + theme(text = element_text(size = 15)) +
  theme_gray()+
  theme(text = element_text(size = 15), legend.position="none")+
  scale_fill_manual(values = c('white', 'lavenderblush4')) +
  scale_x_discrete(labels = c('StableMate',  'SR'))  +
  scale_y_continuous(breaks = c(0, 0.2,0.4,0.6,0.8,1))  + xlab('Methods') +
  ylab('Balanced accuracy (SB)')

# Figure S2C
ggboxplot(
  sumNSB_record, x = "method", y = "SumNSB", fill = "method",
  facet = 'size', add = "jitter", color = 'black',  ylim = c(0,1.1)
) + stat_pvalue_manual(sumNSB_test) + theme(text = element_text(size = 15)) +
  theme_gray()+
  theme(text = element_text(size = 15), legend.position="none")+
  scale_fill_manual(values = c('white', 'lavenderblush4')) +
  scale_x_discrete(labels = c('StableMate',  'SR'))  +
  scale_y_continuous(breaks = c(0, 0.2,0.4,0.6,0.8,1))  + xlab('Methods') +
  ylab('Balanced accuracy (NSB)')

# Figure S2D
ggboxplot(
  MSE_record, x = "method", y = "NMSE", fill = "method", add = "jitter", color = 'black',  
) + stat_pvalue_manual(MSE_test) + theme(text = element_text(size = 15)) +
  theme_gray()+
  theme(text = element_text(size = 15), legend.position="none")+
  scale_fill_manual(values = c('white', '#388080', 'lavenderblush4', 'blue' , 'red')) +
  scale_x_discrete(labels = c('StableMate', 'OLS',  'SR', 'LASSO', 'RF')) + xlab('')  + ylab('Negative log10 MSE') +
  facet_wrap(~size, ncol = 1) + xlab('Methods')

# Figure S2E
ggboxplot(
  time_record, x = "method", y = "Seconds", fill = "method", add = "jitter", color = 'black',
) + stat_pvalue_manual(time_test) +
  theme_gray() + 
  theme(text = element_text(size = 15), legend.position="none")+
  scale_fill_manual(values = c('white', 'lavenderblush4')) +
  scale_x_discrete(labels = c('StableMate',  'SR'))   + xlab('Methods')+
  facet_wrap(~size, ncol = 1)


############### Plot that are not included in the manuscript ############### 
# Included tests for F1 scores
ggboxplot(
  F1SB_record, x = "method", y = "F1SB", fill = "method",
  facet = 'size', add = "jitter", color = 'lightgrey',  ylim = c(0,1.1)
) + stat_pvalue_manual(F1SB_test) + theme(text = element_text(size = 15)) +
  theme(text = element_text(size = 15), legend.position="none")+
  scale_fill_manual(values = c('white', 'lavenderblush4')) +
  scale_x_discrete(labels = c('StableMate',  'SR'))  +
  scale_y_continuous(breaks = c(0, 0.2,0.4,0.6,0.8,1))  + xlab('') + ylab('F1: SB')

ggboxplot(
  F1NSB_record, x = "method", y = "F1NSB", fill = "method",
  facet = 'size', add = "jitter", color = 'lightgrey',  ylim = c(0,1.1)
) + stat_pvalue_manual(F1NSB_test) + theme(text = element_text(size = 15)) +
  theme(text = element_text(size = 15), legend.position="none")+
  scale_fill_manual(values = c('#D81B60', '#FFC107')) +
  scale_x_discrete(labels = c('StableMate',  'SR'))  +
  scale_y_continuous(breaks = c(0, 0.2,0.4,0.6,0.8,1))  + xlab('')+ ylab('F1: NSB')
