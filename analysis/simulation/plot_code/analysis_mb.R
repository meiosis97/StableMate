library(ggplot2)
library(reshape2)
library(ggpubr)
library(rstatix)
setwd('/data/gpfs/projects/punim1662/yidi_projects/stablemate')

############### Load results ###############
load('./simulation/r_object/mb_50.RData')
BIC_record_50 <- BIC_record
AUC_record_50 <- AUC_record
time_record_50 <- time_record

load('./simulation/r_object/mb_100.RData')
BIC_record_100 <- BIC_record
AUC_record_100 <- AUC_record
time_record_100 <- time_record

load('./simulation/r_object/mb_200.RData')
BIC_record_200 <- BIC_record
AUC_record_200 <- AUC_record
time_record_200 <- time_record


############### Combine and collapse ############### 
#Combine result of different simulation scenarios 
BIC_record <- rbind(BIC_record_50, BIC_record_100, BIC_record_200)
AUC_record <- rbind(AUC_record_50, AUC_record_100, AUC_record_200)
time_record <- rbind(time_record_50, time_record_100, time_record_200)

#Create simulation id and sizes of systems simulated
id <- 1:nrow(BIC_record)
size <- c(rep(50, nrow(BIC_record_50)), rep(100, nrow(BIC_record_100)), rep(200, nrow(BIC_record_200)))

#Collapse the data frame for ggplot
BIC_record <- melt(data.frame(BIC_record, id= id),
                   id.vars = 'id', variable.name =  'method', value.name = 'NBIC')
AUC_record <- melt(data.frame(AUC_record, id= id),
                   id.vars = 'id', variable.name =  'method', value.name = 'AUC')
time_record <- melt(data.frame(time_record, id= id),
                    id.vars = 'id', variable.name =  'method', value.name = 'Seconds')

#Add factor level to x variable (method)
BIC_record$method <- factor(BIC_record$method, levels = c('st2estar', 'st2e'))
AUC_record$method <- factor(AUC_record$method, levels = c('st2estar', 'st2e', 'SR'))
time_record$method <- factor(time_record$method, levels = c('st2estar', 'st2e', 'SR'))

#Add sizes of systems simulated as a metadata
BIC_record$size <- size[BIC_record$id]
AUC_record$size <- size[AUC_record$id]
time_record$size <- size[time_record$id]


############### Stat test############### 
#Paired t-tests
BIC_test <- BIC_record %>%
  group_by(size) %>%  t_test(NBIC ~ method, 
                             ref.group = 'st2estar', paired = TRUE, alternative = 'greater')
AUC_test <- AUC_record %>%
  group_by(size) %>%  t_test(AUC ~ method,
                             ref.group = 'st2estar', paired = TRUE, alternative = 'greater')
time_test <- time_record %>%
  group_by(size) %>%  t_test(Seconds ~ method,
                             ref.group = 'st2estar', paired = TRUE)

#Prepare the p-value var for plot
BIC_test <- BIC_test %>% add_xy_position(x = "method")  %>% add_significance()
AUC_test <- AUC_test %>% add_xy_position(x = "method")
AUC_test$y.position <- c(1.05,1.1,1.05,1.1,1.05,1.1)
time_test <- time_test %>% add_xy_position(x = "method")


############### Plot ############### 
# Figure S1A
ggboxplot(
  AUC_record, x = "method", y = "AUC", fill = "method",
  facet = 'size', add = "jitter", color = 'black',  ylim = c(0,1.1)
) + stat_pvalue_manual(AUC_test) + theme_gray() +
  theme(text = element_text(size = 15), legend.position="none")+
  scale_fill_manual(values = c('white', '#388080', 'lavenderblush4')) +
  scale_x_discrete(labels = c('ST2*', 'ST2', 'SR'))  +
  scale_y_continuous(breaks = c(0, 0.2,0.4,0.6,0.8,1))  + xlab('')

# Figure S1B
ggboxplot(
  BIC_record, x = "method", y = "NBIC", fill = "method",
  facet = 'size', add = "jitter", color = 'black', 
) + stat_pvalue_manual(BIC_test)  + theme_gray() + 
  theme(text = element_text(size = 15), legend.position="none") +
  scale_fill_manual(values = c('white', '#388080'))+
  scale_x_discrete(labels = c('ST2*', 'ST2'))  + xlab('')

# Figure S1C
ggboxplot(
  time_record, x = "method", y = "Seconds", fill = "method",
  facet = 'size', add = "jitter", color = 'black', 
) + stat_pvalue_manual(time_test) + theme(text = element_text(size = 15)) +
  theme_gray() + 
  theme(text = element_text(size = 15), legend.position="none")+
  scale_fill_manual(values = c('white', '#388080', 'lavenderblush4')) +
  scale_x_discrete(labels = c('ST2*', 'ST2', 'SR')) + xlab('Methods')

