# Load stablemate functions
# Some essential package has already been loaded
setwd('/data/gpfs/projects/punim1662/yidi_projects/stablemate')
load("./intro/r_object/introdata.RData")
source('./stablemate.R')
require(ggplot2)
require(ggpubr)

############### Figure 1C ###############
plot(mod, label = c(bquote("X"[1]), bquote("X"[2]),  bquote("X"[3]),
                    bquote("X"[4]), bquote("X"[5]), bquote("X"[6]),
                    bquote("X"[7]), bquote("X"[8]), bquote("X"[9]),
                    bquote("X"[10]), bquote("X"[11]), bquote("X"[12]),
                    bquote("X"[13]), bquote("X"[14]), bquote("X"[15]),
                    bquote("X"[16]), bquote("X"[17]), bquote("X"[18]), bquote("X"[19])), label_size = 5, parse = T)+
  theme(text = element_text(size = 15)) + 
  xlab('Predictivity score') + ylab('Stability score')


############### Figure 1B ###############
p1 <- ggplot() + geom_point(aes(trainX[,'X13'],trainY, col = group)) + 
  geom_smooth(aes(trainX[,'X13'],trainY, col = group), method = 'lm') +
  scale_color_manual('Environments',values = c('black','#388080','orange'),
                     labels = c(bquote("e"[1]), bquote("e"[2]), bquote("e"[3]))) +
  xlab(bquote(paste('X'[3],' (stable predictor)'))) + ylab('Y') + 
  theme(text = element_text(size = 15)) + guides(color = 'none')

p2 <- ggplot() + geom_point(aes(trainX[,'X26'],trainY, col = group)) + 
  geom_smooth(aes(trainX[,'X26'],trainY, col = group), method = 'lm') +
  scale_color_manual('Environments',values = c('black','#388080','orange'),
                     labels = c(bquote("e"[1]), bquote("e"[2]), bquote("e"[3]))) +
  xlab(bquote(paste('X'[15],' (unstable predictor)'))) + ylab('Y') + 
  theme(text = element_text(size = 15), legend.position = 'top') + ylab('') 

ggarrange(p1,p2, common.legend = T)
