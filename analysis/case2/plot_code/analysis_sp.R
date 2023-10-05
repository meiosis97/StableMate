############### The following preloading was used for the entire script (MUST RUN) ############### 
############### All the plot sessions (except for Figure S4C) can be ran independently of each others after running the preloading ############### 
require(dplyr)
require(vegan)
require(ggpubr)
require(pheatmap)
require(RColorBrewer)
require(pROC)
require(rstatix)
require(ggvenn)
setwd('F:/stablemate/stablemate')
source('./stablemate.R')
data_list <- readRDS('./case2/data/processed_data_sp.RDS')

# Convert study names to country names
meta_data <- data_list[[3]]
cohort <- meta_data$cohort
country_table <- c('Austria', 'Italy A', 'Italy B',
                  'USA', 'Germany', 'Japan', 'China', 'France')
names(country_table) <- unique(cohort)
country <- country_table[cohort]

# Set color scheme
cols <- c('black', 'blue', 'red', 'palevioletred1', 'darkgreen', 'lavenderblush4', 'purple3','orange')
names(cols) <- country_table

# Make sample type (Normal or CRC) variable 
sample_type <- meta_data$sample_type
sample_type[sample_type == 'control'] <- 'Normal'

# Create a seperate variable for storing species names
species <- colnames(data_list[[1]])

############### Figure S4A ############### 
load('./case2/r_object/sp/all/all.RData')
# Make plot
plot.SRST2E(mod, label_size = 4) +  theme(text = element_text(size = 15), legend.position = 'bottom') + ylim(0,1)


############### Figure 3A ############### 
load('./case2/r_object/sp/all/all.RData')

# MDS and permutation ANOVA on the full data
mds_all <- X %>% dist %>% cmdscale()
dist_mat_all <- mds_all %>% dist
permanova_all <- adonis2(dist_mat_all ~ country + sample_type)

# Prepare the data for plotting 
mds_all <- data.frame(mds_all) %>%
  mutate(`Sample type` = sample_type) %>% mutate(Cohort = country) 

# Make plot 3A1
p1 <- ggscatter(mds_all, x = "X1", y = "X2", 
                color = 'Sample type',
                palette =  c('orange','black'),
                size = 1.8, 
                ellipse = TRUE,
                ellipse.type = "convex",ellipse.alpha = 0.02,
                repel = TRUE)+ theme_gray() + 
  annotate("text", x=-Inf, y = Inf, label = paste("R2 =", round(permanova_all$R2[2],6)),parse = F, size = 5,
           vjust=1, hjust=0)+
  theme(text = element_text(size=15)) + xlab('') + ylab('PCoA2 on all the 313 species') +
  guides(color = 'none', fill = 'none')

p2 <- ggscatter(mds_all, x = "X1", y = "X2", 
                color = 'Cohort',
                palette =  cols,
                size = 1.8, 
                ellipse = TRUE,
                ellipse.type = "convex", ellipse.alpha = 0.02, 
                repel = TRUE) + theme_gray() + 
  annotate("text", x = -Inf, y = Inf, label = paste("R2 =", round(permanova_all$R2[1],6)),parse = F, size = 5,
           vjust=1, hjust=0)+
  theme(text = element_text(size=15),
        axis.title.x = element_text(hjust = -3)) +
  xlab('PCoA1 on all the 313 species') + ylab('') 

# MDS and permutation ANOVA on the 23 species selected by stablemate
var_sel <- print(mod)$SB_selected
mds_sm <- X[,var_sel] %>% dist %>% cmdscale()
dist_mat_sm <- mds_sm %>% dist
permanova_sm <- adonis2(dist_mat_sm ~  country + sample_type)

# Prepare the data for plotting 
mds_sm <- data.frame(mds_sm) %>%
  mutate(`Sample type` = sample_type) %>% mutate(Cohort = country) 

# Make plot 3A2
p3 <- ggscatter(mds_sm, x = "X1", y = "X2", 
                color = 'Sample type',
                palette =  c('orange','black'),
                size = 1.8, 
                ellipse = TRUE,
                ellipse.type = "convex",ellipse.alpha = 0.02,
                repel = TRUE)+ theme_gray() + 
  annotate("text", x = -Inf, y = Inf, label = "R2 = 0.04773",parse = F, size = 5,
           vjust=1, hjust=0)+
  theme(text = element_text(size=15)) + xlab('') + ylab('PCoA2 on the 23 species selected')  +
  guides(color = 'none', fill = 'none')

p4 <- ggscatter(mds_sm, x = "X1", y = "X2", 
                color = 'Cohort',
                palette =  cols,
                size = 1.8, 
                ellipse = TRUE,
                ellipse.type = "convex", ellipse.alpha = 0.02,
                repel = TRUE) + theme_gray() + 
  annotate("text", x = -Inf, y = Inf, 
           label = "R2 = 0.08861",parse = F, size = 5,
           vjust=1, hjust=0)+ylab('')+
  theme(text = element_text(size=15),
        axis.title.x = element_text(hjust = -20)) +
  xlab('PCoA1 on the 23 species selected') + ylab('') 

# Combine plots
ggarrange(p1,p2,p3,p4, align = 'v', common.legend = T, legend = 'right')

# facet warp plot 
mds_all$nvar <- 'All 313 species'
mds_sm$nvar <- 'StableMate 23 species'
mds_warp <- rbind(mds_all,mds_sm)
dat_text <- data.frame(
  label = c(paste("R2 =", round(permanova_all$R2[2],6)), paste("R2 =", round(permanova_sm$R2[2],6))),
  nvar   = c('All 313 species', 'StableMate 23 species')
)
ggscatter(mds_warp, x = "X1", y = "X2", 
                color = 'Sample type',
                palette =  c('orange','black'),
                size = 1.8, 
                ellipse = TRUE,
                ellipse.type = "convex",ellipse.alpha = 0.02,
                repel = TRUE)+ theme_gray() + facet_wrap(~nvar, ncol = 2, scales = 'free')+
  theme(text = element_text(size=15), legend.position="bottom") + xlab('PCoA1') + ylab('PCoA2')+
  geom_text(data = dat_text, aes(x = -Inf, y = Inf, label = label), hjust = 0, vjust = 1, size = 5)


dat_text <- data.frame(
  label = c(paste("R2 =", round(permanova_all$R2[1],6)), paste("R2 =", round(permanova_sm$R2[1],6))),
  nvar   = c('All 313 species', 'StableMate 23 species')
)
ggscatter(mds_warp, x = "X1", y = "X2", 
          color = 'Cohort',
          palette =  cols,
          size = 1.8, 
          ellipse = TRUE,
          ellipse.type = "convex",ellipse.alpha = 0.02,
          repel = TRUE)+ theme_gray() + facet_wrap(~nvar, ncol = 2, scales = 'free')+
  theme(text = element_text(size=15), legend.position="bottom") + xlab('PCoA1') + ylab('PCoA2') +
  geom_text(data = dat_text, aes(x = -Inf, y = Inf, label = label), hjust = 0, vjust = 1, size = 5)

# Variance parition 
var_source <- c('Cohort','Sample type', 'Residual' )
var_source <- factor(var_source, levels =  c('Residual','Sample type', 'Cohort' ))
col_plot <- data.frame(y = c(rep('All 313 species', 3), rep('StableMate 23 species', 3)),
                       var_source = var_source,
                       R2 = c(permanova_all$R2[-4], permanova_sm$R2[-4]))
ggplot(col_plot) + geom_col(aes(y = y, x = R2,
                        fill =var_source), position = position_stack()) + 
  scale_fill_manual('Variable', values = c('grey','blue','red')) + theme(text = element_text(size = 15)) +
  ylab('') 


############### Figure 3B ############### 
# Create a matrix for storing predictivity scores
score_mat <- data.frame(matrix(NA, ncol = length(species), nrow = length(country_table)))
# Create a indicator matrix for indicating predictive and stable predictors for each cohort.
indi_stabpred <- data.frame(matrix(0, ncol = length(species), nrow = length(country_table)))
# Create a indicator matrix for indicating predictive predictors for each cohort.
indi_pred <- data.frame(matrix(0, ncol = length(species), nrow = length(country_table)))
# Indexing the matrices
colnames(score_mat) <- colnames(indi_stabpred)  <- colnames(indi_pred) <- species
rownames(score_mat) <- rownames(indi_stabpred)  <- rownames(indi_pred)  <- country_table

for(c in names(country_table)){
  
  # Load stablemate results for each cohort
  dir <- paste('./case2/r_object/sp/one_by_one/', c, '.RData', sep = '')
  load(dir)
  
  # Record predictivity scores
  score_mat[country_table[c],] <- mod$prediction_ensemble$selection_prob$marginal$selection_prob[-1]
  
  # Record MB and SB selections
  MB_sel <- print(mod)$MB_selected
  SB_sel <- print(mod)$SB_selected
  indi_pred[country_table[c],MB_sel] <- 1
  indi_stabpred[country_table[c],SB_sel] <- 1
  
}

# Load stablemate results obtained on the pooled data
load('./case2/r_object/sp/all/all.RData')

# Record predictivity scores
score_mat <- rbind(score_mat, mod$prediction_ensemble$selection_prob$marginal$selection_prob[-1])
rownames(score_mat)[9] <- 'Pooled'

# Record MB selections
MB_sel <- print(mod)$MB_selected
indi_pred <- rbind(indi_pred, 0)
indi_pred[9,MB_sel] <- 1
rownames(indi_pred)[9] <- 'Pooled'

# Record SB selections
SB_sel <- print(mod)$SB_selected
indi_stabpred <- rbind(indi_stabpred, 0)
indi_stabpred[9,SB_sel] <- 1
rownames(indi_stabpred)[9] <- 'Pooled'

# Rank predictivity scores
rank_mat <- data.frame(apply(-score_mat,1,rank, ties.method = 'max'), check.names = F)
rank_mat <- rank_mat[,c('Pooled',country_table)]
rank_mat$Mean <- round(rowMeans(rank_mat[,-1]),2)

# Keep only the species that are selected as predictive by at least one cohort
rank_mat <- rank_mat[colSums(indi_pred) > 0,]
indi_stabpred <- indi_stabpred[,colSums(indi_pred) > 0]

# Order species by their mean ranking in predictivity scores
rank_mat <- rank_mat[order(rank_mat$Pooled),]

# Create labels for the heatmap
label_mat <- rank_mat
# Label a cell by asterisk if the corresponding species of the cell is selected as predictive and stable for the corresponding cohort
for(i in colnames(indi_stabpred)){
  
  for(j in rownames(indi_stabpred)){
    
    label_mat[i,j] <- paste(label_mat[i,j], ifelse(indi_stabpred[j,i]==1, '*',''))
    
  }
  
}

# Create row annotation
anno_row <- data.frame(frequency = colSums(indi_stabpred[,rownames(rank_mat)]))
rownames(anno_row) <- rownames(rank_mat)

# Make the heatmap
pheatmap(rank_mat[1:50,], cluster_rows = F, cluster_cols = F, display_numbers = label_mat[1:50,], color =  
           c(colorRampPalette(brewer.pal(n = 4, name = "Blues")[1:2])(nrow(rank_mat)),
             colorRampPalette(brewer.pal(n = 4, name = "Blues")[3:4])(max(rank_mat)-nrow(rank_mat))),angle_col = 45,
         annotation_row = anno_row[1:50,,drop = F], fontsize = 11)

# Make the annotation
ggplot() + geom_col(aes(1:length(species),1:length(species),color = 1:length(species))) +
  geom_col(aes(1:5,1:5,fill = as.factor(1:5))) +
  theme(text = element_text(size = 15),legend.position = 'top')+
  scale_color_gradientn('Rank',colors = c(colorRampPalette(brewer.pal(n = 4, name = "Blues")[1:2])(nrow(rank_mat)),
                                          colorRampPalette(brewer.pal(n = 4, name = "Blues")[3:4])(max(rank_mat)-nrow(rank_mat))))+
  scale_fill_manual('Selection frequency',values = RColorBrewer::brewer.pal(5,'Greens'))


############### Figure S4B ############### 
load('./case2/r_object/sp/one_by_one/FengQ_2015.RData')
# Make plot
plot.SRST2E(mod) +  theme(text = element_text(size = 15)) + ylim(0,1)


############### Figure 3C ############### 
# Create a matrix for storing LODO results
AUC_mat <- data.frame(matrix(NA, ncol = 5, nrow = length(country_table)))
colnames(AUC_mat) <- c('SM-stab','SM-pred','GLM', 'LASSO', 'RF')
rownames(AUC_mat) <- country_table

# Create a function for LODO benchmark
benchmark <- function(){
  
  # Weight samples by cohort-case (disease status) sizes
  n_sample <- length(Ytrain)
  agg_group <- paste(Ytrain, envs)
  w <- 1/(table(agg_group)[agg_group])
  w <- w/sum(w)*n_sample
  
  # Prediction with stable ensemble
  Yhat <- predict(mod, Xtest,  predict_fun = predict_logit)
  auc_smstab <- suppressMessages(auc(as.factor(Ytest), Yhat))
  
  # Prediction with predictive ensemble
  Yhat <- predict(mod, Xtest,  assay = 'prediction_ensemble', predict_fun = predict_logit)
  auc_smpred <- suppressMessages(auc(as.factor(Ytest), Yhat))

  # Train GLM (logistic regression) and test on the testing cohort
  glm_mod <- glm(Ytrain~., data = data.frame(Xtrain), family = 'binomial', weights = w)
  Yhat <- predict(glm_mod, data.frame(Xtest), type = 'response')
  auc_glm <- suppressMessages(auc(as.factor(Ytest), as.numeric(Yhat)))

  # Train Lasso (logistic regression) and test on the testing cohort
  # Training environments were used as folds for cross validation
  glmnet_mod <- cv.glmnet(Xtrain, Ytrain, family = 'binomial',
                         weights = w, foldid = as.numeric(as.factor(envs))) 
  Yhat <- predict(glmnet_mod, Xtest, s = glmnet_mod$lambda.1se)
  auc_lasso <- suppressMessages(auc(as.factor(Ytest), as.numeric(Yhat)))

  # Train RF  and test on the testing cohort
  rf_mod <- randomForest::randomForest(Xtrain, as.factor(Ytrain), weights = w)
  Yhat <- predict(rf_mod, Xtest, type = 'prob')[,2]
  auc_rf <- suppressMessages(auc(as.factor(Ytest), as.numeric(Yhat)))
  
  c(auc_smstab, auc_smpred, auc_glm, auc_lasso, auc_rf)
  
}

# LODO
for(c in names(country_table)){
  
  # Load stablemate results for each cohort
  dir <- paste('./case2/r_object/sp/leave_one_out/no_', c, '.RData', sep = '')
  load(dir)
  
  # Record predictivity scores
  AUC_mat[country_table[c],] <- benchmark()

}

# Prepare the data for boxplot
plot_data <- AUC_mat
plot_data$dataset <- rownames(plot_data)
plot_data <- melt(plot_data)

# Paired t-test
AUC_test <- plot_data   %>%  t_test(value ~ variable, 
                                     ref.group = 'SM-stab', paired = TRUE)
AUC_test$y.position <- c(0.9,0.925,0.95,0.975)

# Make the plot
ggboxplot(
  plot_data, x = "variable", y = "value", color = 'black'
) + stat_pvalue_manual(AUC_test, label = 'p.adj') +
  theme_gray() + geom_point(aes(x= variable, y = value, col = dataset),size=3) +
  ylab('AUC') + xlab('Method')+
  scale_color_manual('Cohort', values = cols)+
  theme(text = element_text(size = 15), legend.position = 'bottom') +
  scale_x_discrete(labels=paste(colnames(AUC_mat), round(colMeans(AUC_mat),3), sep = '\n'))  + 
  guides(color=guide_legend(ncol=2)) + ylim(0.5,1)


############### Figure 3D ############### 
load('./case2/r_object/sp/all/all.RData')

# Weight samples by cohort-case (disease status) sizes
n_sample <- length(Y)
agg_group <- paste(Y, envs)
w <- 1/(table(agg_group)[agg_group])
w <- w/sum(w)*n_sample

# Selection made by stablemate
sm_sel <- print(mod)$SB_selected

# Train RF model and obtain top RF selections
rf_mod <- randomForest::randomForest(X, as.factor(Y), weights = w)
rf_sel <- names(rf_mod$importance[order(rf_mod$importance, decreasing = T),])[1:length(sm_sel)]

# Train Lasso model and obtain top Lasso selections
lasso_mod <- glmnet(X, Y, family = 'binomial', weights = as.numeric(w), nlambda = 1000)
sel.matrix <- (lasso_mod$beta != 0) # Indicator matrix that record lasso selection path
first.entrance <- apply(sel.matrix, 1, which.max) # Calculate at which lambda each predictor is first selected
first.entrance[which(apply(sel.matrix, 1, sum) == 0)] <- Inf # If not selected in the entire path, set the selection order to infinity
first.entrance <- sort(first.entrance) 
lasso_sel <- names(first.entrance[1:length(sm_sel)]) # Get the first several predictors entering the selection path

# Create the Venn diagram
ggvenn(
  list("SB selection" = sm_sel, "RF selection" = rf_sel, "Lasso selection" = lasso_sel), 
  fill_color = c("#0073C2FF", "#EFC000FF", "#CD534CFF"),
  stroke_size = 0.6, set_name_size = 4, text_size = 8, show_percentage = F
)

# Get different subsets of predictors depicted in the Venn diagram
all_sel <- Reduce(intersect, list(sm_sel,rf_sel,lasso_sel)) # Predictors selected by all the three methods
sm_and_lasso_sel <- intersect(sm_sel, lasso_sel) # Predictors selected by stablemate and Lasso
sm_or_lasso_sel <- unique(c(sm_sel, lasso_sel)) # Predictors selected by stablemate or Lasso


############### Figure S4C (Running Figure 3D is required) ############### 
# Benchmark variable selection by training RF models with different predictor subsets and run LODO assessment on the RF models. 
benchmark <- function(){
  # Weight samples by cohort-case (disease status) sizes
  n_sample <- length(Ytrain)
  agg_group <- paste(Ytrain, envs)
  w <- 1/(table(agg_group)[agg_group])
  w <- w/sum(w)*n_sample
  
  res <- c()
  # Train RF with predictors selected by all the three methods
  rfmod <- randomForest::randomForest(Xtrain[,all_sel], as.factor(Ytrain), weights = w)
  Yhat <- predict(rfmod, Xtest, type = 'prob')[,2]
  res[1] <- suppressMessages(auc(as.factor(Ytest), as.numeric(Yhat)))
  
  # Train RF with predictors selected by stablemate and Lasso
  rfmod <- randomForest::randomForest(Xtrain[,sm_and_lasso_sel], as.factor(Ytrain), weights = w)
  Yhat <- predict(rfmod, Xtest, type = 'prob')[,2]
  res[2] <- suppressMessages(auc(as.factor(Ytest), as.numeric(Yhat)))
  
  # Train RF with predictors selected by stablemate or Lasso
  rfmod <- randomForest::randomForest(Xtrain[,sm_or_lasso_sel], as.factor(Ytrain), weights = w)
  Yhat <- predict(rfmod, Xtest, type = 'prob')[,2]
  res[3] <- suppressMessages(auc(as.factor(Ytest), as.numeric(Yhat)))
  
  # Train RF with predictors selected by Lasso
  rfmod <- randomForest::randomForest(Xtrain[,lasso_sel], as.factor(Ytrain), weights = w)
  Yhat <- predict(rfmod, Xtest, type = 'prob')[,2]
  res[4] <- suppressMessages(auc(as.factor(Ytest), as.numeric(Yhat)))
  
  # Train RF with predictors selected by RF
  rfmod <- randomForest::randomForest(Xtrain[,rf_sel], as.factor(Ytrain), weights = w)
  Yhat <- predict(rfmod, Xtest, type = 'prob')[,2]
  res[5] <- suppressMessages(auc(as.factor(Ytest), as.numeric(Yhat)))
  
  # Train RF with predictors selected by 
  rfmod <- randomForest::randomForest(Xtrain[,sm_sel], as.factor(Ytrain), weights = w)
  Yhat <- predict(rfmod, Xtest, type = 'prob')[,2]
  res[6] <- suppressMessages(auc(as.factor(Ytest), as.numeric(Yhat)))
  
  res
}

# Create a matrix for storing LODO results
AUC_mat <- data.frame(matrix(nrow = length(country_table), ncol =6))
colnames(AUC_mat) <- c('RF & Lasso & SM', 'Lasso & SM', 'Lasso | SM', 'Lasso', 'RF','SM')
rownames(AUC_mat) <- country_table

# LODO
for(c in names(country_table)){
  
  # Load stablemate results for each cohort
  dir <- paste('./case2/r_object/sp/leave_one_out/no_', c, '.RData', sep = '')
  load(dir)
  
  # Record predictivity scores
  AUC_mat[country_table[c],] <- benchmark()
  
}

# Prepare the data for boxplot
plot_data <- AUC_mat
plot_data$dataset <- rownames(plot_data)
plot_data <- melt(plot_data)

# Paired t-test
AUC_test <- plot_data   %>%  t_test(value ~ variable, 
                                    ref.group = 'SM', paired = TRUE)
AUC_test$y.position <- c(0.9,0.925,0.95,0.975,1)

# Make the plot
ggboxplot(
  plot_data, x = "variable", y = "value", color = 'black'
) + stat_pvalue_manual(AUC_test, label = 'p.adj') +
  theme_gray() + geom_point(aes(x= variable, y = value, col = dataset),size=3) +
  ylab('AUC') + xlab('Selection')+
  scale_color_manual('Cohort', values = cols)+
  theme(text = element_text(size = 15), legend.position = 'bottom',
        axis.text.x = element_text(angle = 45,hjust = 1)) +
  scale_x_discrete(labels=paste(colnames(AUC_mat), round(colMeans(AUC_mat),3), sep = '\n'))  + 
  guides(color=guide_legend(ncol=2)) + ylim(0.5,1)


############### Figure S4D ############### 
load('./case2/r_object/sp/all/all.RData')

# Make the boxplot of the distribution of Bacteroides xylanisolvens abundance
p1 <- ggplot() + geom_boxplot(aes(country, X[,'Bacteroides xylanisolvens'],col = as.factor(Y)))+
  geom_jitter(aes(country, X[,'Bacteroides xylanisolvens'],col = as.factor(Y)),
              position=position_dodge2(0.8))+
  scale_color_manual('Condition', labels =c('Normal','CRC'), values = c('black','orange'))+
  ylab('Bacteroides xylanisolvens')+
  xlab('') + 
  theme(text = element_text(size= 15), 
        axis.text.x = element_text(angle = 45, vjust = 0.5,size = 0))

# Make the boxplot of the distribution of Prevotella copri abundance
p2 <- ggplot() + geom_boxplot(aes(country, X[,'Prevotella copri'],col = as.factor(Y)),show.legend = F)+
  geom_jitter(aes(country, X[,'Prevotella copri'],col = as.factor(Y)),
              position=position_dodge2(0.8),show.legend = F)+
  scale_color_manual('Condition', labels =c('Normal','CRC'), values = c('black','orange'))+
  ylab('Prevotella copri')+
  xlab('Datasets') + 
  theme(text = element_text(size= 15), 
        axis.text.x = element_text(angle = 45, vjust = 0.5))

# Combine the two plots
ggarrange(p1,p2, align = 'v', ncol = 1, common.legend = T, legend = 'right')


############### Figure S6A ############### 
load('./case2/r_object/sp/all/all.RData')

# Weight samples by cohort-case (disease status) sizes
n_sample <- length(Y)
agg_group <- paste(Y, envs)
w <- 1/(table(agg_group)[agg_group])
w <- w/sum(w)*n_sample

# Train RF model 
rf_mod <- randomForest::randomForest(X, as.factor(Y), weights = w)

# Train Lasso model 
lasso_mod <- glmnet(X, Y, family = 'binomial', weights = as.numeric(w), nlambda = 1000)
sel.matrix <- (lasso_mod$beta != 0) # Indicator matrix that record lasso selection path
first.entrance <- apply(sel.matrix, 1, which.max) # Calculate at which lambda each predictor is first selected
first.entrance[which(apply(sel.matrix, 1, sum) == 0)] <- Inf # If not selected in the entire path, set the selection order to infinity

# Create a matrix storing importance scores of predictors obtained from each method. 
score_mat <- data.frame(matrix(ncol = 3, nrow = length(species)))
colnames(score_mat) <- c('SM','Lasso','RF')
rownames(score_mat) <- species
score_mat$SM <- mod$stable_ensemble$selection_prob$marginal$selection_prob[-1]
score_mat$Lasso <- -first.entrance
score_mat$RF <- rf_mod$importance[,1]

# Rank importance scores within methods
rank_mat <- data.frame(apply(-score_mat, 2, rank, ties.method = 'max'))
rank_mat <- rank_mat[apply(rank_mat,1,min) <= 50,] # Only retain species that are ranked top 50 in at least one cohort
rank_mat$Mean <- round(rowMeans(rank_mat),2) # Calculate mean ranking of species
rank_mat <- rank_mat[order(rank_mat$Mean),] # Order species by their mean ranking

# Make the heatmap
pheatmap(rank_mat, color = c(colorRampPalette(brewer.pal(n = 4, name = "Blues")[1:2])(50),
                             colorRampPalette(brewer.pal(n = 4, name = "Blues")[3:4])(max(rank_mat)-50)),
         display_numbers = rank_mat,cluster_rows = F, cluster_cols = F, angle_col = 45, fontsize = 12)

