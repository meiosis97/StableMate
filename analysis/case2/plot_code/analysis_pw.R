############### The following preloading was used for the entire script (MUST RUN) ############### 
############### All the plot sessions can be ran independently of each others after running the preloading ############### 
require(dplyr)
require(ggpubr)
require(pheatmap)
require(RColorBrewer)
require(pROC)
require(rstatix)
setwd('/data/gpfs/projects/punim1662/yidi_projects/stablemate/')
source('./stablemate.R')
data_list <- readRDS('./case2/data/processed_data_pw.RDS')

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
pathway <- colnames(data_list[[1]])


############### Figure S5A ############### 
# Create a matrix for storing predictivity scores
score_mat <- data.frame(matrix(NA, ncol = length(pathway), nrow = length(country_table)))
# Create a indicator matrix for indicating predictive and stable predictors for each cohort.
indi_stabpred <- data.frame(matrix(0, ncol = length(pathway), nrow = length(country_table)))
# Create a indicator matrix for indicating predictive predictors for each cohort.
indi_pred <- data.frame(matrix(0, ncol = length(pathway), nrow = length(country_table)))
# Indexing the matrices
colnames(score_mat) <- colnames(indi_stabpred)  <- colnames(indi_pred) <- pathway
rownames(score_mat) <- rownames(indi_stabpred)  <- rownames(indi_pred)  <- country_table

for(c in names(country_table)){
  
  # Load stablemate results for each cohort
  dir <- paste('./case2/r_object/pw/one_by_one/', c, '.RData', sep = '')
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
load('./case2/r_object/pw/all/all.RData')

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

# Keep only the pathways that are selected as predictive by at least one cohort
rank_mat <- rank_mat[colSums(indi_pred) > 0,]
indi_stabpred <- indi_stabpred[,colSums(indi_pred) > 0]

# Order pathways by their mean ranking in predictivity scores
rank_mat <- rank_mat[order(rank_mat$Pooled),]

# Create labels for the heatmap
label_mat <- rank_mat
# Label a cell by asterisk if the corresponding pathway of the cell is selected as predictive and stable for the corresponding cohort
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
ggplot() + geom_col(aes(1:length(pathway),1:length(pathway),color = 1:length(pathway))) +
  geom_col(aes(1:4,1:4,fill = as.factor(1:4))) +
  theme(text = element_text(size = 15),legend.position = 'top')+
  scale_color_gradientn('Rank',colors = c(colorRampPalette(brewer.pal(n = 4, name = "Blues")[1:2])(nrow(rank_mat)),
                                          colorRampPalette(brewer.pal(n = 4, name = "Blues")[3:4])(max(rank_mat)-nrow(rank_mat))))+
  scale_fill_manual('Selection frequency',values = RColorBrewer::brewer.pal(5,'Greens'))


############### Figure S4C ############### 
load('./case2/r_object/pw/all/all.RData')

# Make the boxplot of the distribution of L-arginine biosynthesis III abundance
ggplot() + geom_boxplot(aes(country, X[,'PWY-5154: L-arginine biosynthesis III (via N-acetyl-L-citrulline)'], col = as.factor(Y)))+
  geom_jitter(aes(country, X[,'PWY-5154: L-arginine biosynthesis III (via N-acetyl-L-citrulline)'],col = as.factor(Y)),
              position=position_dodge2(0.8))+
  scale_color_manual('Sample type', labels =c('Normal','CRC'), values = c('black','orange'))+
  ylab('L-arginine biosynthesis III \n (via N-acetyl-L-citrulline)')+
  xlab('Dataset') + 
  theme(text = element_text(size= 15), 
        axis.text.x = element_text(angle = 45, vjust = 0.5), axis.title.y = element_text(hjust = 1),legend.position = 'right' )


############### Figure S4D ############### 
load('./case2/r_object/pw/one_by_one/FengQ_2015.RData')
# Make plot
plot.SRST2E(mod) +  theme(text = element_text(size = 15)) + ylim(0,1)


############### Figure S4B ############### 
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
  dir <- paste('./case2/r_object/pw/leave_one_out/no_', c, '.RData', sep = '')
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


############### Figure S6B ############### 
load('./case2/r_object/pw/all/all.RData')

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
score_mat <- data.frame(matrix(ncol = 3, nrow = length(pathway)))
colnames(score_mat) <- c('SM','Lasso','RF')
rownames(score_mat) <- pathway
score_mat$SM <- mod$stable_ensemble$selection_prob$marginal$selection_prob[-1]
score_mat$Lasso <- -first.entrance
score_mat$RF <- rf_mod$importance[,1]

# Rank importance scores within methods
rank_mat <- data.frame(apply(-score_mat, 2, rank, ties.method = 'max'))
rank_mat <- rank_mat[apply(rank_mat,1,min) <= 50,] # Only retain pathways that are ranked top 50 in at least one cohort
rank_mat$Mean <- round(rowMeans(rank_mat),2) # Calculate mean ranking of pathways
rank_mat <- rank_mat[order(rank_mat$Mean),] # Order pathways by their mean ranking

# Make the heatmap
pheatmap(rank_mat, color = c(colorRampPalette(brewer.pal(n = 4, name = "Blues")[1:2])(50),
                             colorRampPalette(brewer.pal(n = 4, name = "Blues")[3:4])(max(rank_mat)-50)),
         display_numbers = rank_mat,cluster_rows = F, cluster_cols = F, angle_col = 45, fontsize = 12)

