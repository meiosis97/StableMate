# Load stablemate functions
# Some essential package has already been loaded
setwd('/data/gpfs/projects/punim1662/yidi_projects/stablemate')
source('./stablemate.R')
source('./simulation/downloaded_file/st2.R')

############### Parallel computing setting, to be ran in slurm system ###############
cl <- startMPIcluster(verbose=TRUE) # By default will start one fewer slave
# Than elements in tasks requested (to account for master)
registerDoMPI(cl)
niter <- 100


############### Iteration ###############
mod_list <- foreach(i=1:niter, .verbose = T) %dopar% {
  
  require(gmat)
  require(igraph)
  require(Matrix)
  require(splines)
  require(StabilizedRegression)
  require(pROC)
  require(progress)
  require(caret)
  require(glmnet)
  require(randomForest)
  
  # A matrix storing mean squared errors of predictions made in the testing environment
  MSE_record <- data.frame(OLS = NA, SM_stab = NA, SM_pred = NA, SR = NA, Lasso = NA, RF = NA)
  # A matrix storing running time
  time_record <- data.frame(SM = NA, SR = NA)
  # A matrix storing sensitivities in predicting stable blanket predictors
  sensSB_record <- data.frame(SM = NA, SR = NA)
  # A matrix storing specificities in predicting non-stable blanket predictors 
  specSB_record <- data.frame(SM = NA, SR = NA)
  # A matrix storing sensitivities in predicting stable blanket predictors 
  sensNSB_record <- data.frame(SM = NA, SR = NA)
  # A matrix storing specificities in predicting non-stable blanket predictors 
  specNSB_record <- data.frame(SM = NA, SR = NA)
  # A matrix storing F1 scores in predicting stable blanket predictors 
  F1SB_record <- data.frame(SM = NA, SR = NA)
  # A matrix storing F1 scores in predicting non-stable blanket predictors 
  F1NSB_record <- data.frame(SM = NA, SR = NA)
  # A matrix storing AUC scores calculated on importance scores of predictors in distinguishing predictors in stable blankets 
  AUC_record <- data.frame(SM = NA, SR = NA)
  
  
  ############### Simulation ###############
  idxSb = NULL
  idxMb= NULL 
  while(length(idxSb)==0 | length(idxSb) == length(idxMb)){
    {  nNode <- 50 # Number of variables
    propE <- 0.1 # Proportion of edges being connected
    nI <- 30 # Number of intervention
    
    dag <- rgraph(nNode, propE, dag = T, ordered = TRUE) # Sample a random graph of 50 variables
    mat <- igraph::as_adjacency_matrix(dag) # Extract the adjacency matrix of the graph
    mat <- rbind(matrix(ncol = nNode, nrow = nI,0),mat) # Expand the adjacency matrix by including interventions
    mat <- cbind(matrix(ncol = nI, nrow = nrow(mat),0),mat)
    
    
    nNode <- nrow(mat) # Recalculate the number of variables after inclusion of interventions
    idxY <- nI + 25 # Set the id of the response
    idxX <- (1:nNode)[-c(1:nI,idxY)] # Set the id of the predictors
    idxItv <- sample(idxX, nI) # Randomly chose predictors to intervene, the chosen predictors' id are recorded
    for(i in 1:nI) mat[i,idxItv[i]] = 1 # Adjust the adjacency matrix accordingly by connecting interventions with their targeted predictors
    dag <- graph_from_adjacency_matrix(mat) # Create a graph from the expanded adjancency matrix
    
    # Classify nodes
    
    idxPa <- which(mat[,idxY]==1) # Find causal parents of the response
    idxCh <- which(mat[idxY,]==1) # Find child nodes of the response
    # Find spouses of the response (parents of children)
    idxChPa <- c()
    if(length(idxCh)>0){
      
      for(i in idxCh) idxChPa <- c(idxChPa,which(mat[,i]==1)) # For each child, find its parent
      idxChPa <- unique(idxChPa)
      idxChPa <- idxChPa[idxChPa!=idxY]
      
    }
    idxMb <- unique(c(idxPa,idxCh,idxChPa)) # The Markov blanket (predictive predictors) includes parents, children and spouses
    idxMb <- idxMb[!idxMb%in%(1:nI)] # Interventions should not be in the Markov blanket
    
    # Find children that are intervened 
    idxChItv <- idxCh[colSums(matrix(mat[1:nI,idxCh] , nrow = nI))>0]
    
    # Find descendants of children that are intervened 
    idxChItvDec <- c()
    if(length(idxChItv)>0){
      
      for(i in idxChItv) idxChItvDec <- c(idxChItvDec,subcomponent(dag,i, mode = 'out')) # For each intervened child, find its descendants
      idxChItvDec <- unique(idxChItvDec)
      
    }
    # Find children that are not intervened and are also not descendants of intervened children 
    idxChUnItv <- idxCh[!idxCh%in%idxChItvDec] 
    
    # Find parents of unintervened children
    idxChUnItvPa <- c()
    if(length(idxChUnItv)>0){
      
      for(i in idxChUnItv) idxChUnItvPa <- c(idxChUnItvPa,which(mat[,i]==1)) # For each unintervened child, find its parents
      idxChUnItvPa <- unique(idxChUnItvPa)
      idxChUnItvPa <- idxChUnItvPa[idxChUnItvPa!=idxY]
      
    }
    
    # The stable blanket (predictive and stable predictors) includes parents, unintervened children and parents of unintervened children
    idxSb <- unique(c(idxPa,idxChUnItv,idxChUnItvPa))
    
    # Label nodes by whether they are parents, the Markov blanket or the stable blanket of the response
    nodeClass <- rep('Others',nNode)
    nodeClass[idxY] <- 'Response'
    nodeClass[idxX] <- 'Predictors'
    nodeClass[idxMb] <- 'MB'
    nodeClass[idxSb] <- 'SB'
    nodeClass[idxPa] <- 'PA'
    nodeClass[1:nI] <- 'Intervention'
    nodeClass <- factor(nodeClass, levels = c('Intervention','Predictors','Response', 'PA','SB','MB'))
    
    
    # Plot the network with nodes colored by their class
    plot(dag, vertex.color = RColorBrewer::brewer.pal(6,'Set2')[nodeClass],
         vertex.size = 3,
         vertex.label.color="white",
         edge.width=0.5,                               
         edge.arrow.size=0.1,                            
         edge.arrow.width=1,                           
         edge.lty="solid",                            
         edge.curved=0.3,   vertex.label.cex = 0.01
    )
    
    # Each variable is generated from a linear combination of another variable
    # Let Z be the vector of variable, then the data generation process of the system can be completely
    # described by the following equation: Z = BZ + U, where B is the matrix of path coefficients and U is the vector of 
    # random disturbance. Rearrange the equation we get Z = inv(I-B)U. Therefore, to generate Z, we just need to generate U
    # and multiply it with the matrix inv(I-B)
    
    # Create a matrix storing path coefficients.
    # Simulate path coefficients according to equation (13) in Supplementary Methods 7.2.1
    B <- mat
    B <- t(B)
    B[B!=0] <- sample(c(1,-1),sum(B), replace = T) * runif(sum(B),0.25,2)
    Binv <- as.matrix(solve(diag(nNode)-B))
    
    # Simulate data
    N <- 300 # Number of samples in each environment
    data <- list() # Create a list for storing data
    sigma <- runif(nNode,0.1,0.5) # Simulate standard deviations of U
    nE <- 4 # Number of environments
    mu <- matrix(nrow = nE, ncol = nNode) # Create a matrix for storing means of U
    # Refer to equation (15) in Supplementary Methods 7.2.1 for simulating means of U
    # Simulate means of disturbances of the response and predictors, which do not change across environments.
    muXY <-  runif(nNode - nI,-2,2) 
    # For each environment, simulate means of disturbances of interventions, which change across environments
    for(i in 1:nE){
      
      mu[i,] <- c(sample(c(1,-1),nI, replace = T) * runif(nI,2,10), 
                  muXY) 
      
    }
    # Generate data (Z) by generating random disturbances and multiply them with inv(I-B)
    for(i in 1:nE){
      
      data[[i]] <- t(Binv %*% replicate(N, rnorm(nNode, mu[i,], sigma)))
      colnames(data[[i]]) <- 1:nNode
      
    }
    
    # Split training data and testing data
    # Data from the first three environments are used for training
    trainX <- do.call(rbind,data[1:3])[,-c(1:nI,idxY)] 
    colnames(trainX) <- paste('X', colnames(trainX), sep = '')
    trainY <- do.call(rbind,data[1:3])[,idxY]
    group <- factor(rep(c('e1','e2','e3'),each =300)) # Environment variable
    
    # Data from the forth environment are used for testing
    testX <- data[[4]][,-c(1:nI,idxY)]
    colnames(testX) <- paste('X', colnames(testX), sep = '')
    testY <- data[[4]][,idxY]
    
    # Data from all environments
    allX <- rbind(trainX,testX)
    allY <- c(trainY, testY)
    allG <- factor(rep(c('e1','e2','e3','e4'),each =300)) # Environment variable
    
    }
  }
  
  
  ############### Model fitting ###############
  # Benchmark methods' selection accuracy in selecting predictors in stable blankets, 
  # as well as benchmark methods' prediction performace in extrapolatting the response 
  # in an unseen environment
  # Four environments are simulated, the first three are used as the training set, and the 
  # fourth is used as the testing set
  
  # Identify predictors that genuinely belong to the Markov blanket
  MB <- paste('X',idxMb,sep='')
  # Identify predictors that genuinely belong to the stable blanket
  SB <- paste('X',idxSb,sep='')
  # Identify predictors that genuinely belong to the non-stable blanket
  NSB <- MB[!MB%in%SB]
  # Check for each predictor whether they belong to the stable (non-stable) blanket or not
  inSB<- factor(colnames(trainX) %in% SB, level = c(FALSE,TRUE))
  inNSB<- factor(colnames(trainX) %in% NSB, level = c(FALSE,TRUE))
  
  
  ################################ OLS ######################################
  mod_ols <- lm(trainY ~., data = data.frame(trainX)) # Model training
  yhat_ols <- predict(mod_ols, data.frame(testX)) # Testing
  MSE_record[1, 'OLS'] <- -mean((testY-yhat_ols)^2) # Calculate MSE
  
  
  ################################ Lasso ######################################
  mod_lasso <- cv.glmnet(trainX, trainY, foldid = as.numeric(as.factor(group))) # Model training
  yhat_lasso <-  predict(mod_lasso, testX, s = mod_lasso$lambda.1se) # Testing
  MSE_record[1, 'Lasso'] <- -mean((testY-yhat_lasso)^2) # Calculate MSE
  
  
  ################################ RF ######################################
  mod_rf <- randomForest(trainX, trainY) # Model training
  yhat_rf <- predict(mod_rf, testX) # Testing
  MSE_record[1, 'RF'] <- -mean((testY-yhat_rf)^2) # Calculate MSE
  
  
  ################################ stablemate ######################################
  time_record[1, 'SM'] <- system.time({
    # Model training
    mod_sm <- SRST2E(trainY, trainX, group, K = 100,
                         subsample_method = 'all')
    
  })[3]
  
  # Obtain predictivity and stability scores of predictors, which are set as importance scores of predictors
  # Predictivity and stability scores are calculated based on Equation (1) in Supplementary Method 7.1.3
  scores_sm <- mod_sm$stable_ensemble$selection_prob$marginal$selection_prob[-1]
  AUC_record[1, 'SM'] <- auc(inSB, scores_sm)
  
  # Predict in the testing environment based on the SM model trained with stable predictors, calculate MSE
  yhat_sm_stab <- predict.SRST2E(mod_sm, testX)
  MSE_record[1, 'SM_stab'] <- -mean((testY-yhat_sm_stab)^2)
  
  # Predict in the testing environment based on the SM model trained without considering stable predictors, calculate MSE
  yhat_sm_pred <- predict.SRST2E(mod_sm, testX, assay = 'prediction_ensemble')
  MSE_record[1, 'SM_pred'] <- -mean((testY-yhat_sm_pred)^2)
  
  # Extract selections made by stablemate
  make_selection <- print.SRST2E 
  selection_sm <- make_selection(mod_sm, sigthresh = 1-exp(-3))
  # Check for each predictor whether they are selected into stable blanket (non-stable blanket) or not 
  inSB_sm <- factor(colnames(trainX) %in% 
                            selection_sm$SB_selected, level = c(FALSE,TRUE))
  inNSB_sm <- factor(colnames(trainX) %in% 
                             selection_sm$NSB_selected, level = c(FALSE,TRUE))
  # Calculate sensitivities and specificities of stable blanket and non-stable blanket selections
  F1SB_record[1, 'SM'] <-  confusionMatrix(table(inSB_sm, inSB), mode = "everything", positive="TRUE")$byClass['F1']
  F1NSB_record[1, 'SM'] <-  confusionMatrix(table(inNSB_sm, inNSB), mode = "everything", positive="TRUE")$byClass['F1']
  sensSB_record[1, 'SM'] <- sensitivity(table(inSB_sm, inSB), positive = 'TRUE')
  specSB_record[1, 'SM'] <-  specificity(table(inSB_sm, inSB), negative = 'FALSE')
  sensNSB_record[1, 'SM'] <- sensitivity(table(inNSB_sm, inNSB), positive = 'TRUE')
  specNSB_record[1, 'SM'] <-  specificity(table(inNSB_sm, inNSB), negative = 'FALSE')
  
  
  ################################## SR ########################################
  time_record[1, 'SR'] <- system.time({
    # Model training, SRanalysis is used for variable selection but not prediction
    mod_sr <- SRanalysis(trainX, trainY, group, num_reps = 100, prescreen_types = 'lasso', pred_scores = 'bic',
                         pars = list(stab_test  = 'exact', B = 1000, prescreen_size = 10,
                                     alpha_stab = exp(-3), alpha_pred = exp(-3))) 
    
  })[3]
  
  # Model traniing for prediction
  pred_mod_sr <- StabilizedRegression(trainX, trainY, group,
                                      pars = list(stab_test  = 'exact', B = 1000, prescreen_size = 10,
                                                  alpha_stab = exp(-3), alpha_pred = exp(-3),
                                                  prescreen_type = 'lasso', pred_score = 'bic'))
  
  # Obtain selection probability of predictors over the stability selection procedure embedded in
  # SRanalysis. Selection probabilities are set as importance scores of predictors.
  scores_sr <- mod_sr$results$SR$selection_probs
  AUC_record[1, 'SR'] <- auc(inSB, scores_sr) # Calculate AUC scores 
  
  # Predict in the testing environment, calculate MSE
  yhat_sr <- predict(pred_mod_sr, testX)
  MSE_record[1, 'SR'] <- -mean((testY-yhat_sr)^2)
  
  # Check for each predictor whether they are selected into stable blanket (non-stable blanket) or not 
  inSB_sr <- factor(mod_sr$results$SR$selection_probs > mod_sr$results$SR$siglevel & 
                      mod_sr$results$SRdiff$selection_probs < mod_sr$results$SRdiff$siglevel, level = c(FALSE,TRUE))
  inNSB_sr <- factor(mod_sr$results$SR$selection_probs < mod_sr$results$SR$siglevel & 
                       mod_sr$results$SRdiff$selection_probs > mod_sr$results$SRdiff$siglevel, level = c(FALSE,TRUE))
  # Calculate sensitivities and specificities of stable blanket and non-stable blanket selections
  F1SB_record[1, 'SR'] <-  confusionMatrix(table(inSB_sr, inSB), mode = "everything", positive="TRUE")$byClass['F1']
  F1NSB_record[1, 'SR'] <-  confusionMatrix(table(inNSB_sr, inNSB), mode = "everything", positive="TRUE")$byClass['F1']
  sensSB_record[1, 'SR'] <- sensitivity(table(inSB_sr, inSB), positive = 'TRUE')
  specSB_record[1, 'SR'] <-  specificity(table(inSB_sr, inSB), negative = 'FALSE')
  sensNSB_record[1, 'SR'] <- sensitivity(table(inNSB_sr, inNSB), positive = 'TRUE')
  specNSB_record[1, 'SR'] <-  specificity(table(inNSB_sr, inNSB), negative = 'FALSE')
  
  
  list(MSE_record = MSE_record, time_record = time_record, 
       sensSB_record = sensSB_record, specSB_record = specSB_record,
       sensNSB_record = sensNSB_record, specNSB_record = specNSB_record,
       F1SB_record = F1SB_record, F1NSB_record = F1NSB_record, AUC_record = AUC_record) # Save the result
  
}


############### Save data ###############
# Combine results obtained from parallel computing
MSE_record <- do.call(rbind, lapply(mod_list, function(x) x$MSE_record))
time_record <- do.call(rbind, lapply(mod_list, function(x) x$time_record))
sensSB_record <-  do.call(rbind, lapply(mod_list, function(x) x$sensSB_record))
specSB_record <- do.call(rbind, lapply(mod_list, function(x) x$specSB_record))
sensNSB_record <- do.call(rbind, lapply(mod_list, function(x) x$sensNSB_record))
specNSB_record <-  do.call(rbind, lapply(mod_list, function(x) x$specNSB_record))
F1SB_record <- do.call(rbind, lapply(mod_list, function(x) x$F1SB_record))
F1NSB_record <- do.call(rbind, lapply(mod_list, function(x) x$F1NSB_record))
AUC_record <- do.call(rbind, lapply(mod_list, function(x) x$AUC_record))

save.image('./simulation/r_object/sb_50.RData')
