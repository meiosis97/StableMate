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
  
  # A matrix storing mean BIC scores of selections over variable selection ensembles
  BIC_record <- data.frame(st2e = NA, SM = NA) 
  # A matrix storing running time
  time_record <- data.frame(st2e = NA, SM = NA, SR = NA)
  # A matrix storing AUC scores calculated on importance scores of predictors in distinguishing predictors in Markov blankets 
  AUC_record <- data.frame(st2e = NA, SM = NA, SR = NA)
  
  ############### Simulation ###############
  idxMb = NULL
  while(length(idxMb)==0){
    {  nNode <- 100 # Number of variables
    propE <- 0.05 # Proportion of edges being connected
    nI <- 60 # Number of intervention
    
    dag <- rgraph(nNode, propE, dag = T, ordered = TRUE) # Sample a random graph of 50 variables
    mat <- igraph::as_adjacency_matrix(dag) # Extract the adjacency matrix of the graph
    mat <- rbind(matrix(ncol = nNode, nrow = nI,0),mat) # Expand the adjacency matrix by including interventions
    mat <- cbind(matrix(ncol = nI, nrow = nrow(mat),0),mat)
    
    
    nNode <- nrow(mat) # Recalculate the number of variables after inclusion of interventions
    idxY <- nI + 50 # Set the id of the response
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
  # Benchmark methods' selection accuracy in selecting predictors in Markov blankets
  # Identification of Markov blankets can be done in one environment, hence we 
  # only train models in the fourth environment (the testing environment)
  
  N <- length(testY) # Number of samples
  X <- testX # Predictor matrix, which is collected from the forth environment
  colnames(X) <- paste('X', 1:ncol(testX), sep = '') # Clean up predictor names
  y <- matrix(testY) # Response
  # Identify predictors that genuinely belong to the Markov blanket
  MB <- paste('X',idxMb,sep='')
  MB <- paste('X', which(colnames(testX) %in% MB), sep = '')
  # Check for each predictor whether they belong to the Markov blanket or not
  inMB <- colnames(X) %in% MB
  
  
  ############### st2 ###############
  time_record[1, 'st2e'] <- system.time({
    # Model training
    mod_st2e <- stst(100, X, 
                     y, 0.5, exp(2.45))
    
  })[3]
  
  # Obtain selection frequencies of predictors, which are set as importance scores of predictors
  scores_st2e <- mod_st2e[[2]]
  BIC_record[1, 'st2e'] <- -mean(mod_st2e[[1]]) # Calculate mean BIC scores of selections
  AUC_record[1, 'st2e'] <- auc(inMB, scores_st2e) # Calculate AUC scores 
  
  
  ############### stablemate ###############
  time_record[1, 'SM'] <- system.time({
    # Model fitting
    mod_sm <- ST2E(y, X, K = 100,
                   subsample_method = 'all')
    
  })[3]
  
  # Obtain predictivity scores of predictors, which are set as importance scores of predictors
  # Predictivity scores are calculated based on Equation (1) in Supplementary Method 7.1.3
  scores_sm <- mod_sm$selection_prob$marginal$selection_prob[-1]
  BIC_record[1, 'SM'] <- -mean(mod_sm$scores) # Calculate mean BIC scores of selections
  AUC_record[1, 'SM'] <- auc(inMB, scores_sm) # Calculate AUC scores 
  
  
  ############### stabilized regression ###############
  time_record[1, 'SR'] <- system.time({
    # Model fitting
    mod_sr <- SRanalysis(X, y, NA, num_reps = 100, prescreen_types = 'lasso', pred_scores = 'bic',
                         pars = list(stab_test  = 'exact', B = 2000, prescreen_size = 20,
                                     alpha_stab = 1-exp(-3), alpha_pred = 1-exp(-3)), verbose = 1)
    
  })[3]
  
  # Obtain selection probability of predictors over the stability selection procedure embedded in
  # SRanalysis. Selection probabilities are set as importance scores of predictors.
  scores_sr <- mod_sr$results$SR$selection_probs
  AUC_record[1, 'SR'] <- auc(inMB, scores_sr) # Calculate AUC scores 
  
  list(BIC_record = BIC_record, time_record = time_record, AUC_record = AUC_record) # Save the result
  
}


############### Save data ###############
# Combine results obtained from parallel computing
BIC_record <- do.call(rbind, lapply(mod_list, function(x) x$BIC_record))
AUC_record <- do.call(rbind, lapply(mod_list, function(x) x$AUC_record))
time_record <- do.call(rbind, lapply(mod_list, function(x) x$time_record))

# Save the result
save.image('./simulation/r_object/mb_100.RData')