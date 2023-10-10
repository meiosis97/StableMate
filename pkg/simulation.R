sim_scm <- function(nSample = 300, nE = 4, nNode = 50, propE = 0.1, nI = 30, idxY = 25){
  require(gmat)
  require(igraph)
  require(Matrix)
  idxSb = NULL
  idxMb= NULL 
  while(length(idxSb)==0 | length(idxSb) == length(idxMb)){
    { nNode <- nNode # Number of variables
    propE <- propE # Proportion of edges being connected
    nI <- nI # Number of intervention
    
    dag <- rgraph(nNode, propE, dag = T, ordered = TRUE) # Sample a random graph of 50 variables
    mat <- igraph::as_adjacency_matrix(dag) # Extract the adjacency matrix of the graph
    mat <- rbind(matrix(ncol = nNode, nrow = nI,0),mat) # Expand the adjacency matrix by including interventions
    mat <- cbind(matrix(ncol = nI, nrow = nrow(mat),0),mat)
    
    
    nNode <- nrow(mat) # Recalculate the number of variables after inclusion of interventions
    idxY <- nI + idxY # Set the id of the response
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
    N <- nSample # Number of samples in each environment
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
    trainX <- do.call(rbind,data[1:(nE-1)])[,-c(1:nI,idxY)] 
    colnames(trainX) <- paste('X', colnames(trainX), sep = '')
    trainY <- do.call(rbind,data[1:(nE-1)])[,idxY]
    group <- factor(rep(c('e1','e2','e3'),each =nSample)) # Environment variable
    
    # Data from the forth environment are used for testing
    testX <- data[[nE]][,-c(1:nI,idxY)]
    colnames(testX) <- paste('X', colnames(testX), sep = '')
    testY <- data[[nE]][,idxY]
    
    # Data from all environments
    allX <- rbind(trainX,testX)
    allY <- c(trainY, testY)
    allG <- factor(rep(paste('e', 1:nE, sep = ''),each =nSample)) # Environment variable
    
    }
  }
  
  
  # Plot the network with nodes colored by their class
  plot(dag, vertex.color = RColorBrewer::brewer.pal(6,'Set2')[nodeClass],
       vertex.size = 5,
       vertex.label.color="white",
       edge.width=0.5,                               
       edge.arrow.size=0.1,                            
       edge.arrow.width=1,                           
       edge.lty="solid",                            
       edge.curved=0.3,   vertex.label.cex = 0.01
  )
  
  return(list(trainX = trainX, trainY = trainY, testX = testX, testY = testY, envs = group, idxMb = idxMb, idxSb = idxSb, idxNsb = setdiff(idxMb, idxSb)))
  
}
