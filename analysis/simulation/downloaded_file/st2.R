###################################################
#
# R codes for Stochastic Stepwise Ensembles for Variable Selection
# Lu Xin, University of Waterloo
# June 2011
#
###################################################
#
# The R codes are separated to two parts.
# The first part contains the functions that implement the ST2E algorithm.
# The main function that implements the ST2E algorithm is stst().
# This function requires two other functions add() and del().
#
# The second part contains the example for how to use stst(), tuning example
# and the R codes for running the other examples in the paper.
#
###################################################


######################################################
#              
# First Part: functions
#
######################################################

make.formula<-function(binstr)
{
  # makes formula for regression
  # must make sure dimnames(X)[[2]] = "X1","X2",...
  ind<-which(binstr==1)
  rht <- paste("X", ind, sep="", collapse="+")
  ans<-paste("y ~ ", rht)
  return(as.formula(ans))
}


add<-function(V,V.l,adddata,add.gs,add.cn)
{
  # performs Forward Step of ST2E
  
  n.data<-dim(adddata)[1]
  p.data<-dim(adddata)[2]-1
  AIC.add<-10000000000000000000
  
  if(V.l==0)
  {
    groupsize<-sample(1:round(p.data*add.gs),1) 
    candidate<-choose(p.data,groupsize)^(1/add.cn)
    
    for(i.add in 1:candidate)
    {
      aa.a<-rep(0,p.data)
      ch.a<-sample(1:p.data,groupsize)
      aa.a[ch.a]<-1
      AIC.t<-BIC(lm(make.formula(aa.a),data=adddata))
      if(AIC.t<AIC.add)
      {
        AIC.add<-AIC.t
        c1.a<-ch.a
      }
    }
    
    return(c(AIC.add,length(c1.a),c1.a))
  }
  
  if(V.l>0)
  {
    groupsize<-sample(1:round((p.data-V.l)*add.gs),1) 
    candidate<-choose((p.data-V.l),groupsize)^(1/add.cn)
    
    for(i.add in 1:candidate)
    {
      aa.a<-rep(0,p.data)
      ch.a<-sample((1:p.data)[-V],groupsize)
      aa.a[c(V,ch.a)]<-1
      AIC.t<-BIC(lm(make.formula(aa.a),data=adddata))
      if(AIC.t<AIC.add)
      {
        AIC.add<-AIC.t
        c1.a<-c(V,ch.a)
      }
    }
    
    return(c(AIC.add,length(c1.a),c1.a))
  }
  
}


del<-function(V,deldata,del.gs,del.cn)
{
  #performs Backward Step of ST2E
  
  p.data<-dim(deldata)[2]-1
  AIC.del<-100000000000000000
  groupsize<-sample(1:round(length(V)*del.gs),1)
  
  if(groupsize==length(V))
  {
    return(c(BIC(lm(deldata[,p.data+1]~1)),0,0))
  }
  
  else if (groupsize<length(V))
  {
    candidate<-choose(length(V),groupsize)^(1/del.cn)
    for(i.del in 1:candidate)
    {
      aa.d<-rep(0,p.data)
      ch.d<-sample(V,length(V)-groupsize)
      aa.d[ch.d]<-1
      AIC.t<-BIC(lm(make.formula(aa.d),data=deldata))
      if(AIC.t<AIC.del)
      {
        AIC.del<-AIC.t
        c1.d<-ch.d
      }
    }
    
    return(c(AIC.del,length(c1.d),c1.d))
    
  }
}


stst<-function(T,X,Y,gs.lambda,cn.kappa)   
{
  
  # main function of ST2E 
  
  # The function returns a vector with length p+2, where p is the number
  # of predictors. The first number records the average within-ensemble
  # variation and the second number records the mean strength of the
  # ensemble (Section 2.3.2). The remaining p numbers record the number of   
  # times that the p variables are selected, respectively.  
  
  # T         = ensemble size
  # X         = predictor matrix
  # Y         = response vector
  # gs.lambda = parameter lambda (set to be 0.5 throughout this papaer)
  # cn.kappa  = parameter kappa
  
  stst.p<-dim(X)[2] 
  mydata<-data.frame(X,Y)
  
  forest<-rep(0,stst.p) 
  AIC.n<-BIC(lm(Y~1)) 
  k.aic<-rep(0,T)
  div<-matrix(0,T,stst.p)
  
  for(i in 1:T)
  {
    tree<-rep(0,stst.p)
    AIC<-AIC.n
    V<-0     
    V.l<-0   
    ind<-1   
    
    while(ind==1)
    {
      ind1<-0  
      ind2<-0  
      
      if(V.l<stst.p)
      {
        f<-add(V,V.l,mydata,gs.lambda,cn.kappa)
        
        if(f[1]<AIC)
        {
          
          AIC<-f[1]
          V.l<-f[2]
          V<-f[-(1:2)]
          ind1<-1
        }
        
      }
      
      if(V.l>0)
      {
        b<-del(V,mydata,gs.lambda,cn.kappa)
        
        if(b[1]<AIC)
        {    
          AIC<-b[1]
          V.l<-b[2]
          V<-b[-(1:2)]
          ind2<-1
        }
        
      }
      
      if(ind1==0&&ind2==0)
      {ind<-0}
      
    }
    
    k.aic[i]<-AIC
    if(V.l>0)
    {
      tree[V]<-1
    }
    forest<-forest+tree
    div[i,]<-tree
  }
  
  return(list(k.aic,forest))
  
}
