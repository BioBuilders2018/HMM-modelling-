statdist <- function(gamma){
  m = dim(gamma)[1]
  matrix(1,1,m) %*% solve(diag(1,m) - gamma + matrix(1,m,m))
}

norm.HMM.generate_sample <-function(n,m,mu,sigma,gamma,delta=NULL){
    if(is.null(delta))delta<-solve(t(diag(m)-gamma+1),rep(1,m))
    mvect <- 1:m
    state <- numeric(n)
    state[1] <- sample(mvect,1,prob=delta)
    for (i in 2:n)
      state[i]<-sample(mvect,1,prob=gamma[state[i-1],])
    x <- rnorm(n,mu[state],sigma[state])
    x
  }

norm.HMM.lalphabeta<-function(x,m,mu,sigma,gamma,delta=NULL){
  if(is.null(delta))delta<-solve(t(diag(m)-gamma+1),rep(1,m))
  n <- length(x)
  lalpha <- lbeta<-matrix(NA,m,n)
  # allprobs <- outer(x,lambda,dpois) # <- OLD
  allprobs <- matrix(NA,n,m) # <- NEW
  for(i in 1:m){ allprobs[,i] <- dnorm(x,mu[i],sigma[i])} # <- NEW
  foo <- delta*allprobs[1,]
  sumfoo <- sum(foo)
  lscale <- log(sumfoo)
  foo <- foo/sumfoo
  lalpha[,1] <- log(foo)+lscale
  for (i in 2:n)
  {
    foo <- foo%*%gamma*allprobs[i,]
    sumfoo <- sum(foo)
    lscale <- lscale+log(sumfoo)
    foo <- foo/sumfoo
    lalpha[,i] <- log(foo)+lscale
  }
  lbeta[,n] <- rep(0,m)
  foo <- rep(1/m,m)
  lscale <- log(m)
  for (i in (n-1):1)
  {
    foo <- gamma%*%(allprobs[i+1,]*foo)
    lbeta[,i] <- log(foo)+lscale
    sumfoo <- sum(foo)
    foo <- foo/sumfoo
    lscale <- lscale+log(sumfoo)
  }
  list(la=lalpha,lb=lbeta)
}

norm.HMM.EM <- function(x,m,mu,sigma,gamma,delta,maxiter=1000,tol=1e-6,...){
  n <- length(x)
  # lambda.next <- lambda # <- OLD
  mu.next <- mu # <- NEW
  sigma.next <- sigma # <- NEW
  gamma.next <- gamma
  delta.next <- delta
  for (iter in 1:maxiter)
  {
    # lallprobs <- outer(x,lambda,dpois,log=TRUE) # <- OLD
    lallprobs <- matrix(NA,n,m) # <- NEW
    for(i in 1:m){ lallprobs[,i] <- log(dnorm(x,mu[i],sigma[i]))} # <- NEW
    # fb <- pois.HMM.lalphabeta(x,m,lambda,gamma,delta=delta) # <- OLD
    fb <- norm.HMM.lalphabeta(x,m,mu,sigma,gamma,delta=delta) # <- NEW
    la <- fb$la
    lb <- fb$lb
    c <- max(la[,n])
    llk <- c+log(sum(exp(la[,n]-c)))
    for (j in 1:m)
    {
      for (k in 1:m)
      {
        gamma.next[j,k] <- gamma[j,k]*sum(exp(la[j,1:(n-1)]+
                                                lallprobs[2:n,k]+lb[k,2:n]-llk))
      }
      # lambda.next[j] <- sum(exp(la[j,]+lb[j,]-llk)*x)/ # <- OLD
      # sum(exp(la[j,]+lb[j,]-llk)) # <- OLD
      mu.next[j] <- sum(exp(la[j,]+lb[j,]-llk)*x)/ # <- NEW
        sum(exp(la[j,]+lb[j,]-llk)) # <- NEW
      sigma.next[j] <- sqrt(sum(exp(la[j,]+lb[j,]-llk) # <- NEW
                                *(x-mu.next[j])^2)/ # <- NEW
                              sum(exp(la[j,]+lb[j,]-llk))) # <- NEW
    }
    gamma.next <- gamma.next/apply(gamma.next,1,sum)
    delta.next <- exp(la[,1]+lb[,1]-llk)
    delta.next <- delta.next/sum(delta.next)
    crit <- sum(abs(mu-mu.next)) + # <- NEW
      sum(abs(sigma-sigma.next)) + # <- NEW
      # sum(abs(lambda-lambda.next)) + # <- OLD
      sum(abs(gamma-gamma.next)) +
      sum(abs(delta-delta.next))
    if(crit<tol)
    {
      np <- m*m+m-1
      AIC <- -2*(llk-np)
      BIC <- -2*llk+np*log(n)
      return(list(mu=mu,sigma=sigma,gamma=gamma,delta=delta,
                  mllk=-llk,AIC=AIC,BIC=BIC))
    }
    mu <- mu.next # <- NEW
    sigma <- sigma.next # <- NEW
    # lambda <- lambda.next # <- OLD
    gamma <- gamma.next
    delta <- delta.next
  }
  print(paste("No convergence after",maxiter,"iterations"))
  NA
}

norm.HMM.state_probs <- function(x,m,mu,sigma,delta=NULL,...){                                                           
    if(is.null(delta))delta<-solve(t(diag(m)-gamma+1),rep(1,m)) 
    n          <- length(x)                                    
    fb         <- norm.HMM.lalphabeta(x,m,mu,sigma,gamma,        
                                      delta=delta)                               
    la         <- fb$la                                        
    lb         <- fb$lb                                        
    c          <- max(la[,n])                                  
    llk        <- c+log(sum(exp(la[,n]-c)))                   
    stateprobs <- matrix(NA,ncol=n,nrow=m)                     
    for (i in 1:n) stateprobs[,i]<-exp(la[,i]+lb[,i]-llk)     
    stateprobs                                                
  }      
    
norm.HMM.local_decoding <-function(x,m,mu,sigma,gamma,delta=NULL,...){                                                           
    stateprobs <-                                              
    norm.HMM.state_probs(x,m,mu,sigma,gamma,delta=delta)     
    n=length(x)
    ild <- rep(NA,n)                                           
    for (i in 1:n) ild[i]<-which.max(stateprobs[,i])           
    ild                                                        
  }          

norm.HMM.viterbi<-function(x,m,mu,sd,gamma,delta=NULL,...){                                                           
  if(is.null(delta))delta<-solve(t(diag(m)-gamma+1),rep(1,m))   
  n         <- length(x)
  normprobs <- matrix(NA,n,m) # <- NEW
  for(i in 1:m){ normprobs[,i] <- log(dnorm(x,mu[i],sigma[i]))}
  #poisprobs <- outer(x,lambda,dpois)                         
  xi        <- matrix(0,n,m)                                
  foo       <- delta*normprobs[1,]                          
  xi[1,]    <- foo/sum(foo)                                 
  for (i in 2:n)                                             
  {                                                        
    foo    <- apply(xi[i-1,]*gamma,2,max)*normprobs[i,]      
    xi[i,] <- foo/sum(foo)                                   
  }                                                        
  iv<-numeric(n)                                             
  iv[n]     <-which.max(xi[n,])                             
  for (i in (n-1):1)                                         
    iv[i] <- which.max(gamma[,iv[i+1]]*xi[i,])               
  iv                                                         
}                                                          

norm.HMM.state_prediction <-function(x,m,mu,sigma,gamma,delta=NULL,H=1,...){                                                           
    if(is.null(delta))delta<-solve(t(diag(m)-gamma+1),rep(1,m))  
    n          <- length(x)                                    
    fb         <- norm.HMM.lalphabeta(x,m,                     
                                      mu,sigma,gamma,delta=delta)                  
    la         <- fb$la                                       
    c          <- max(la[,n])                                 
    llk        <- c+log(sum(exp(la[,n]-c)))                    
    statepreds <- matrix(NA,ncol=H,nrow=m)                     
    foo1       <- exp(la[,n]-llk)                              
    foo2       <- diag(m)                                      
    for (i in 1:H)                                             
    {                                                        
      foo2           <- foo2%*%gamma                           
      statepreds[,i] <- foo1%*%foo2                            
    }                                                        
    statepreds                                                 
  }      

norm.HMM.conditionals <-function(x,m,mu,sigma,gamma,delta=NULL,xrange=NULL,...){                                                           
    if(is.null(delta))                                        
      delta  <- solve(t(diag(m)-gamma+1),rep(1,m))             
    if(is.null(xrange))                                        
      xrange <-qnorm(0.001,min(mu),min(sigma)):                       
        qnorm(0.999,max(mu),max(sigma))                        
    n      <- length(x)                                        
    fb     <- norm.HMM.lalphabeta(x,m,mu,sigma,gamma,delta=delta)  
    la     <- fb$la                                            
    lb     <- fb$lb                                           
    la     <- cbind(log(delta),la)                             
    lafact <- apply(la,2,max)                                  
    lbfact <- apply(lb,2,max)                                  
    w      <- matrix(NA,ncol=n,nrow=m)                         
    for (i in 1:n)                                             
    {                                                        
      foo   <- (exp(la[,i]-lafact[i])%*%gamma)*               
        exp(lb[,i]-lbfact[i])                          
      w[,i] <- foo/sum(foo)                                    
    }               
    allprobs <- matrix(NA,n,m) # <- NEW
    for(i in 1:m){ allprobs[,i] <- log(dnorm(x,mu[i],sigma[i]))}
    cdists   <- allprobs%*%w                                   
    list(xrange=xrange,cdists=cdists)                         
  }          


# example data
# case m=2
m= 2
mu = c(-3,3)
sigma = c(1,2)
gamma= rbind(c(0.9,0.1),c(0.1,0.9))
delta= statdist(gamma)
n=200
# Generate synthetic data
x <- norm.HMM.generate_sample(n,m,mu,sigma,gamma,delta)
# Fit normal-HMM with the EM algorith,
res <- norm.HMM.EM(x,m,mu,sigma,gamma,delta)    
res


