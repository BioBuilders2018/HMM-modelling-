# binominal HMM
statdist <- function(gamma){
  m = dim(gamma)[1]
  matrix(1,1,m) %*% solve(diag(1,m) - gamma + matrix(1,m,m))
}

binom.HMM.generate_sample <-function(n,m,N,p,gamma,delta=NULL){                                                           
  if(is.null(delta))delta<-solve(t(diag(m)-gamma+1),rep(1,m))    
  mvect <- 1:m                                               
  state <- numeric(n)                                        
  state[1] <- sample(mvect,1,prob=delta)                     
  for (i in 2:n)                                             
    state[i]<-sample(mvect,1,prob=gamma[state[i-1],])        
  x <- rbinom(n,size=N[state],prob=p[state])                      
  x                                                          
} 
 
m= 2
N = c(2,2)
p = c(0.9,0.1)
gamma= rbind(c(0.9,0.1),c(0.1,0.9))
delta= statdist(gamma)
n=200

x=binom.HMM.generate_sample(n=n,m=m,N=N,p=p,gamma = gamma)

binom.HMM.lalphabeta<-function(x,m,N,p,gamma,delta=NULL){
  if(is.null(delta))delta<-solve(t(diag(m)-gamma+1),rep(1,m))
  n <- length(x)
  lalpha <- lbeta<-matrix(NA,m,n)
  # allprobs <- outer(x,lambda,dpois) # <- OLD
  allprobs <- matrix(NA,n,m) # <- NEW
  for(i in 1:m){ allprobs[,i] <- dbinom(x,N[i],p[i])} # <- NEW
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

binom.HMM.lalphabeta(x,m,N,p,gamma)

binom.HMM.state_probs <-function(x,m,N,p,gamma,delta=NULL,...){                                                           
    if(is.null(delta))delta<-solve(t(diag(m)-gamma+1),rep(1,m)) 
    n          <- length(x)                                    
    fb         <- binom.HMM.lalphabeta(x,m,N,p,gamma,delta=delta)                              
    la         <- fb$la                                        
    lb         <- fb$lb                                        
    c          <- max(la[,n])                                  
    llk        <- c+log(sum(exp(la[,n]-c)))                   
    stateprobs <- matrix(NA,ncol=n,nrow=m)                     
    for (i in 1:n) stateprobs[,i]<-exp(la[,i]+lb[,i]-llk)     
    stateprobs                                                
  }    

binom.HMM.local_decoding <-function(x,m,N,p,gamma,delta=NULL,...){                                                           
    stateprobs <-                                              
    binom.HMM.state_probs(x,m,N,p,gamma,delta=delta)     
    n=length(x)
    ild <- rep(NA,n)                                           
    for (i in 1:n) ild[i]<-which.max(stateprobs[,i])           
    ild                                                        
  }       

binom.HMM.state_prediction <- function(x,m,N,p,gamma,delta=NULL,H=1,...){                                                           
    if(is.null(delta))delta<-solve(t(diag(m)-gamma+1),rep(1,m))  
    n          <- length(x)                                    
    fb         <- binom.HMM.lalphabeta(x,m,                     
                                      N,p,gamma,delta=delta)                  
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

binom.HMM.EM <- function(x,m,N,p,gamma,delta,maxiter=1000,tol=1e-6,...){
  n <- length(x)
  # lambda.next <- lambda # <- OLD
  p.next <- p # <- NEW
  gamma.next <- gamma
  delta.next <- delta
  for (iter in 1:maxiter)
  {
    #lallprobs <- outer(x,lambda,dpois,log=TRUE) # <- OLD
    #lallprobs <- outer(x,N,p,dbinom,log=TRUE) # <- OLD
    lallprobs <- matrix(NA,n,m) # <- NEW
    for(i in 1:m){ lallprobs[,i] <- log(dbinom(x,N,p[i]))} # <- NEW
    # fb <- pois.HMM.lalphabeta(x,m,lambda,gamma,delta=delta) # <- OLD
    fb <- binom.HMM.lalphabeta(x,m,N,p,gamma,delta=delta) # <- NEW
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
      p.next[j] <- sum(exp(la[j,]+lb[j,]-llk)*x)/ # <- NEW
        sum(exp(la[j,]+lb[j,]-llk)*n) # <- NEW
    }
    gamma.next <- gamma.next/apply(gamma.next,1,sum)
    delta.next <- exp(la[,1]+lb[,1]-llk)
    delta.next <- delta.next/sum(delta.next)
    crit <- sum(abs(p-p.next)) + # <- NEW
      #sum(abs(N-N.next)) + # <- NEW
      # sum(abs(lambda-lambda.next)) + # <- OLD
      sum(abs(gamma-gamma.next)) +
      sum(abs(delta-delta.next))
    if(crit<tol)
    {
      np <- m*m+m-1
      AIC <- -2*(llk-np)
      BIC <- -2*llk+np*log(n)
      return(list(N=N,p=p,gamma=gamma,delta=delta,
                  mllk=-llk,AIC=AIC,BIC=BIC))
    }
    #N <- N.next # <- NEW
    p <- p.next # <- NEW
    # lambda <- lambda.next # <- OLD
    gamma <- gamma.next
    delta <- delta.next
  }
  print(paste("No convergence after",maxiter,"iterations"))
  NA
}

binom.HMM.viterbi<-function(x,m,N,p,gamma,delta=NULL,...){                                                           
  if(is.null(delta))delta<-solve(t(diag(m)-gamma+1),rep(1,m))   
  n         <- length(x)
  binomprobs <- matrix(NA,n,m) # <- NEW
  for(i in 1:m){ binomprobs[,i] <- log(dbinom(x,N[i],p[i]))}
  #poisprobs <- outer(x,lambda,dpois)                         
  xi        <- matrix(0,n,m)                                
  foo       <- delta*binomprobs[1,]                          
  xi[1,]    <- foo/sum(foo)                                 
  for (i in 2:n)                                             
  {                                                        
    foo    <- apply(xi[i-1,]*gamma,2,max)*binomprobs[i,]      
    xi[i,] <- foo/sum(foo)                                   
  }                                                        
  iv<-numeric(n)                                             
  iv[n]     <-which.max(xi[n,])                             
  for (i in (n-1):1)                                         
    iv[i] <- which.max(gamma[,iv[i+1]]*xi[i,])               
  iv                                                         
}                                                          


binom.HMM.viterbi(x,m,N,p,gamma)  

binom.HMM.EM(x,m,N,p,gamma,delta)
binom.HMM.state_prediction(x,m,N,p,gamma,H=20)
binom.HMM.local_decoding(x,m,N,p,gamma)

# the mabye pile
hmmBinom(x,size=N, prob=p)
library(HiddenMarkov)
Pi=matrix(c(0.7,0.3,0.2,0.8),2,2,byrow=TRUE); delta=c(0.3,0.7)
n=100; pn=list(size=rep(5,n)); pm=list(prob=c(0.3,0.8))
myhmm=dthmm(NULL,Pi=gamma,delta=delta,distn="binom",pn=N,pm=p)
myhmm
x=simulate(myhmm,n)
mod=BaumWelch(x)
mod

Viterbi(x)


