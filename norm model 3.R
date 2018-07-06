setwd("C:/Users/mathies/Desktop/igem HMM")
source("A3.R") # Hmm likelihood functions

# generate random data from guessing
gamma= rbind(c(0.7,0.3),c(0.7,0.3)) # state transition matrix
sigma=c(2,1) # sigma2= std of hyphae braching expansion length, sigma1= std of -"- main branch 
mu=c(4,2) # mu1= main branching mean expansion length, mu2= branching expasion length 


y=norm.HMM.generate_sample(n=20,m=2,mu=mu,sigma=sigma,gamma=gamma)

## 2 - state 
## Initial values
m <- 2 # number of states
mu0 <- quantile(y,c(0.25,0.75)) # inital mu guess
sigma0<-rep(sd(y),m)
gamma0 <- matrix(0.05,ncol=m,nrow=m) # gamma matrix guess
diag(gamma0) <- 1-(m-1)*gamma0[1,1]
delta= statdist(gamma0) # stat dist of the gamma matrix best guess

fit2=norm.HMM.EM(x=y,m,mu0,sigma0,gamma0,delta)

fit2 # results lm= state mean,sd=state std, gamma=transition matrix, delta=statinary distribution

## 3 - state 
## Initial values
m <- 3
mu0 <- quantile(y,c(0.25,0.5,0.75))
gamma0 <- matrix(0.05,ncol=m,nrow=m)
sigma0<-rep(sd(y),m)
diag(gamma0) <- 1-(m-1)*gamma0[1,1]
delta=statdist(gamma0)

## optimize
fit3=norm.HMM.EM(x=y,m,mu0,sigma0,gamma0,delta)
fit3


## 4 - state 
## Initial values
m <- 4
mu0 <- quantile(y,c(0.2,0.4,0.6,0.8))
gamma0 <- matrix(0.025,ncol=m,nrow=m)
sigma0<-rep(sd(y),m)
diag(gamma0) <- 1-(m-1)*gamma0[1,1]
delta=statdist(gamma0)

## optimize
fit4 <- norm.HMM.EM(x=y,m,mu0,sigma0,gamma0,delta)
fit4

AIC <- c(fit2$AIC,fit3$AIC,fit4$AIC)
ll <-  -c(fit2$mllk,fit3$mllk,fit4$mllk)
m <- c(2,3,4)
df <- m + (m^2-m) ## lambda + gamma

# report  no of states, df, AIC and negative loglikelihood 
res=data.frame(m, df, AIC, ll) 
names(res)=c("states","df","AIC","likelihood")

print(res)

print(fit2) # best fit

# prediction

states0=norm.HMM.local_decoding(x=y,m=2,mu=fit2$mu,sigma = fit2$sigma,gamma = fit2$gamma)

X=data.frame(states0)
cbind(y,X)

sidebranch <- function(states0,substate) {
  
  budd=which(states0==substate)
  Y=matrix(NA,nrow = length(states0),ncol=length(budd))
  
  
  
  
  for (i in 1:length(budd)){
    Y[,i]= c(rep(NA,length(states0)-length(norm.HMM.generate_sample(n=length(states0)-budd[i]+1,m=2,mu = fit2$mu,fit2$sigma,gamma = fit2$gamma)))
             ,norm.HMM.generate_sample(n=length(states0)-budd[i]+1,m=2,mu = fit2$mu,fit2$sigma,gamma = fit2$gamma))
    
  }
  
  return(y=data.frame(Y))
}

sidebranch_fit <- function(x,m) {
  # mu inital guess
  mu0=matrix(NA,nrow = length(x),ncol=2)
  for (i in 1:length(x)){
  mu0[i,1] <- quantile(na.omit(x[,i]),c(0.25))
  } 
  for (i in 1:length(x)){
    mu0[i,2] <- quantile(na.omit(x[,i]),c(0.75))
  }
  
  # sigma inital guess
  sigma0=matrix(NA,nrow = length(x),ncol=2)
  for (i in 1:length(x)){
  sigma0[i,1]=sd(na.omit(x[,i]))
  }
  for (i in 1:length(x)){
    sigma0[i,2]=sd(na.omit(x[,i]))
  }
  
  gamma0 <- matrix(0.05,ncol=2,nrow=2) # gamma matrix guess
  diag(gamma0) <- 1-(2-1)*gamma0[1,1]
  delta= statdist(gamma0)
  
  fit=list()
  for (i in 1:length(x)){
    fit[[i]]=norm.HMM.EM(na.omit(x[,i]),m=m,mu0[i,],sigma0[i,],gamma0,delta)
  }
  
  return(fit)
}


sbr=sidebranch(X$states0,substate = 2)
data.frame(y,sbr)
overall=sidebranch_fit(x=sbr,m=2) # antal fits= antal columns i substate dvs et fit pr column

sidebranch_states<- function(x,m){
  Y=matrix(NA,nrow = length(states0),ncol=length(sbr))
  
  for (i in 1:length(sbr)){
    Y[,i]=c(rep(NA,length(X$states0)-length(norm.HMM.local_decoding(x=na.omit(sbr[,i]),m=2,mu=overall[[i]]$mu, sigma=overall[[i]]$sigma,gamma=overall[[i]]$gamma))),
          norm.HMM.local_decoding(x=na.omit(sbr[,i]),m=2,mu=overall[[i]]$mu,sigma=overall[[i]]$sigma ,gamma=overall[[i]]$gamma))
  }
  
  
  return(y=data.frame(Y))
  
}

mainBr.to.sideBr1_states=sidebranch_states(x=sbr,m=2)


mainBr.sideBr1=data.frame(states0,mainBr.to.sideBr1_states)


sidebranch1 <- function(states,substate) {
  
  budd=list()
  for ( i in 1:length(states) ){
  budd[[i]]=which(states[,i]==substate)
  }
  
  #Y=matrix(NA,nrow = length(states0),ncol=length(budd[[1]])-1)
  Y=rep(list(list()),length(budd))
  if ( length(budd)==0){
    return(rep(NA,length(states)))
  }
  
  # matrix of all the generated new budding data
  for (k in 1:length(budd)){
  for (i in 1:(length(budd[[k]])-1)){
    Y[[k]][[i]]= c(rep(NA,length(states[,1])-length(norm.HMM.generate_sample(n=length(states0)-budd[[k]][i]+1,m=2,mu =  overall[[k]]$mu, overall[[k]]$sigma,gamma =  overall[[k]]$gamma)))
             ,norm.HMM.generate_sample(n=length(states[,1])-budd[[k]][i]+1,m=2,mu = overall[[k]]$mu, overall[[k]]$sigma,gamma =  overall[[k]]$gamma))
  }
  }
  
  Y
  return(y=(Y))
}

branch_list <- function(states0,substate) {

sbr=sidebranch1(states0,substate) 
  
sbr_list=list()
  for (i in 1:length(sbr)){
    sbr_list[[i]]=sbr[,i]
  }
  return(sbr_list)
  
    
}


# return list of list with branching points for each branch
QQ=mainBr.to.sideBr1_states
sbr1=sidebranch1(QQ,substate=2)

# tomorrow fix states of side branchen !!!!!

# prediction passed the data info- work in progress


# do this later:

# viterbi algorithmen
norm.HMM.viterbi(na.omit(sbr[,3]),m=2,mu=overall[[3]]$mu,sd=overall[[3]]$sigma,gamma=overall[[3]]$gamma)

# state pred
norm.HMM.state_prediction(x,m=2,mu=fit2$mu,
                          sigma =fit2$sigma,gamma=fit2$gamma,H=100)


# other stuff - the mabye pile

lb=1:length(y) # braching length
k1=12 # kinetic param k1
k2=18 # kinetic param k2
Kt=15 # kinetic param Kt
s0= 0 # starting point  

pattern <- function(states) {
  Q=c()
  if (any(states==1)) {
    for (i in 1:length(X$states)){
      Q[i]=k1+k2*lb[i]/(lb[i]+Kt)
    }
    
    b1=c()
    for (i in 1:length(X$states)){
      b1[i]= s0 + Q[i]
    }
    
  } 
  
  budd=which(2 == X)
  
  ###############################
  Y=matrix(nrow = length(states), ncol=length(budd)) # budding matrix for alle buds
  
  for (i in 1:length(budd)) {
    Y[budd[i],i]=b1[budd[i]]
  }
  y=data.frame(b1,Y)
  return(y)
}

br=pattern(states = X$states0) # identifies braching points

plot(lb,br$b1,ylim=c(0,50),type="l",ylab = "extension rate unit/h", xlab="brach length unit",
     main="hyphae extesion rate")
lines(br$X1,type="p",col="2")
lines(br$X2,type="p",col="3")
lines(br$X3,type="p",col="4")
legend("topright", inset=.02,
       c("1st branch","2nd branch", "3rd brach"), col=c(2,3,4),lty=c(NA,1,2),pch=c(1,1,1),lwd =c(NA,NA,NA),
       cex=1,ncol=1)

