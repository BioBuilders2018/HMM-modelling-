

dat1=rpois(100,lambda = 10) # hyphae data main branch, lambda mean rate pr hour
dat2=rpois(0.05*100,lambda = 5/100) # hayphae side braches for main branch, lambda mean side 
                                    # braching rate pr hour 

# example : lambda=10 vs lambda=0.05 => braching is 200 times less likely to occur

dat=sample(c(dat1,dat2)) # randomly sample braching and none branching


setwd("E:/mathies/Desktop/stat teori") # set this to your own
source("A1.R") # Hmm likelihood function
y=dat

## 2 - state 
## Initial values
m <- 2 # number of states
lambda0 <- quantile(y,c(0.25,0.75)) # inital lambda guess
gamma0 <- matrix(0.05,ncol=m,nrow=m) # gamma matrix guess
diag(gamma0) <- 1-(m-1)*gamma0[1,1]

fit2 <- pois.HMM.mle(dat,m,lambda0,gamma0)
fit2 # results lambda= state mean, gamma=transition matrix, delta=statinary distribution

## 3 - state 
## Initial values
m <- 3
lambda0 <- quantile(y,c(0.25,0.5,0.75))
gamma0 <- matrix(0.05,ncol=m,nrow=m)
diag(gamma0) <- 1-(m-1)*gamma0[1,1]

## optimize
fit3 <- pois.HMM.mle(y,m,lambda0,gamma0)
fit3


## 4 - state 
## Initial values
m <- 4
lambda0 <- quantile(y,c(0.2,0.4,0.6,0.8))
gamma0 <- matrix(0.025,ncol=m,nrow=m)
diag(gamma0) <- 1-(m-1)*gamma0[1,1]

## optimize
fit4 <- pois.HMM.mle(y,m,lambda0,gamma0)
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

pois.HMM.lalphabeta <-function(x,m,lambda ,gamma ,delta=NULL){
  if(is.null(delta))delta <-solve(t(diag(m)-gamma+1),rep(1,m))
  n <- length(x)
  lalpha <- lbeta <- matrix(NA,m,n)
  allprobs <- outer(x,lambda ,dpois)
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
  list(la=lalpha ,lb=lbeta)
}  

pois.HMM.state_probs <-function(x,m,lambda ,gamma ,delta=NULL ,...) {
  if(is.null(delta))delta <-solve(t(diag(m)-gamma+1),rep(1,m))
  n <- length(x)
  fb <- pois.HMM.lalphabeta(x,m,lambda ,gamma , delta=delta)
  la <- fb$la
  lb <- fb$lb
  c <- max(la[,n])
  llk <- c+log(sum(exp(la[,n]-c)))
  stateprobs <- matrix(NA ,ncol=n,nrow=m)
  for (i in 1:n) stateprobs[,i]<-exp(la[,i]+lb[,i]-llk)
  stateprobs
}



pois.HMM.local_decoding <-function(x,m,lambda ,gamma ,delta=NULL ,...){
  stateprobs <- pois.HMM.state_probs(x,m,lambda ,gamma ,delta=delta)
  n=length(x)
  ild <- rep(NA,n)
  for (i in 1:n) ild[i]<-which.max(stateprobs[,i])
  ild
}

pois.HMM.lalphabeta(x=y,m=2,lambda=fit2$lambda ,gamma=fit2$gamma) 

pois.HMM.state_probs(x=y,m=2,lambda=fit2$lambda ,gamma=fit2$gamma)

states=pois.HMM.local_decoding(x=dat,m=2,lambda = fit2$lambda,gamma = fit2$gamma)

X=data.frame(states)
cbind(dat,X)
lb=1:length(dat) # braching length
k1=12 # kinetic param k1
k2=18 # kinetic param k2
Kt=15 # kinetic param Kt
s0= 0 # starting point  

pattern <- function(states) {
  Q=c()
  if (any(states==2)) {
    for (i in 1:length(X$states)){
      Q[i]=k1+k2*lb[i]/(lb[i]+Kt)
    }
    
    b1=c()
    for (i in 1:length(X$states)){
      b1[i]= s0 + Q[i]
    }
    
  } 
  
  budd=which(1 == X)
  
  ###############################
  Y=matrix(nrow = length(states), ncol=length(budd)) # budding matrix for alle buds
  
  for (i in 1:length(budd)) {
    Y[budd[i],i]=b1[budd[i]]
  }
y=data.frame(b1,Y)
return(y)
}

br=pattern(states = X$states) # identifies braching points

plot(lb,br$b1,ylim=c(0,50),type="l",ylab = "extension rate unit/h", xlab="brach length unit",
     main="hyphae extesion rate")
lines(br$X1,type="p",col="2")
lines(br$X2,type="p",col="3")
lines(br$X3,type="p",col="4")
legend("topright", inset=.02,
       c("1st branch","2nd branch", "3rd brach"), col=c(2,3,4),lty=c(NA,1,2),pch=c(1,1,1),lwd =c(NA,NA,NA),
       cex=1,ncol=1)

