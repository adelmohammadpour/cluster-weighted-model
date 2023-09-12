#################################################################################################################
##      Mixture of Gaussian Cluster-Weighted Model for 2*2 (two cluster and two subcluster for each cluster)  ###
#################################################################################################################


FUMCWM <- function(
  Y,			                        # Y: a numerical vector for the response variable
  X,			                        # X: a numerical vector for the unique covariate
  initial = "random.soft",        # initial values: "random.soft", "random.hard", "manual" , "kmeans"
  threshold = 1.0e-02,            # stopping value 
  lab                             # lab: a numerical vector for the lable of observation
)
{
  ckv<-c(2,2)               
  iter.max = 1000
  k<-length(ckv)               # Number of cluster
  ck<-sum(ckv)                 # Number of sub cluster
  n <- length(Y)               # sample size

  muX1 <- array(0,c(k),dimnames=list(paste("comp.",1:k,sep="")))       # mean of X in each external group. main group
  sdX1 <- array(0,c(k),dimnames=list(paste("comp.",1:k,sep="")))       # sd of X in each external group
  sdX2 <- array(0,c(ck),dimnames=list(paste("comp.",1:ck,sep="")))     # sd of X in each internal group
  PX2  <- array(0,c(n,ck),dimnames=list(1:n,paste("comp.",1:ck,sep=""))) 
 
  betaT2 <- array(0,c(ck,2),dimnames=list(paste("comp.",1:ck,sep=""),paste("beta.",0:1,sep=""))) # (nk x (2)) matrix of coefficients (plus intercept) for each internal group
  sdY2   <- array(0,c(ck),dimnames=list(paste("comp.",1:ck,sep="")))
  PY2    <- array(0,c(n,ck),dimnames=list(1:n,paste("comp.",1:ck,sep="")))
  muY2   <- array(0,c(n,ck),dimnames=list(1:n,paste("comp.",1:ck,sep=""))) 
  muX2 <- array(0,c(ck),dimnames=list(paste("comp.",1:ck,sep="")))       # mean of X in each internal group. Total group

  
  ########  initial values start
  
  if(initial=="random.soft"){
  post <- array(runif(n*ck),c(n,ck))   # soft
  post <- post/rowSums(post)                 
  ztr    <- post
  }
  
  if(initial=="random.hard"){
   post <- t(rmultinom(n, size = 1, prob=rep(1/4,4)))  # hard                
    ztr    <- post
  }
  
  if(initial=="manual"){
    post <- initial              
    ztr    <- post
  }
  
  if(initial=="kmeans"){
    km1<-kmeans(Y,2,nstart=50)
    clu1<-km1$clus
    ay1<-Y[clu1==1]
    ay2<-Y[clu1==2]
    
    km11<-kmeans(ay1,2,nstart=50) 
    km12<-kmeans(ay2,2,nstart=50)
    clu22<-c(km11$clus,km12$clus+2)
    
    zkm<-matrix(0,n,ck)
    for (i1 in 1:n){
      zkm[i1,clu22[i1]]<-1
    }
    ztr    <- zkm
  } 
  
  ########  initial values End

  z3<- ztr
  z1<-matrix(0,n,2)
  z1[,1]<-z3[,1]+z3[,2]
  z1[,2]<-z3[,3]+z3[,4]
  
  prior1 <- colMeans(z1)
  prior2 <- c(colSums(z3[,1:2])/sum(z3[,1:2]),colSums(z3[,3:4])/sum(z3[,3:4]))
  prior3 <- colMeans(z3)
  
  # ------------ #
  # EM algorithm #
  # ------------ #
  
  # Preliminary definition of convergence criterions
  
  check     <- 0
  iteration <- 1
  loglik    <- NULL
  aloglik   <- NULL
  aloglik   <- c(0,0)
  a         <- NULL
  a         <- c(0,0)
ca<-0


  while(check<1){
    ca<-ca+1

    for(j in 1:ck){
      weights=z3[,j]
#       if(j<3){
#         weightsi=z3[,1]+z3[,2]
#         mu.obs    <- weighted.mean(x=X,w=weightsi)
#         covx <- sqrt(sum(weightsi*(X-mu.obs)^2)/sum(weightsi))
#       }
#       if(j>2){
#         weightsi=z3[,3]+z3[,4]
#         mu.obs    <- weighted.mean(x=X,w=weightsi) 
#         covx <- sqrt(sum(weightsi*(X-mu.obs)^2)/sum(weightsi))
#       }
      mu.obs    <- weighted.mean(x=X,w=weights)
      covx <- sqrt(sum(weights*(X-mu.obs)^2)/sum(weights))
#      mu.obs    <- c(weighted.mean(x=X[,1],w=weights),weighted.mean(x=X[,2],w=weights),weighted.mean(x=X[,3],w=weights))  
     # t(((weights*(X-mu.obs))%*%((weights*(X-mu.obs)/sum(weights)
#       va1<-sum(weights*(X[,1]-mu.obs[1])^2)/(sum(weights)-1)  # var x1
#       va2<-sum(weights*(X[,2]-mu.obs[2])^2)/(sum(weights)-1)  # var x2
#       va3<-sum(weights*(X[,3]-mu.obs[3])^2)/(sum(weights)-1)  # var x3
#       co1<-sum(weights*(X[,1]-mu.obs[1])*(X[,2]-mu.obs[2]))/(sum(weights)-1)  # cov 1,2
#       co2<-sum(weights*(X[,1]-mu.obs[1])*(X[,3]-mu.obs[3]))/(sum(weights)-1)  # cov 1,3
#       co3<-sum(weights*(X[,2]-mu.obs[2])*(X[,3]-mu.obs[3]))/(sum(weights)-1)  # cov 2,3
#       covx<-matrix(c(va1,co1,co2,co1,va2,co3,co2,co3,va3),3)
     # sdX2[j] <- sqrt(sum(weights*(X-mu.obs)^2)/sum(weights))
   #   muX2[j,]  <-   mu.obs
#      PX2[,j]  <- dmvnorm(X,mean=mu.obs, covx)    
       PX2[,j]  <- dnorm(X,mean=mu.obs, covx)
    }
    

    # --- #
    # Y|x #
    # --- #

    for(j in 1:ck){
      weights=z3[,j]
      modelY      <- lm(Y ~ X,weights=weights) 
      beta.obs    <- modelY$coefficients
      sigma.obs   <- summary(modelY)$sigma
      betaT2[j,]   <-  beta.obs
      sdY2[j]      <-  sigma.obs 
      muY2[,j]     <- beta.obs[1]+beta.obs[2]*X
 #    muY2[,j]     <- beta.obs[1]+beta.obs[2]*X[,1]+beta.obs[3]*X[,2]+beta.obs[4]*X[,3]
      PY2[,j]      <- dnorm(Y,mean=muY2[,j],sd=sdY2[j])              
    }

    # -------------- # 
    # log-likelihood # 
    # -------------- #

priortt<-c(rep(prior1[1],ckv[1]),rep(prior1[2],ckv[2]))
pt1<-matrix(rep(priortt,n),n,ck,byrow=TRUE)
pt2<-matrix(rep(prior2,n),n,ck,byrow=TRUE) 
   pt3<-pt1*pt2
 llvalues    <- sum(log(rowSums(pt3*PY2*PX2)))  
loglik[iteration] <- llvalues

    # ---------------------------------------------- #
    #               Aitken's criteria                #
    # ---------------------------------------------- #

    if(iteration>2){
      if (loglik[iteration]==loglik[iteration-1])check <- 1
      if (loglik[iteration]!=loglik[iteration-1]){
      a[iteration-1]      <- (loglik[iteration]-loglik[iteration-1])/(loglik[iteration-1]-loglik[iteration-2])     
aloglik[iteration]  <- loglik[iteration-1]+(1/(1-a[iteration-1])*(loglik[iteration]-loglik[iteration-1]))
 if(abs(aloglik[iteration]-loglik[iteration])<threshold){  
        check <- 1 }
    }
    }
    
    if(iteration==iter.max | ck==1) check <- 1
    iteration <- iteration + 1
    # ++++++ #
    # E-Step #
    # ++++++ #

    
  z.num31  <- matrix(c(prior1[1],prior1[1],prior1[2],prior1[2]),n,ck,byrow=TRUE)*matrix(rep(prior2,n),n,ck,byrow=TRUE)*PY2*PX2   # (n x ck)
  z.num3  <- z.num31/matrix(rep(rowSums(z.num31),4),ncol=4)  # (n x ck)
 
    ######## Start 3 Mode #####
    
    ####### S Mode 1 ## 12 34
    
    z.den31  <- rowSums(z.num31[,1:2])
    z.den32  <- rowSums(z.num31[,3:4])
    z2<- cbind(z.num31[,1:2]/matrix(rep(z.den31,2),ncol=2) ,z.num31[,3:4]/matrix(rep(z.den32,2),ncol=2))
    z.den311  <- rowSums(z.num31)
    z1<- cbind((rowSums(z.num31[,1:2]))/matrix(rep(z.den311,1),ncol=1) ,(rowSums(z.num31[,3:4]))/matrix(rep(z.den311,1),ncol=1))
    z2[z2=="NaN"]<-0      
    z31<-matrix(0,n,ck,byrow=TRUE)
      
      tt1<-apply(z1,1,which.max)
      tt21<-apply(z2[,1:2],1,which.max)
      tt22<-apply(z2[,3:4],1,which.max)
      tt1[tt1==2]<--1
      
      
      z31[,1][tt21*tt1==1]<-1
      z31[,2][tt21*tt1==2]<-1
      z31[,3][tt22*tt1==-1]<-1
      z31[,4][tt22*tt1==-2]<-1
      
      z11<-matrix(0,n,k,byrow=TRUE)
      z11[,1]<-z31[,1]+z31[,2]
      z11[,2]<-z31[,3]+z31[,4]
      
      prior11 <- colMeans(z11)
      prior21 <- c(colSums(z31[,1:2])/sum(z31[,1:2]),colSums(z31[,3:4])/sum(z31[,3:4]))
      prior31 <- colMeans(z31)
   
    priortt<-c(rep(prior11[1],ckv[1]),rep(prior11[2],ckv[2]))
    pt1<-matrix(rep(priortt,n),n,ck,byrow=TRUE)
    pt2<-matrix(rep(prior21,n),n,ck,byrow=TRUE)  
    pt3<-pt1*pt2
 
    llvalues1    <- sum(log(rowSums(pt3*PY2*PX2)))  
    
    ####### S Mode 2 ## 13 24
    
    z.den31  <- rowSums(z.num31[,c(1,3)])
    z.den32  <- rowSums(z.num31[,c(2,4)])
    z2<- cbind(z.num31[,c(1,3)]/matrix(rep(z.den31,2),ncol=2) ,z.num31[,c(2,4)]/matrix(rep(z.den32,2),ncol=2))
    
    z.den311  <- rowSums(z.num31)
    z1<- cbind((rowSums(z.num31[,c(1,3)]))/matrix(rep(z.den311,1),ncol=1) ,(rowSums(z.num31[,c(2,4)]))/matrix(rep(z.den311,1),ncol=1))
    z2[z2=="NaN"]<-0    
    
    
    z32<-matrix(0,n,ck,byrow=TRUE)

    tt1<-apply(z1,1,which.max)
    tt21<-apply(z2[,c(1,2)],1,which.max)
    tt22<-apply(z2[,c(3,4)],1,which.max)
    tt1[tt1==2]<--1
    
    
    z32[,1][tt21*tt1==1]<-1
    z32[,2][tt21*tt1==2]<-1
    z32[,3][tt22*tt1==-1]<-1
    z32[,4][tt22*tt1==-2]<-1
    
    z12<-matrix(0,n,k,byrow=TRUE)
    z12[,1]<-z32[,1]+z32[,2]
    z12[,2]<-z32[,3]+z32[,4]
    
    prior12 <- colMeans(z12)
    prior22 <- c(colSums(z32[,1:2])/sum(z32[,1:2]),colSums(z32[,3:4])/sum(z32[,3:4]))
    prior32 <- colMeans(z32)
    
    priortt<-c(rep(prior12[1],ckv[1]),rep(prior12[2],ckv[2]))
    pt1<-matrix(rep(priortt,n),n,ck,byrow=TRUE)
    pt2<-matrix(rep(prior22,n),n,ck,byrow=TRUE)  
    pt3<-pt1*pt2
    llvalues2    <- sum(log(rowSums(pt3*PY2*PX2))) 
    
    ####### S Mode 3 ## 14 23
    
    z.den31  <- rowSums(z.num31[,c(1,4)])
    z.den32  <- rowSums(z.num31[,c(2,3)])
    z2<- cbind(z.num31[,c(1,4)]/matrix(rep(z.den31,2),ncol=2) ,z.num31[,c(2,3)]/matrix(rep(z.den32,2),ncol=2))
    
    z.den311  <- rowSums(z.num31)
    z1<- cbind((rowSums(z.num31[,c(1,4)]))/matrix(rep(z.den311,1),ncol=1) ,(rowSums(z.num31[,c(2,3)]))/matrix(rep(z.den311,1),ncol=1))
    z2[z2=="NaN"]<-0    
    
    z33<-matrix(0,n,ck,byrow=TRUE)

    tt1<-apply(z1,1,which.max)
    tt21<-apply(z2[,c(1,2)],1,which.max)
    tt22<-apply(z2[,c(3,4)],1,which.max)
    tt1[tt1==2]<--1
    
    
    z33[,1][tt21*tt1==1]<-1
    z33[,2][tt21*tt1==2]<-1
    z33[,3][tt22*tt1==-1]<-1
    z33[,4][tt22*tt1==-2]<-1
    
    z13<-matrix(0,n,k,byrow=TRUE)
    z13[,1]<-z33[,1]+z33[,2]
    z13[,2]<-z33[,3]+z33[,4]
    
    prior13 <- colMeans(z13)
    prior23 <- c(colSums(z33[,1:2])/sum(z33[,1:2]),colSums(z33[,3:4])/sum(z33[,3:4]))
    prior33 <- colMeans(z33)
    
    priortt<-c(rep(prior13[1],ckv[1]),rep(prior13[2],ckv[2]))
    pt1<-matrix(rep(priortt,n),n,ck,byrow=TRUE)
    pt2<-matrix(rep(prior23,n),n,ck,byrow=TRUE)  
    pt3<-pt1*pt2
    llvalues3    <- sum(log(rowSums(pt3*PY2*PX2))) 
    ####### End 3 Mode #### 
    SM<-which.max(c( llvalues1 , llvalues2 , llvalues3 ))
    
    if(SM==1){
      z3<-z31
      z1<-z11
      prior1<-  prior11 
      prior2<-  prior21 
      prior3<- prior31 
    }
    if(SM==2){
      z3<-z32
      z1<-z12
      prior1<-  prior12 
      prior2<-  prior22 
      prior3<- prior32 
    }
    if(SM==3){
      z3<-z33
      z1<-z13
      prior1<-  prior13 
      prior2<-  prior23 
      prior3<- prior33 
    }
    
  } 

 finalloglik <- loglik[iteration-1] 

  # ----------------------------------------------------------------------- #
  #                  The EM-algorithm is finished                           #
  # ----------------------------------------------------------------------- #
  
Ze1<-apply(z1,1,which.max)
Ze2<-apply(z2,1,which.max)
Ze3<-apply(z3,1,which.max)
AR<-ARI(Ze3,lab)

  return(
    list(
      muX  =muX1,      # Mean of X in each sub cluster
      sdX   =sdX1,     # standard deviation X in each sub cluster
      Prior1=prior1,   # propbability for each cluster
      Prior2=prior2,   # propbability for each sub cluster
      Prior3=prior3,   # propbability for each sub cluster*propbability for each cluster
      Fll   = loglik,  # vector of likelihood values
      Fn    =finalloglik, # Final likelihood value
      Ze1    = Ze1,    # label observaion for cluster
      Ze2    =Ze2,     # label observaion for sub cluster in each cluster
      Ze3   =Ze3,      # label observaion for sub cluster
      AR    =AR$ari    # ARI criteria value
    )
  )
}
# library(MASS)
# library("cwm")
# library(bikm1)
# library(optimr)
# library(optimx)
# library("fossil")
# library("fitdistrplus")
# library(mvtnorm)

