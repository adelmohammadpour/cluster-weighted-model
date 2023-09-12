########################################################################################
##     Random generation  for 2*2 (two cluster and two subcluster for each cluster)  ###
########################################################################################
RG3 <- function(n,muX,sdX,beta,sdY){
#n: vector of sample size for each subcluster   
#muX: mean vector of X variable
#sdX: variance covariance matrix of X   
#beta: coefficient vector for each subcluster 
#sdY: standard diviotion vector of Y|X for each subcluster 

  nn<-sum(n)
  pi1<-(n[1]+n[2])/(nn)
  pi11<-n[1]/(n[1]+n[2])
  pi22<-n[3]/(n[3]+n[4])
  zg0<-rbinom(nn,1,1-pi1)-2       # cluster 1 & 2   -2==1  -1==2
  zg1<-zg0
  for (ig1 in 1:nn){
    ifelse(zg1[ig1]==-2,zg1[ig1]<-rbinom(1,1,1-pi11)+1,zg1[ig1]<-rbinom(1,1,1-pi22)+1+2)  
  }

  X<-matrix(0,nn,3)
  Y<-matrix(0,nn,1)
  sdmv<-sdX
   muXmv<-muX
  for (ig2 in 1:nn){
    if (zg1[ig2]==1){X[ig2,]<-mvrnorm(1,mu=muXmv[zg1[ig2],],sdmv)}
    if (zg1[ig2]==2){X[ig2,]<-mvrnorm(1,mu=muXmv[zg1[ig2],],sdmv)}
    if (zg1[ig2]==3){X[ig2,]<-mvrnorm(1,mu=muXmv[zg1[ig2],],sdmv)}
    if (zg1[ig2]==4){X[ig2,]<-mvrnorm(1,mu=muXmv[zg1[ig2],],sdmv)}
  } 
  for (i1 in 1:nn){
    Y[i1]<-rnorm(1,mean=(beta[zg1[i1],1]*X[i1,1]+beta[zg1[i1],2]*X[i1,2]+beta[zg1[i1],3]*X[i1,3]),sd=sdY[zg1[i1]])
  }
  
 
zt2<-matrix(0,nn,4)
for (i1 in 1:nn){
  zt2[i1,zg1[i1]]<-1
}

  return(
    list(
      X     = X,
      Y     = Y,
      zg0  = zg0,
      zg1  = zg1,
      zt2   =zt2,
      muX   =muX,
      beta  =beta
    )
  )
  
}
