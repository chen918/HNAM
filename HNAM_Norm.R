set.seed(99999)
if (!require(mvtnorm)) {
  install.packages("mvtnorm")
  library(mvtnorm)
}
if (!require(invgamma)) {
  install.packages("invgamma")
  library(invgamma)
}
if (!require(sna)) {
  install.packages("sna")
  library(sna)
}
if (!require(coda)) {
  install.packages("coda")
  library(coda)
}
if (!require(truncnorm)) {
  install.packages("truncnorm")
  library(truncnorm)
}

ddensity=0.8

rho=0.2

itl=500
medianbeta_m = matrix(0,nrow=3,ncol=itl)
medianzeta_m = matrix(0,nrow=4,ncol=itl)
medianSigma_v = c()
medianomega_v = c()
medianRho_v = c()
varRho_v = c()
coverage = 0
ubv = c()
lbv = c()



for (k in 1:itl){
  
  ################
  #Generating the data
  #level 2 model
  node=50 #Num of nodes
  X=matrix(rnorm(node*3),ncol=3)
  BBeta=matrix(c(1,1.5,2),nrow =3) #true values of beta
  oomega=1 #true value of omega
  Tao=matrix(rnorm(node, mean = 0, sd = sqrt(oomega)),ncol = 1)
  
  W=rgraph(node, m=1, tprob=ddensity, mode="graph", diag = F, replace=FALSE,
           tielist=NULL, return.as.edgelist= F )
  
  W[lower.tri(W)][which(W[lower.tri(W)]!=0)]=rgamma(n = length(which(W[lower.tri(W)]!=0)), shape = 0.1, scale = 2000)
  
  W[upper.tri(W)]<-t(W)[upper.tri(W)]
  
  W=t(apply(W, 1, function(x) (x/sum(x))))
  
  W[which(is.na(W))]=1/(node-1)
  
  diag(W)<-0
  
  DDelta= solve(diag(x=1,nrow=node)-rho*W)%*%(X%*%BBeta+Tao)
  ####Level 1 model
  Numpat=30 #Num of observations attributed to each node
  N = node*Numpat
  B = matrix(0,nrow = N,ncol = node) #attributing matrix
  for (i in 1:node){
    B[((i-1)*Numpat+1):(i*Numpat),i]=1
  }
  ssigma=1 #true value of sigma
  epi=matrix(rnorm(N, mean=0, sd=sqrt(ssigma)),ncol=1)
  Z=cbind(1,matrix(rnorm(N*3),ncol=3))
  ZZeta=matrix(c(0.5,1,1.5,2),nrow=4) #true value of zeta
  
  y=Z%*%ZZeta+B%*%DDelta+epi #y
  
  # ordered eigenvalues of W
  Lamda = c(sort(eigen(W)$values, decreasing = T))
  
  #####
  #MCMC
  num=20000 # number of iterations
  Rho = rep(0, num)
  Sigma = rep(0, num)
  omega = rep(0, num)
  zetamt = matrix(0, 4 , num)
  betamt = matrix(0, 3, num)
  detamt = matrix(0, node, num)
  #initial values
  Rho[1]=rho
  
  Sigma[1]=ssigma
  omega[1]=oomega
  zetamt[,1]=ZZeta
  betamt[,1]=BBeta
  detamt[,1]=DDelta
  ####
  
  count = 0
  
  i = 2
  while(i<=num){
    
    Rho[i]=rtruncnorm(1,a=1/min(Lamda), b=1/max(Lamda), mean=(t(detamt[,i-1])%*%t(W)%*%(detamt[,i-1]-X%*%betamt[,i-1]))/(omega[i-1]*sum(Lamda^2)+t(detamt[,i-1])%*%t(W)%*%W%*%detamt[,i-1]),
                      sd=sqrt(omega[i-1]/(omega[i-1]*sum(Lamda^2)+t(detamt[,i-1])%*%t(W)%*%W%*%detamt[,i-1]))
    )
    
    proportion = abs(det(diag(1,node)-Rho[i]*W)/det(diag(1,node)-Rho[i-1]*W))*exp(
      -1/(2*omega[i-1])*t((diag(1,node)-Rho[i]*W)%*%detamt[,i-1]-X%*%betamt[,i-1])%*%((diag(1,node)-Rho[i]*W)%*%detamt[,i-1]-X%*%betamt[,i-1])+
        1/(2*omega[i-1])*t((diag(1,node)-Rho[i-1]*W)%*%detamt[,i-1]-X%*%betamt[,i-1])%*%((diag(1,node)-Rho[i-1]*W)%*%detamt[,i-1]-X%*%betamt[,i-1])
      -1/2*((Rho[i-1]-(t(detamt[,i-1])%*%t(W)%*%(detamt[,i-1]-X%*%betamt[,i-1]))/(omega[i-1]*sum(Lamda^2)+t(detamt[,i-1])%*%t(W)%*%W%*%detamt[,i-1]))^2/(omega[i-1]/(omega[i-1]*sum(Lamda^2)+t(detamt[,i-1])%*%t(W)%*%W%*%detamt[,i-1])))+
        1/2*((Rho[i]-(t(detamt[,i-1])%*%t(W)%*%(detamt[,i-1]-X%*%betamt[,i-1]))/(omega[i-1]*sum(Lamda^2)+t(detamt[,i-1])%*%t(W)%*%W%*%detamt[,i-1]))^2/(omega[i-1]/(omega[i-1]*sum(Lamda^2)+t(detamt[,i-1])%*%t(W)%*%W%*%detamt[,i-1])))
      -1/2*((Rho[i]-0.36)/0.7)^2+1/2*((Rho[i-1]-0.36)/0.7)^2
    )
    
    proportion = min(1, Re(proportion))
    
    if (runif(1) > proportion){
      Rho[i]=Rho[i-1]
      count = count+1
    }
    
    Sigma[i] = rinvgamma(n=1, shape = (N-1)/2, rate = t(y-Z%*%zetamt[,i-1]-B%*%detamt[,i-1])%*%(y-Z%*%zetamt[,i-1]-B%*%detamt[,i-1])/2
    )
    
    omega[i] = rinvgamma(n=1, shape = (node-1)/2, rate = t((diag(1,node)-Rho[i]*W)%*%detamt[,i-1]-X%*%betamt[,i-1])%*%((diag(1,node)-Rho[i]*W)%*%detamt[,i-1]-X%*%betamt[,i-1])/2
    )
    
    betamt[,i] = rmvnorm(1, mean = solve(t(X)%*%X)%*%t(X)%*%(diag(1,node)-Rho[i]*W)%*%detamt[,i-1],
                         sigma = solve(t(X)%*%X)*omega[i])
    
    zetamt[,i] = rmvnorm(1, mean = solve(t(Z)%*%Z)%*%t(Z)%*%(y-B%*%detamt[,i-1]), sigma = solve(t(Z)%*%Z)*Sigma[i])
    
    
    Zz = y-Z%*%zetamt[,i]
    H = (sqrt(Sigma[i])/sqrt(omega[i]))*(diag(1,node)-Rho[i]*W)
    DeT = (sqrt(Sigma[i])/sqrt(omega[i]))*X%*%betamt[,i]
    D = t(H)%*%H+t(B)%*%B
    C = (t(DeT)%*%H+t(Zz)%*%B)
    
    detamt[,i] = rmvnorm(1, mean = solve(D)%*%t(C), sigma=Sigma[i]*solve(D)
    )
    
    
    i = i+1
    
  }
  
  medianbeta_m[,k] = apply(betamt[,-(1:2000)],MARGIN=1, FUN=median)
  medianzeta_m[,k] = apply(zetamt[,-(1:2000)], MARGIN = 1, FUN = median)
  medianSigma_v[k] = median(Sigma[-(1:2000)])
  medianomega_v[k] = median(omega[-(1:2000)])
  medianRho_v[k] = median(Rho[-(1:2000)])
  
  lb=quantile(Rho[-(1:2000)],0.025)
  ub=quantile(Rho[-(1:2000)],0.975)
  lbv = c(lbv,lb)
  ubv = c(ubv,ub)
  
  if(rho<=ub & rho>=lb){
    coverage=coverage+1
  }
  
  Rate<-(num-count-1)/(num-1) 
  
  cat("k = ", k,"\n")
}

apply(medianbeta_m,MARGIN=1, FUN=mean)

apply(medianzeta_m,MARGIN=1, FUN=mean)

mean(medianSigma_v)

mean(medianomega_v)

mean(medianRho_v)

mean(medianRho_v)-rho

sum((medianRho_v-rho)^2)/itl

coverage/itl

mean(lbv)

mean(ubv)

mean(ubv)-mean(lbv)
