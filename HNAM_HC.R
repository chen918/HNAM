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
medianEtavar_v = c()
medianRho_v = c()
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
  etavar=1 #true value of eta_var
  Tao=matrix(rnorm(node, mean = 0, sd = sqrt(etavar)),ncol = 1)
  AA = 25^2
  kksi=rnorm(1,0,sqrt(AA))
  W=rgraph(node, m=1, tprob=ddensity, mode="graph", diag = F, replace=FALSE,
           tielist=NULL, return.as.edgelist= F )
  
  W[lower.tri(W)][which(W[lower.tri(W)]!=0)]=rgamma(n = length(which(W[lower.tri(W)]!=0)), shape = 0.1, scale = 2000)
  
  W[upper.tri(W)]<-t(W)[upper.tri(W)]
  
  W=t(apply(W, 1, function(x) (x/sum(x))))
  
  W[which(is.na(W))]=1/(node-1)
  
  diag(W)<-0
  
  Eta= solve(diag(x=1,nrow=node)-rho*W)%*%(X%*%BBeta/kksi+Tao)
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
  
  y=Z%*%ZZeta+kksi*B%*%Eta+epi #y
  
  # ordered eigenvalues of W
  Lamda = c(sort(eigen(W)$values, decreasing = T))
  
  #####
  #MCMC
  num=20000 # number of iterations
  Rho = rep(0, num)
  Sigma = rep(0, num)
  Etavar = rep(0, num)
  zetamt = matrix(0, 4 , num)
  betamt = matrix(0, 3, num)
  Etamt = matrix(0, node, num)
  si = rep(0, num)
  
  #initial values
  Rho[1]=rho
  si[1]=kksi
  Sigma[1]=ssigma
  Etavar[1]=etavar
  zetamt[,1]=ZZeta
  betamt[,1]=BBeta
  Etamt[,1]=Eta
  ####
  
  count = 0
  countK = 0
  i = 2
  while(i<=num){
    
    Rho[i]=rtruncnorm(1,a=1/min(Lamda), b=1/max(Lamda), mean=(t(Etamt[,i-1])%*%t(W)%*%(Etamt[,i-1]-X%*%betamt[,i-1]/si[i-1]))/(Etavar[i-1]*sum(Lamda^2)+t(Etamt[,i-1])%*%t(W)%*%W%*%Etamt[,i-1]),
                      sd=sqrt(Etavar[i-1]/(Etavar[i-1]*sum(Lamda^2)+t(Etamt[,i-1])%*%t(W)%*%W%*%Etamt[,i-1]))
                      
    )
    
    proportion = (det(diag(1,node)-Rho[i]*W)/det(diag(1,node)-Rho[i-1]*W))*exp(
      -1/(2*Etavar[i-1])*t((diag(1,node)-Rho[i]*W)%*%Etamt[,i-1]-X%*%betamt[,i-1]/si[i-1])%*%((diag(1,node)-Rho[i]*W)%*%Etamt[,i-1]-X%*%betamt[,i-1]/si[i-1])+
        1/(2*Etavar[i-1])*t((diag(1,node)-Rho[i-1]*W)%*%Etamt[,i-1]-X%*%betamt[,i-1]/si[i-1])%*%((diag(1,node)-Rho[i-1]*W)%*%Etamt[,i-1]-X%*%betamt[,i-1]/si[i-1])
      -1/2*((Rho[i-1]-(t(Etamt[,i-1])%*%t(W)%*%(Etamt[,i-1]-X%*%betamt[,i-1]/si[i-1]))/(Etavar[i-1]*sum(Lamda^2)+t(Etamt[,i-1])%*%t(W)%*%W%*%Etamt[,i-1]))^2/(Etavar[i-1]/(Etavar[i-1]*sum(Lamda^2)+t(Etamt[,i-1])%*%t(W)%*%W%*%Etamt[,i-1])))+
        1/2*((Rho[i]-(t(Etamt[,i-1])%*%t(W)%*%(Etamt[,i-1]-X%*%betamt[,i-1]/si[i-1]))/(Etavar[i-1]*sum(Lamda^2)+t(Etamt[,i-1])%*%t(W)%*%W%*%Etamt[,i-1]))^2/(Etavar[i-1]/(Etavar[i-1]*sum(Lamda^2)+t(Etamt[,i-1])%*%t(W)%*%W%*%Etamt[,i-1])))
    )
    
    proportion = min(1, Re(proportion))
    
    if (runif(1) > proportion){
      Rho[i]=Rho[i-1]
      count = count+1
    }
    
    Sigma[i] = rinvgamma(n=1, shape = (N-1)/2, rate = t(y-Z%*%zetamt[,i-1]-si[i-1]*B%*%Etamt[,i-1])%*%(y-Z%*%zetamt[,i-1]-si[i-1]*B%*%Etamt[,i-1])/2
    )
    
    Etavar[i] = rinvgamma(n=1, shape = (node+1)/2, rate = (t((diag(1,node)-Rho[i]*W)%*%Etamt[,i-1]-X%*%betamt[,i-1]/si[i-1])%*%((diag(1,node)-Rho[i]*W)%*%Etamt[,i-1]-X%*%betamt[,i-1]/si[i-1])+1)/2
    )
    
    betamt[,i] = rmvnorm(1, mean = si[i-1]*solve(t(X)%*%X)%*%t(X)%*%(diag(1,node)-Rho[i]*W)%*%Etamt[,i-1],
                         sigma = si[i-1]^2*solve(t(X)%*%X)*Etavar[i])
    zetamt[,i] = rmvnorm(1, mean = solve(t(Z)%*%Z)%*%t(Z)%*%(y-si[i-1]*B%*%Etamt[,i-1]), sigma = solve(t(Z)%*%Z)*Sigma[i])
    
    
    Zz = y-Z%*%zetamt[,i]
    H = (sqrt(Sigma[i])/sqrt(Etavar[i])/si[i-1])*(diag(1,node)-Rho[i]*W)
    DeT = (sqrt(Sigma[i])/sqrt(Etavar[i])/si[i-1])*X%*%betamt[,i]
    D = t(H)%*%H+t(B)%*%B
    C = (t(DeT)%*%H+t(Zz)%*%B)
    
    Etamt[,i] = rmvnorm(1, mean = solve(D)%*%t(C)/si[i-1], sigma=Sigma[i]*solve(D)/si[i-1]^2
    )
    
    si[i] = rnorm(1, mean = t(Etamt[,i])%*%t(B)%*%(y-Z%*%zetamt[,i])/(t(Etamt[,i])%*%t(B)%*%B%*%Etamt[,i]+Sigma[i]/AA), sd = sqrt(Sigma[i]/(t(Etamt[,i])%*%t(B)%*%B%*%Etamt[,i]+Sigma[i]/AA)))
    
    proportionKsi = exp(
      -1/(2*Etavar[i])*t((diag(1,node)-Rho[i]*W)%*%Etamt[,i]-X%*%betamt[,i]/si[i])%*%((diag(1,node)-Rho[i]*W)%*%Etamt[,i]-X%*%betamt[,i]/si[i])+
        1/(2*Etavar[i])*t((diag(1,node)-Rho[i]*W)%*%Etamt[,i]-X%*%betamt[,i]/si[i-1])%*%((diag(1,node)-Rho[i]*W)%*%Etamt[,i]-X%*%betamt[,i]/si[i-1])
      -1/(2*Sigma[i])*t(y-Z%*%zetamt[,i]-si[i]*B%*%Etamt[,i])%*%(y-Z%*%zetamt[,i]-si[i]*B%*%Etamt[,i])+
        1/(2*Sigma[i])*t(y-Z%*%zetamt[,i]-si[i-1]*B%*%Etamt[,i])%*%(y-Z%*%zetamt[,i]-si[i-1]*B%*%Etamt[,i])
      -1/(2*AA)*si[i]^2+
        1/(2*AA)*si[i-1]^2
      -1/2*((si[i-1]-t(Etamt[,i])%*%t(B)%*%(y-Z%*%zetamt[,i])/(t(Etamt[,i])%*%t(B)%*%B%*%Etamt[,i]+Sigma[i]/AA))^2/(Sigma[i]/(t(Etamt[,i])%*%t(B)%*%B%*%Etamt[,i]+Sigma[i]/AA)))+
        1/2*((si[i]-t(Etamt[,i])%*%t(B)%*%(y-Z%*%zetamt[,i])/(t(Etamt[,i])%*%t(B)%*%B%*%Etamt[,i]+Sigma[i]/AA))^2/(Sigma[i]/(t(Etamt[,i])%*%t(B)%*%B%*%Etamt[,i]+Sigma[i]/AA)))
    )
    
    proportionKsi = min(1, Re(proportionKsi))
    
    if (runif(1) > proportionKsi){
      si[i]=si[i-1]
      countK = countK+1
    }
    
    
    i = i+1
    
  }
  
  medianbeta_m[,k] = apply(betamt[,-(1:2000)],MARGIN=1, FUN=median)
  medianzeta_m[,k] = apply(zetamt[,-(1:2000)], MARGIN = 1, FUN = median)
  medianSigma_v[k] = median(Sigma[-(1:2000)])
  medianEtavar_v[k] = median(Etavar[-(1:2000)])
  medianRho_v[k] = median(Rho[-(1:2000)])
  
  
  lb=quantile(Rho[-(1:2000)],0.025)
  ub=quantile(Rho[-(1:2000)],0.975)
  lbv = c(lbv,lb)
  ubv = c(ubv,ub)
  
  if(rho<=ub & rho>=lb){
    coverage=coverage+1
  }
  
  Rate<-(num-count-1)/(num-1)
  RateK<-(num-countK-1)/(num-1)
  
  print(Rate)
  print(RateK)
  cat("k = ", k,"\n")
}

apply(medianbeta_m,MARGIN=1, FUN=mean)

apply(medianzeta_m,MARGIN=1, FUN=mean)

mean(medianSigma_v)

mean(medianEtavar_v)

mean(medianRho_v)

mean(medianRho_v)-rho

sum((medianRho_v-rho)^2)/itl

coverage/itl

mean(lbv)

mean(ubv)

mean(ubv)-mean(lbv)
