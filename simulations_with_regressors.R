# File for the simulation in the paper Consistent Estimation in 'Large Heterogeneous Panels with Multifactor Structure 
# and Endogeneity' By G. Forchini, B. Jiang and B. Peng
# Simulations for the model where one of the explanatory variables is endogenous
# 
# 
# The following packages are required: mnormt, MASS, compiler,doParallel, foreach,doRNG
#install.packages("mnormt")
#install.packages("foreach")
#install.packages("doRNG")
#install.packages("doParallel")
#install.packages("MASS")
#install.packages("compiler")
#install.packages("plyr")

# Load the required packages

library(MASS)
library(doParallel)
library(foreach)
library(doRNG)

library(plyr)
library(mnormt)
library(compiler)
# Set seed for random numbe
set.seed(2002)

# set working directory
#setwd("~/GF")
#
#


# Define functions used in simulations
################## Factors ################
F_gen<-cmpfun(function(Time,m){
  F<-matrix(0,nrow=m,ncol=Time+50) # F is mx(T+50)
  innovation <-matrix(rnorm(m*(Time+50),0,1.),nrow=m,ncol=Time+50)
  for (i in 2:ncol(F)){
    F[,i]<-.8*F[,i-1]+innovation[,i] # Autoregression introduced
  }
  return(F[,-c(1:50)]) # This returns an matrix m x T. The first 50 observations are ignored
  
})



################ D #######################
D_gen<-cmpfun(function(Time){
  D<-rep(0,Time+50) # Generate vector of Time+50 
  innovation <-rnorm(Time+50,0,1)
  for (i in 2:length(D)){
    D[i]<-0.5*D[i-1]+innovation[i] # D is generated as an AR model. The autoregressive parameter could be changed
  }
  return(matrix(D[-c(1:50)],ncol=Time,nrow=1)) # This return a vector of dimension 1XT (first 50 obs discarded)
})


################## V #######################
V2_gen<-cmpfun(function(Time){
  V<-matrix(0,ncol=Time+50,nrow=2)# Generate a T+50 vector of 0
  innovation <-matrix(rnorm(2*(Time+50),0,.3),nrow=2, ncol=Time+50)
  
  for (i in 2:(Time+50)){
    V[,i]<-.2*V[,i-1]+innovation[,i] # eq (25) V generated using AR process
  }
  return(matrix(V[,-c(1:50)],ncol=Time,nrow=2))  # This returns a matrix of dimentions 1xT. First 50 obs discarded
} )



######################### Epsilon ####################
Epsilon_gen<-cmpfun(function(Time,rho_e1e2){
  a<-1.5 # Variance of first component
  b<-1. # Variance of second component
  Omega1<-matrix(c(a,rho_e1e2*sqrt(a*b),rho_e1e2*sqrt(a*b),b),nrow=2,ncol=2) # Covariance matrix
  E<-t(rmnorm(Time,rep(0,2) , Omega1)) # Generate Time pairs of bivariate normals with mean zero and cov=Omega1
  return(E)  # This returns a matrix of dimensions 2xT
} )



########## The factor loadings  ##############

gG_gen<-cmpfun(function(m,rho_gG1,rho_gG2,rho_G1G2){
  # This generates the factor loadings
  Omega <- .1*matrix(c(1,0 ,rho_gG1,rho_gG2,
                       0 ,1,rho_gG1,rho_gG2,
                       rho_gG1,rho_gG1,1,rho_G1G2,
                       rho_gG2,rho_gG2,rho_G1G2,1),nrow=4,ncol=4) # eq (23)
  
  gG<-rmnorm(m, rep(0,4), Omega)
  return(gG)   # This returns a matrix of dimensions m x 4
} )



######################### random observation simulation #####################
i_gen<-cmpfun(function(Time,cc,rho_e1e2,rho_gG1,rho_gG2,rho_G1G2,m,D){
  F<-F_gen(Time,m)
  V_2i<-V2_gen(Time)
  Epsilon_i<-Epsilon_gen(Time,rho_e1e2)
  Pi_2i<-cc+rnorm(1,0 ,sqrt(5)) # eq (29)
  Pi_1i<-rnorm(1,1 ,sqrt(.5)) # eq (29)
  alpha_2i<-runif(1, min = -1, max = 2)
  lambda_i<-runif(1, min = 1, max = 2)
  A_2i<-runif(1, min = -1, max = 2)
  A_1i<-runif(1, min = -1, max = 2)
  gG_i<-gG_gen(m,rho_gG1,rho_gG2,rho_G1G2 )
  G2_i<-matrix(gG_i[,4],nrow=m,ncol=1)
  G1_i<-matrix(gG_i[,3],nrow=m,ncol=1)
  g<-gG_i[,c(1,2)]
  beta1_i<-rmnorm(1,0 , sqrt(.2))
  beta2_i<-1+rmnorm(1,0 , sqrt(2))
  ue2_i<-t(g)%*%F+Epsilon_i # eq (26)
  X_1i<-A_1i* D+t(G1_i)%*%F+V_2i[1,]
  X_2i<-A_2i* D+t(G2_i)%*%F+V_2i[2,]  # eq (24)
  y_2i<-alpha_2i*D+Pi_1i*X_1i+Pi_2i*X_2i+ue2_i[2,] # eq (28)
  y_1i<-lambda_i*D+beta1_i*X_1i+beta2_i*y_2i+ue2_i[1,] # eq (30)
  #return(list(t(y_1i),t(y_2i),t(X_1i),t(X_2i)))
  return(cbind(t(y_1i),t(y_2i),t(X_1i),t(X_2i)))
})







###################### sample generation and calculation of the statistics of interest ###################
simulations<-cmpfun(function(N,Time,cc,rho_e1e2,rho_gG1,rho_gG2,rho_G1G2){
  m<-2 # Set number of factors m
  D<-D_gen(Time)
  Sample_simulated<-replicate(n =N, i_gen(Time,cc,rho_e1e2,rho_gG1,rho_gG2,rho_G1G2,m,D),simplify=FALSE)
  
  #Construct individual estimators 
  
  #Pesaran CCEMG and CCEP estimators  
  
  #Cross sectional means
  y1_bar<-rowMeans(sapply(Sample_simulated,function(x) x[,1]))
  y2_bar<-rowMeans(sapply(Sample_simulated,function(x) x[,2]))
  X1_bar<-rowMeans(sapply(Sample_simulated,function(x) x[,3]))
  X2_bar<-rowMeans(sapply(Sample_simulated,function(x) x[,4]))
  # define the matrix of means
  X<-cbind(rep(0,T),y1_bar,y2_bar,X1_bar,t(D))
  
  # define the projection matrix on the means
  M<-diag(Time)-X%*%ginv(t(X)%*%X)%*%t(X)
  
  # Project M (y1,y2,X1)  
  M_y1y2X1<-lapply(Sample_simulated, function(x){ M%*%x[,c(1,2,3)]})
  
  
  #CCEMG and CCEP estimators
  MG<-rowMeans(sapply(M_y1y2X1,function(x){solve(t(x[,c(2,3)])%*%x[,c(2,3)])%*%t(x[,c(2,3)])%*%x[,1]}))[1]
  P<-(solve(matrix(rowMeans(sapply(M_y1y2X1,function(x){t(x[,c(2,3)])%*%x[,c(2,3)]})),nrow=2,ncol = 2))%*%
        rowMeans(sapply(M_y1y2X1,function(x){(t(x[,c(2,3)])%*%x[,1])})))[1]
  
  # IV estimators
  H_1<-cbind(rep(1,T),y1_bar,X1_bar,X2_bar,t(D)) #
  H_2<-cbind(rep(1,T),y2_bar,X1_bar,X2_bar,t(D)) #
  
  M_1<-diag(Time)-H_1%*%solve(t(H_1)%*%H_1)%*%t(H_1) #M_{H_1}
  M_2<-diag(Time)-H_2%*%solve(t(H_2)%*%H_2)%*%t(H_2) #M_{H_2}
  
  M1_y1X1X2<-lapply(Sample_simulated, function(x){ M_1%*%x[,c(1,3,4)]})
  pi_1_MG<-rowMeans(sapply(M1_y1X1X2,function(x){solve(t(x[,c(2,3)])%*%x[,c(2,3)])%*%t(x[,c(2,3)])%*%x[,1]}))
  pi_1_P<-(solve(matrix(rowMeans(sapply(M1_y1X1X2,function(x){t(x[,c(2,3)])%*%x[,c(2,3)]})),nrow=2,ncol = 2))%*%
             rowMeans(sapply(M1_y1X1X2,function(x){(t(x[,c(2,3)])%*%x[,1])})))
  
  
  M2_y2X1X2<-lapply(Sample_simulated, function(x){ M_2%*%x[,c(2,3,4)]})
  pi_2_MG<-rowMeans(sapply(M2_y2X1X2,function(x){solve(t(x[,c(2,3)])%*%x[,c(2,3)])%*%t(x[,c(2,3)])%*%x[,1]}))
  pi_2_P<-(solve(matrix(rowMeans(sapply(M2_y2X1X2,function(x){t(x[,c(2,3)])%*%x[,c(2,3)]})),nrow=2,ncol = 2))%*%
             rowMeans(sapply(M2_y2X1X2,function(x){(t(x[,c(2,3)])%*%x[,1])})))
  
  
  IV_MG<-pi_1_MG[2]/pi_2_MG[2]
  IV_P<-pi_1_P[2]/pi_2_P[2]
  
  H<-cbind(rep(1,T),y1_bar,y2_bar,X1_bar,X2_bar,t(D))
  M<-diag(Time)-H%*%solve(t(H)%*%H)%*%t(H)
  M_y1X1X2<-lapply(Sample_simulated, function(x){ M%*%x[,c(1,3,4)]})
  M_y2X1X2<-lapply(Sample_simulated, function(x){ M%*%x[,c(2,3,4)]})
  
  pi_1_MGM<-rowMeans(sapply(M_y1X1X2,function(x){solve(t(x[,c(2,3)])%*%x[,c(2,3)])%*%t(x[,c(2,3)])%*%x[,1]}))
  pi_1_PM<-(solve(matrix(rowMeans(sapply(M_y1X1X2,function(x){t(x[,c(2,3)])%*%x[,c(2,3)]})),nrow=2,ncol = 2))%*%
              rowMeans(sapply(M_y1X1X2,function(x){(t(x[,c(2,3)])%*%x[,1])})))
  
  pi_2_MGM<-rowMeans(sapply(M_y2X1X2,function(x){solve(t(x[,c(2,3)])%*%x[,c(2,3)])%*%t(x[,c(2,3)])%*%x[,1]}))
  pi_2_PM<-(solve(matrix(rowMeans(sapply(M_y2X1X2,function(x){t(x[,c(2,3)])%*%x[,c(2,3)]})),nrow=2,ncol = 2))%*%
              rowMeans(sapply(M_y2X1X2,function(x){(t(x[,c(2,3)])%*%x[,1])})))
  
  
  IV_MGM<-pi_1_MGM[2]/pi_2_MGM[2]
  IV_PM<-pi_1_PM[2]/pi_2_PM[2]
  
  return(c(MG,P,IV_MG,IV_P,IV_MGM,IV_PM))
})



################### Replication over T and N ##############


TN<-cmpfun(function(rho_e1e2,rho_gG1,rho_gG2,rho_G1G2,N,Time,cc,Nrep){
  
  col1<-c('T','N','WI',
          'rho_e1e2','rho_gG1','rho_gG2','rho_G1G2',
          'bias_MG','bias_P','bias_IV_MG','bias_IV_P','bias_IV_MGJ','bias_IV_PJ',
          'MSE_MG','MSE_P','MSE_IV_MG','MSE_IV_P','MSE_IV_MG.M','MSE_IV_P.M',
          'medianbias_MG','medianbias_P','medianbias_IV_MG','medianbias_IV_P','medianbias_IV_MGJ',
          'medianbias_IV_PJ',
          ' MAE_MG',' MAE_P',' MAE_IV_MG',' MAE_IV_P',' MAE_IV_MG.M',' MAE_IV_P.M')
  
  
  
  res<-data.frame(Name=col1)
  results<-replicate(Nrep,simulations(N,Time,cc,rho_e1e2,rho_gG1,rho_gG2,rho_G1G2))
  
  results<-results-1
  # use median or mean below depending on what is sought
  bias_MG<-mean(results[1,])
  bias_P<-mean(results[2,])
  bias_IV_MG<-mean(results[3,])
  bias_IV_P<-mean(results[4,])
  bias_IV_MGJ<-mean(results[5,])
  bias_IV_PJ<-mean(results[6,])
  MSE_MG<-mean(results[1,]^2)
  MSE_P<-mean(results[2,]^2)
  MSE_IV_MG<-mean(results[3,]^2)
  MSE_IV_P<-mean(results[4,]^2)
  MSE_IV_MGJ<-mean(results[5,]^2)
  MSE_IV_PJ<-mean(results[6,]^2)
  
  medianbias_MG<-median(results[1,])
  medianbias_P<-median(results[2,])
  medianbias_IV_MG<-median(results[3,])
  medianbias_IV_P<-median(results[4,])
  medianbias_IV_MGJ<-median(results[5,])
  medianbias_IV_PJ<-median(results[6,])
  MAE_MG<-median((abs(results[1,])))
  MAE_P<-median((abs(results[2,])))
  MAE_IV_MG<-median((abs(results[3,])))
  MAE_IV_P<-median((abs(results[4,])))
  MAE_IV_MGJ<-median((abs(results[5,])))
  MAE_IV_PJ<-median((abs(results[6,])))
  
  col2<-c(Time,N,cc,rho_e1e2,rho_gG1,rho_gG2,rho_G1G2,bias_MG,bias_P,bias_IV_MG,bias_IV_P,bias_IV_MGJ,bias_IV_PJ,
          MSE_MG,MSE_P,MSE_IV_MG,MSE_IV_P,MSE_IV_MGJ,MSE_IV_PJ,
          medianbias_MG,medianbias_P,medianbias_IV_MG,medianbias_IV_P,medianbias_IV_MGJ,medianbias_IV_PJ,
          MAE_MG,MAE_P,MAE_IV_MG,MAE_IV_P,MAE_IV_MGJ,MAE_IV_PJ)
  
  res<-cbind(res,col2)
  
  return(res)  
})

Time<-100
N<-100
Nrep<-100
rho_e1e2<-.2
rho_gG1<-.1
rho_gG2<-.1
rho_G1G2<-.3
cc<-3.
TN(rho_e1e2,rho_gG1,rho_gG2,rho_G1G2,N,Time,cc,Nrep)

# matrix containing the combination of parameters considered
x<- expand.grid(c(0.2,0.6), c(25,50,75,100),c(25,50,75,100),c(1.,5.))

# Detect the number of cores and make them available for parallel calculations
cores.det<-detectCores()
cores.det
cl <- makeCluster(cores.det)
registerDoParallel(cl)
# Now the actual simulation is done in parallel

registerDoRNG(2005)

tab<-foreach(i = 1:64, .combine = cbind, .packages=c('mnormt','MASS','compiler')) %dopar%  
  TN(x[i,1],.1,.1,.3,x[i,2] , x[i,3],x[i,4],10000)

tab1<-tab[,seq(1,64,1)] 
t<-tab1[,c(1,seq(2,64,2))] # table containing all results
rownames(t)<-t[,1]
A<-as.data.frame(t(t[,-1]))
A<-A[order(A$T),]

# Extract the exact results in two data frames
A02.1<-A[A$rho_e1e2==0.2,]
A06.1<-A[A$rho_e1e2==0.6,]

tab1<-tab[,seq(65,128,1)] 
t<-tab1[,c(1,seq(2,64,2))] # table containing all results
rownames(t)<-t[,1]
A<-as.data.frame(t(t[,-1]))
A<-A[order(A$T),]

# Extract the exact results in two data frames
A02.2<-A[A$rho_e1e2==0.2,]
A06.2<-A[A$rho_e1e2==0.6,]

#tab1<-tab[,seq(129,192,1)] 
#t<-tab1[,c(1,seq(2,64,2))] # table containing all results
#rownames(t)<-t[,1]
#A<-as.data.frame(t(t[,-1]))
#A<-A[order(A$T),]

# Extract the exact results in two data frames
#A02.3<-A[A$rho_e1e2==0.2,]
#A06.3<-A[A$rho_e1e2==0.6,]

# Save the results
#write.csv(tab,"results_simulation5.csv")
write.csv(A02.1,"ResultsR02-13.csv")
write.csv(A06.1,"ResultsR06-13.csv")
write.csv(A02.2,"ResultsR02-23.csv")
write.csv(A06.2,"ResultsR06-23.csv")
#write.csv(A02.3,"ResultsR02-33.csv")
#write.csv(A06.3,"ResultsR06-33.csv")