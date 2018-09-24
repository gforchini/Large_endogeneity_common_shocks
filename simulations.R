# File for the simulation in the paper Consistent Estimation in 'Large Heterogeneous Panels with Multifactor Structure 
# and Endogeneity' By G. Forchini, B. Jiang and B. Peng
# Simulations for the model where all explanatory variables are endogenous


# The following packages are required: mnormt, MASS, compiler,doParallel, foreach,doRNG
#install.packages("mnormt")
#install.packages("foreach")
#install.packages("doRNG")
#install.packages("doParallel")
#install.packages("MASS")
#install.packages("compiler")

# Load the required packages

library(MASS)
library(doParallel)
library(foreach)
library(doRNG)

library(mnormt)
library(compiler)
# Set seed for random number
set.seed(2002)

# set working directory
#setwd("~/GF")


# Define functions used in simulations

################## Factors ################


F_gen<-cmpfun(function(Time,m){
  F<-matrix(0,nrow=m,ncol=Time+50) # F is mx(T+50)
  innovation <-matrix(rnorm(m*(Time+50),0,1),nrow=m,ncol=Time+50)
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
  V<-rep(0,Time+50)# Generate a T+50 vector of independent standard normals
  innovation <-rnorm(Time+50,0,.3)
  for (i in 2:length(V)){
    V[i]<-.2*V[i-1]+innovation[i] # eq (25) V generated using AR process
  }
  return(matrix(V[-c(1:50)],ncol=Time,nrow=1))  # This returns a matrix of dimentions 1xT. First 50 obs discarded
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




gG_gen<-cmpfun(function(m,rho_gG){
  # This generates the factor loadings
  Omega <-.1* matrix(c(1,      0       ,rho_gG,
                       0,      1       ,rho_gG,
                       rho_gG, rho_gG  ,1)     ,nrow=3,ncol=3) # eq (23). This is a 3x3 covariance matrix
  # Notice that rho_gG cannot take any value as the covariance matrix must be positive semidefinite
  #a<-1
  #b<-2
  #c<-9;
  #Omega <- .1*matrix(c(a,      -.1*sqrt(a*b)       ,-rho_gG*sqrt(a*c),
  #                  -.1*sqrt(a*b),      b       ,rho_gG*sqrt(b*c),
  #                  -rho_gG*sqrt(a*c), rho_gG*sqrt(b*c)  ,c)     ,nrow=3,ncol=3)
  
  #m_gG<-runif(3*m,min=-2,max=2) # Generate 3m vector of uniform rv between -2  and 2
  gG<-rmnorm(m, rep(0,3), Omega) # generate  mx3 matrix of multivariate normal with
  # mean 0 and covariance Omega
  return(gG)  
} )



######################### random observation simulation #####################
i_gen<-cmpfun(function(Time,cc,rho_e1e2,rho_gG,m,D){ #D is an input as it does not depend on i
  F<-F_gen(Time,m) # Generate factors
  V_2i<-V2_gen(Time) # Generate V2
  Epsilon_i<-Epsilon_gen(Time,rho_e1e2) # Generate Epsilon
  Pi_2i<-cc+rnorm(1,0 ,sqrt(2)) # eq (29) # Generate Pi_2i
  alpha_2i<-runif(1, min = -1, max = 2) # Generate coefficient alpha_i
  lambda_i<-runif(1, min = 1, max = 2) # Generate coefficient lambda_i
  A_2i<-runif(1, min = -1, max = 2) # generate A_2
  gG_i<-gG_gen(m,rho_gG) # Generated factor loadings
  G2_i<-matrix(gG_i[,3],nrow=m,ncol=1) #G2
  g<-gG_i[,c(1,2)] # g
  beta_i<- rmnorm(1,1 , sqrt(2)) # betai~N(1,.2)
  ue2_i<-t(g)%*%F+Epsilon_i # eq (26)
  X_2i<-A_2i* D+t(G2_i)%*%F+V_2i  # eq (24)
  y_2i<-alpha_2i*D+Pi_2i*X_2i+ue2_i[2,] # eq (28)
  y_1i<-lambda_i*D+beta_i*y_2i+ue2_i[1,] # eq (30)
  return(cbind(t(y_1i),t(y_2i),t(X_2i))) # This return a  Tx3 matrix [y1,y2,X2] for an individual i
})



###################### sample generation and calculation of the statistics of interest ###################
simulations<-cmpfun(function(N,Time,cc,rho_e1e2,rho_gG){
  m<-2 # Set number of factors m
  D<-D_gen(Time) # Generate D for all individual (it only depends on time)
  Sample_simulated<-replicate(n =N, i_gen(Time,cc,rho_e1e2,rho_gG,m,D),simplify=FALSE) # sample generated for N individuals
  
  #Construct individual estimators 
  
  #Pesaran CCEMG and CCEP estimators  
  #Cross sectional means
  y1_bar<-rowMeans(sapply(Sample_simulated,function(x) x[,1]))
  y2_bar<-rowMeans(sapply(Sample_simulated,function(x) x[,2]))
  X2_bar<-rowMeans(sapply(Sample_simulated,function(x) x[,3]))
  # define the matrix of means
  X<-cbind(rep(0,T),y1_bar,y2_bar,t(D))
  # define the projection matrix on the means
  M<-diag(Time)-X%*%ginv(t(X)%*%X)%*%t(X)
  
  # Project M (y1,y2,X2)  - notice there is no x1 in our model
  M_y1y2<-lapply(Sample_simulated, function(x){ M%*%x[,c(1,2)]})
  #CCEMG and CCEP estimators
  MG<-mean(sapply(M_y1y2,function(x){(t(x[,2])%*%x[,1])/(t(x[,2])%*%x[,2])}))
  P<-mean(sapply(M_y1y2,function(x){(t(x[,2])%*%x[,1])}))/mean(sapply(M_y1y2,function(x){t(x[,2])%*%x[,2]}))
  
  #HL estimators
  X_HL<-cbind(rep(1,T),X2_bar,t(D)) # Argument of P matrix
  P_HL<-X_HL%*%solve(t(X_HL)%*%X_HL)%*%t(X_HL) # P matrix
  # create new list in which y2 is replaced by P_HL y2
  M_y1P_HLy2X2<-lapply(Sample_simulated, function(x){ cbind(x[,1],P_HL%*%x[,2])})
  # apply ccemg and ccep estimator to this new list
  HL_MG<-mean(sapply(M_y1P_HLy2X2,function(x){(t(x[,2])%*%x[,1])/(t(x[,2])%*%x[,2])}))
  HL_P<-mean(sapply(M_y1P_HLy2X2,function(x){(t(x[,2])%*%x[,1])}))/mean(sapply(M_y1P_HLy2X2,function(x){t(x[,2])%*%x[,2]}))
  
  
  # IV estimators
  H_1<-cbind(rep(1,T),y1_bar,X2_bar,t(D)) #
  H_2<-cbind(rep(1,T),y2_bar,X2_bar,t(D)) #
  M_1<-diag(Time)-H_1%*%solve(t(H_1)%*%H_1)%*%t(H_1)
  M_2<-diag(Time)-H_2%*%solve(t(H_2)%*%H_2)%*%t(H_2)
  
  M1y1_M1X2<-lapply(Sample_simulated, function(x){ cbind(M_1%*%x[,1],M_1%*%x[,3])})
  M2y2_M2X2<-lapply(Sample_simulated, function(x){ cbind(M_2%*%x[,2],M_2%*%x[,3])})
  
  # Estimation of reduced form parameters MG
  pi_1_MG<-mean(sapply(M1y1_M1X2,function(x){(t(x[,2])%*%x[,1])/(t(x[,2])%*%x[,2])}))
  pi_2_MG<-mean(sapply(M2y2_M2X2,function(x){(t(x[,2])%*%x[,1])/(t(x[,2])%*%x[,2])}))
  
  
  # Estimation of reduced form parameters P
  pi_1_P<-mean(sapply(M1y1_M1X2,function(x){(t(x[,2])%*%x[,1])}))/mean(sapply(M1y1_M1X2,function(x){t(x[,2])%*%x[,2]}))
  pi_2_P<-mean(sapply(M2y2_M2X2,function(x){(t(x[,2])%*%x[,1])}))/mean(sapply(M2y2_M2X2,function(x){t(x[,2])%*%x[,2]}))
  
  #IV_MG
  IV_MG<-pi_1_MG/pi_2_MG
  #IVP 
  IV_P<-pi_1_P/pi_2_P
  
  # IV estimators when the RF parameters are estimated jointly
  H<-cbind(rep(1,T),y1_bar,y2_bar,X2_bar,t(D)) #
  MM<-diag(Time)-H%*%solve(t(H)%*%H)%*%t(H)
  
  MMy1_MMy2_MMX2<-lapply(Sample_simulated, function(x){MM%*%x})
  
  # Estimation of reduced form parameters MG
  pi_MGM<-rowMeans(sapply(MMy1_MMy2_MMX2,function(x){(t(x[,3])%*%x[,c(1,2)])/c(t(x[,3])%*%x[,3])}))
  
  pi_PM<-rowMeans(sapply(MMy1_MMy2_MMX2,function(x){(t(x[,3])%*%x[,c(1,2)])}))/mean(sapply(MMy1_MMy2_MMX2,function(x){t(x[,3])%*%x[,3]}))
  
  
  #IV_MG
  IV_MGM<-pi_MGM[1]/pi_MGM[2]
  #IVP 
  IV_PM<-pi_PM[1]/pi_PM[2]
  return(c(MG,P,HL_MG, HL_P,IV_MG,IV_P,IV_MGM,IV_PM))
  #return(list(MG,P,HL_MG, HL_P,IV_MG,IV_P,IV_MGM,IV_PM))
})



################### Replication over T and N ##############


TN<-cmpfun(function(cc,rho_e1e2,rho_gG,N,Time,Nrep){
  m<-2
  D<-D_gen(Time)
  col1<-c('T','N','WI','rho_e1e2','rho_gG','bias_MG','bias_P','bias_HL_MG','bias_HL_P','bias_IV_MG','bias_IV_P','bias_IV_MGJ','bias_IV_PJ',
          'MSE_MG','MSE_P','MSE_HL_MG','MSE_HL_P','MSE_IV_MG','MSE_IV_P','MSE_IV_MGJ','MSE_IV_PJ',
          'medianbias_MG','medianbias_P','medianbias_HL_MG','medianbias_HL_P','medianbias_IV_MG','medianbias_IV_P','medianbias_IV_MGJ',
          'medianbias_IV_PJ',
          ' MEA_MG',' MEA_P',' MEA_HL_MG',' MEA_HL_P',' MEA_IV_MG',' MEA_IV_P',
          ' MEA_IV_MGJ',' MEA_IV_PJ')
  res<-data.frame(Name=col1)
  results<-replicate(Nrep,simulations(N,Time,cc,rho_e1e2,rho_gG),simplify = TRUE)
  
  results<-results-1
  # use median or mean below depending on what is sought
  bias_MG<-mean(results[1,])
  bias_P<-mean(results[2,])
  bias_HL_MG<-mean(results[3,])
  bias_HL_P<-mean(results[4,])
  bias_IV_MG<-mean(results[5,])
  bias_IV_P<-mean(results[6,])
  bias_IV_MGM<-mean(results[7,])
  bias_IV_PM<-mean(results[8,])
  MSE_MG<-mean((results[1,])^2)
  MSE_P<-mean((results[2,])^2)
  MSE_HL_MG<-mean((results[3,])^2)
  MSE_HL_P<-mean((results[4,])^2)
  MSE_IV_MG<-mean((results[5,])^2)
  MSE_IV_P<-mean((results[6,])^2)
  MSE_IV_MGM<-mean((results[7,])^2)
  MSE_IV_PM<-mean((results[8,])^2)
  
  medianbias_MG<-median(results[1,])
  medianbias_P<-median(results[2,])
  medianbias_HL_MG<-median(results[3,])
  medianbias_HL_P<-median(results[4,])
  medianbias_IV_MG<-median(results[5,])
  medianbias_IV_P<-median(results[6,])
  medianbias_IV_MGM<-median(results[7,])
  medianbias_IV_PM<-median(results[8,])
  MEA_MG<-median(abs(results[1,]))
  MEA_P<-median(abs(results[2,]))
  MEA_HL_MG<-median(abs(results[3,]))
  MEA_HL_P<-median(abs(results[4,]))
  MEA_IV_MG<-median(abs(results[5,]))
  MEA_IV_P<-median(abs(results[6,]))
  MEA_IV_MGM<-median(abs(results[7,]))
  MEA_IV_PM<-median(abs(results[8,]))
  
  col2<-c(Time,N,cc,rho_e1e2,rho_gG,bias_MG,bias_P,bias_HL_MG,bias_HL_P,bias_IV_MG,bias_IV_P,bias_IV_MGM,
          bias_IV_PM,
          MSE_MG,MSE_P,MSE_HL_MG,MSE_HL_P,MSE_IV_MG,MSE_IV_P,MSE_IV_MGM,MSE_IV_PM,
          medianbias_MG,medianbias_P,medianbias_HL_MG,medianbias_HL_P,medianbias_IV_MG,medianbias_IV_P,
          medianbias_IV_MGM,medianbias_IV_PM,
          MEA_MG,MEA_P,MEA_HL_MG,MEA_HL_P,MEA_IV_MG,MEA_IV_P,
          MEA_IV_MGM,MEA_IV_PM)
  res<-cbind(res,col2)
  return(res)  
})

Time<-100
N<-100
Nrep<-10
rho_e1e2<-.2
rho_gG<-.1

cc<-3.
TN(cc,rho_e1e2,rho_gG,N,Time,Nrep)
  
  
# matrix containing the combination of parameters considered
#x<- expand.grid(c(0.2,0.6), c(25,50,75,100),c(25,50,75,100))
x<- expand.grid(c(0.2,0.6), c(25,50,75,100),c(25,50,75,100),c(1.,5.))

# Detect the number of cores and make them available for parallel calculations
cores.det<-detectCores()
cores.det
cl <- makeCluster(cores.det)
registerDoParallel(cl)


# Now the actual simulation is done in parallel

registerDoRNG(2005)



#tab<-foreach(i = 1:32, .combine = cbind, .packages=c('mnormt','MASS','compiler'),
#             .export=c('TN','simulations','F_gen','i_gen','D_gen','x','cc',
#                       'V2_gen','Epsilon_gen','gG_gen')) %dopar% TN(cc,x[i,1],.3 ,x[i,2] , x[i,3],100)

tab<-foreach(i = 1:64, .combine = cbind, .packages=c('mnormt','MASS','compiler')) %dopar% TN(x[i,4],x[i,1],.3 ,x[i,2] , x[i,3],10000)

tab1<-tab[,seq(1,64,1)] 
t<-tab1[,c(1,seq(2,64,2))] # table containing all results
rownames(t)<-t[,1]
A<-as.data.frame(t(t[,-1]))
A<-A[order(A$T),]

# Extract the exact results in two data frames
A02.1<-A[A$rho_e1e2==0.2,]
A06.1<-A[A$rho_e1e2==0.6,]

tab1<-tab[,seq(65,64*2,1)] 
t<-tab1[,c(1,seq(2,64,2))] # table containing all results
rownames(t)<-t[,1]
A<-as.data.frame(t(t[,-1]))
A<-A[order(A$T),]

# Extract the exact results in two data frames
A02.2<-A[A$rho_e1e2==0.2,]
A06.2<-A[A$rho_e1e2==0.6,]

#tab1<-tab[,seq(129,64*3,1)] 
#t<-tab1[,c(1,seq(2,64,2))] # table containing all results
#rownames(t)<-t[,1]
#A<-as.data.frame(t(t[,-1]))
#A<-A[order(A$T),]

# Extract the exact results in two data frames
#A02.3<-A[A$rho_e1e2==0.2,]
#A06.3<-A[A$rho_e1e2==0.6,]

# Save the results
#write.csv(tab,"results_simulation5.csv")
write.csv(A02.1,"Results02-12.csv")
write.csv(A06.1,"Results06-12.csv")
write.csv(A02.2,"Results02-22.csv")
write.csv(A06.2,"Results06-22.csv")
#write.csv(A02.3,"Results02-32.csv")
#write.csv(A06.3,"Results06-32.csv")
