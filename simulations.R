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
library(mnormt)
library(MASS)
library(compiler)
library(doParallel)
library(foreach)
library(doRNG)

# Set seed for random numbe
set.seed(2002)

# Detet the number of cores and make them available for parallel calculations
cores.det<-detectCores()
cores.det
cl <- makeCluster(cores.det)
registerDoParallel(cl)

# Define functions used in simulations

################## Factors ################
factors_generation<-cmpfun(function(Time,m){
  C <- matrix(rnorm(m^2,0,1),nrow=m,ncol=m)
  F<-t(C)%*%matrix(rnorm(m*(Time+50),0,1),nrow=m,ncol=Time+50)
  for (i in 2:ncol(F)){
    F[,i]<-F[,i]+0.8*F[,i-1] 
    # The autoregresive parameter could be changed
  }
  return(F[,-c(1:50)])
  # This returns an matrix Time x m
})




################ D #######################
D_generation<-cmpfun(function(Time){
  c1<-rmnorm(1,0,1)
  D<-c1*rnorm(Time+50,0,1)
  for (i in 2:length(D)){
    D[i]<-D[i]+0.5*D[i-1] # The autoregressive parameter could be changed
  }
  return(matrix(D[-c(1:50)],ncol=Time,nrow=1))
  # This return a vector of dimension Time
})



################## V #######################
V2_generation<-cmpfun(function(Time){
  V<-t(rmnorm(Time+50,0 , 1))
  for (i in 2:ncol(V)){
    V[,i]<-V[,i]+0.2*V[,i-1] # eq (25)
  }
  return(V[,-c(1:50)])
  # This produces a matrix of dimentions Time x 1
} )



######################### Epsilo ####################
Epsilon_generation<-cmpfun(function(Time,rho_e1e2){
  a<-1.5
  b<-1
  Omega1<-matrix(c(a,rho_e1e2*sqrt(a*b),rho_e1e2*sqrt(a*b),b),nrow=2,ncol=2)
  E<-t(rmnorm(Time,rep(0,2) , Omega1))
  return(E)
  # This returns a matrix of dimensions Time x 2
} )



########## The factor loadings are generated in two different ways ##############

gG_generation<-cmpfun(function(m,rho_gG){
  # This generates the factor loadings
  Omega <- matrix(c(1,0,rho_gG,
                    0,1,rho_gG,
                    rho_gG,rho_gG,1),nrow=3,ncol=3) # eq (23)
  cov_gG<-kronecker(Omega, diag(m), FUN = "*", make.dimnames = FALSE)
  m_gG<-runif(3*m,min=-2,max=2)
  gG<-matrix(rmnorm(9, m_gG, cov_gG),nrow=m,ncol=3)
  return(gG)
  # This returns a matrix of dimensions m x 4
} )



######################### random observation simulation #####################
individual_generation<-cmpfun(function(Time,rho_e1e2,rho_gG,m,D){
  F<-factors_generation(Time,m)
  cc<-1.1 # Strength of the instruments
  V_2i<-V2_generation(Time)
  Epsilon_i<-Epsilon_generation(Time,rho_e1e2)
  Pi_2i<-sqrt(cc)+rnorm(1,0 ,.05) # eq (29)
  alpha_2i<-runif(1, min = -1, max = 2)
  lambda_i<-runif(1, min = 1, max = 2)
  A_2i<-runif(1, min = -1, max = 2)
  gG_i<-gG_generation(m,rho_gG)
  G2_i<-matrix(gG_i[,3],nrow=m,ncol=1)
  g<-gG_i[,c(1,2)]
  beta_i<-1+rmnorm(1,0 , .2)
  ue2_i<-t(g)%*%F+Epsilon_i # eq (26)
  X_2i<-A_2i* D+t(G2_i)%*%F+V_2i  # eq (24)
  y_2i<-alpha_2i*D+Pi_2i*X_2i+ue2_i[2,] # eq (28)
  y_1i<-lambda_i*D+beta_i*y_2i+ue2_i[1,] # eq (30)
  return(list(t(y_1i),t(y_2i),t(X_2i)))
})


###################### sample generation and calculation of the statistics of interest ###################
simulations<-cmpfun(function(N,Time,rho_e1e2,rho_gG){
  m<-2 # Set number of factors m
  D<-D_generation(Time)
  Sample_simulated<-replicate(n =N, individual_generation(Time,rho_e1e2,rho_gG,m,D))
  #Construct individual estimators 

  #Pesaran CCEMG and CCEP estimators  
  #means
  y1_bar<-rowMeans(sapply(Sample_simulated[1,], unlist))
  y2_bar<-rowMeans(sapply(Sample_simulated[2,], unlist))
  X2_bar<-rowMeans(sapply(Sample_simulated[3,], unlist))
  # define the matrix of means
  X<-cbind(y1_bar,y2_bar,X2_bar,t(D))
  # define the projection matrix on the means
  M<-diag(Time)-X%*%ginv(t(X)%*%X)%*%t(X)
  # Project y1 and y2  - notice there is no x1 in our model
  M_times_y1_i<-lapply(Sample_simulated[1,], function(x){ M%*%x})
  M_times_y2_i<-lapply(Sample_simulated[2,], function(x){ M%*%x})
  
  # The CCEMG estimator is calculated
  MG<-mean(mapply(function(x, y) (t(x)%*%y)/(t(x)%*%x), M_times_y2_i, M_times_y1_i))
  
  # The CCEP estimator is calculated
  P<-sum(mapply(function(x, y)(t(x)%*%y), M_times_y2_i, M_times_y1_i))/
    sum(sapply(M_times_y2_i,function(x)(t(x)%*%x)))
  
  #HL estimators
  X_HL<-cbind(rep(1,T),X2_bar,t(D)) # Argument of P matrix
  P_HL<-X_HL%*%solve(t(X_HL)%*%X_HL)%*%t(X_HL) # P matrix
  XX<-cbind(y1_bar,P_HL%*%y2_bar,t(D)) #Argument of the M matrix
  M_HL<-diag(Time)-XX%*%ginv(t(XX)%*%XX)%*%t(XX) # M matrix
  M_HL_times_y1_i<-lapply(Sample_simulated[1,], function(x){ M_HL%*%x})
  M_HL_times_y2_i<-lapply(Sample_simulated[2,], function(x){ M_HL%*%P_HL%*%x})
  HL_MG<-mean(mapply(function(x, y) solve(t(x)%*%x)%*%(t(x)%*%y), M_HL_times_y2_i, M_HL_times_y1_i))
  HL_P<-sum(mapply(function(x, y)(t(x)%*%y), M_HL_times_y2_i, M_HL_times_y1_i))/sum(sapply(M_HL_times_y2_i,function(x)(t(x)%*%x)))
  
  # IV estimators
  H_1<-cbind(rep(1,T),y1_bar,X2_bar,t(D)) #
  H_2<-cbind(rep(1,T),y2_bar,X2_bar,t(D)) #
  M_1<-diag(Time)-H_1%*%solve(t(H_1)%*%H_1)%*%t(H_1)
  M_2<-diag(Time)-H_2%*%solve(t(H_2)%*%H_2)%*%t(H_2)
  M_1_times_y1_i<-lapply(Sample_simulated[1,], function(x){ M_1%*%x})
  M_1_times_X2_i<-lapply(Sample_simulated[3,], function(x){ M_1%*%x})
  M_2_times_X2_i<-lapply(Sample_simulated[3,], function(x){ M_2%*%x})
  M_2_times_y2_i<-lapply(Sample_simulated[2,], function(x){ M_2%*%x})
  
  # Estimation of reduced form parameters MG
  Pi_1_MG<-sum(mapply(function(x, y) solve(t(x)%*%x)%*%(t(x)%*%y), M_1_times_X2_i, M_1_times_y1_i))
  Pi_2_MG<-sum(mapply(function(x, y) solve(t(x)%*%x)%*%(t(x)%*%y), M_2_times_X2_i, M_2_times_y2_i))
 
  #IV_MG
  IV_MG<-(t(Pi_2_MG)%*%Pi_1_MG)/(t(Pi_2_MG)%*%Pi_2_MG)
  
  # Estimation of reduced form parameters P
  A<-rowSums(mapply(function(x, y) x*y, M_1_times_X2_i, M_1_times_y1_i))
  B<-rowSums(sapply(M_1_times_X2_i,function(x) x^2))
  Pi_1_P<-(t(B)%*%A)/(t(B)%*%B)
  
  A<-rowSums(mapply(function(x, y) x*y, M_2_times_X2_i, M_2_times_y2_i))
  B<-matrix(rowSums(sapply(M_2_times_X2_i,function(x)x^2))) 
  Pi_2_P<-(t(B)%*%A)/(t(B)%*%B)
  
  #IVP
  IV_P<-Pi_1_P/Pi_2_P
  
  return(list(MG,P,HL_MG, HL_P,IV_MG,IV_P))
})



################### Replication over T and N ##############


TN<-cmpfun(function(rho_e1e2,rho_gG,N,Time,Nrep){
  
  col1<-c('T','N','rho_e1e2','rho_gG','bias_MG','bias_P','bias_HL_MG','bias_HL_P','bias_IV_MG','bias_IV_P','MSE_MG','MSE_P','MSE_HL_MG','MSE_HL_P','MSE_IV_MG','MSE_IV_P')
  res<-data.frame(Name=col1)
  results<-replicate(Nrep,simulations(N,Time,rho_e1e2,rho_gG))
  # use median or mean below depending on what is sought
  bias_MG<-median(unlist(results[1,]))-1
  bias_P<-median(unlist(results[2,]))-1
  bias_HL_MG<-median(unlist(results[3,]))-1
  bias_HL_P<-median(unlist(results[4,]))-1
  bias_IV_MG<-median(unlist(results[5,]))-1
  bias_IV_P<-median(unlist(results[6,]))-1
  MSE_MG<-median((unlist(results[1,])-1)^2)
  MSE_P<-median((unlist(results[2,])-1)^2)
  MSE_HL_MG<-median((unlist(results[3,])-1)^2)
  MSE_HL_P<-median((unlist(results[4,])-1)^2)
  MSE_IV_MG<-median((unlist(results[5,])-1)^2)
  MSE_IV_P<-median((unlist(results[6,])-1)^2)
  col2<-c(Time,N,rho_e1e2,rho_gG,bias_MG,bias_P,bias_HL_MG,bias_HL_P,bias_IV_MG,bias_IV_P,MSE_MG,MSE_P,MSE_HL_MG,MSE_HL_P,MSE_IV_MG,MSE_IV_P)
  res<-cbind(res,col2)
  return(res)  
})


# matrix containing the combination of parameters considered
x<- expand.grid(c(0.2,0.6), c(25,50,75,100),c(25,50,75,100))


# Now the actual simulation is done in parallel

registerDoRNG(2005)

tab<-foreach(i = 1:32, .combine = cbind, .packages=c('mnormt','MASS','compiler'),
             .export=c('TN','simulations','factors_generation','individual_generation','D_generation','x',
                       'V2_generation','Epsilon_generation','gG_generation')) %dopar% TN(x[i,1],.3 ,x[i,2] , x[i,3],10000)

t<-tab[,c(1,seq(2,64,2))] # table containing all results

rownames(t)<-t[,1]
A<-as.data.frame(t(t[,-1]))
A<-A[order(A$T),]

# Extract the exact results in two data frames
A02<-A[A$rho_e1e2==0.2,]
A06<-A[A$rho_e1e2==0.6,]

# Save the results
write.csv(tab,"results_simulation5.csv")
write.csv(A02,"Results02.csv")
write.csv(A06,"Results06.csv")
