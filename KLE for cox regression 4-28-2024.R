# New Insights into Multicollinearity in Cox Proportional Hazard Models: The Kibria-Lukman Estimator and its Application
# This code is associated with the article referenced as ---------------------. 
# Please ensure to cite it if you utilize this code.

library(survival)
library(survminer)
source("MSEandMAE.R")

#number of observations (n)
#number of covariates  (p)
#parameter of baseline function  (tau)
#censoring probability   (CP) 
#correlation value    (rho)
#replication number   (m)

CoxKL<- function(n, p, tau, CP, rho, m){
  #true value of the parameters
  Tbeta<- c(-0.5, 0.8, -0.4, 0.9)
  #Tbeta<- c(-0.5, 0.8, -0.4, 0.9, 0.5, -0.5)
  
  MPLE<- matrix(ncol = p, nrow = m)
  Liu1<- matrix(ncol = p, nrow = m)
  RE1<- matrix(ncol = p, nrow = m)
  RE2<- matrix(ncol = p, nrow = m)
  KLE1<- matrix(ncol = p, nrow = m)
  KLE2<- matrix(ncol = p, nrow = m)
  KLE3<- matrix(ncol = p, nrow = m)
  KLE4<- matrix(ncol = p, nrow = m)
  KLE5<- matrix(ncol = p, nrow = m)
  
  for (u in 1:m) {
    
    #generating the observations of covariates
    Z<- matrix(rnorm(n*(p+1)), nrow=n, ncol=(p+1))
    
    X<- matrix(0, nrow=n, ncol=p)
    for (j in 1:p) {
      X[,j]<- sqrt(1-rho^2)*Z[,j]+rho*Z[,5]
    }
    
    # Generate survival times and censoring indicators
    time <- rexp(n, tau)*exp(X %*% Tbeta)
    status <- sample(c(0, 1), n, replace = TRUE, prob = c(CP, 1 - CP))
    gdata <- data.frame(time, status, X)
    
    #Calculate the maximum partial likelihood estimator
    MDcox<- coxph(Surv(time, status) ~ X, data=gdata)
    MPLE[u,]<- coef(MDcox)
    
    # deriving information matrix 
    InDTD<-vcov(MDcox) 
    DTD<- solve(InDTD)
    
    # eigenvalues and eigenvectors of DTD
    lam<- eigen(DTD)$values
    W<- eigen(DTD)$vector
    I<- diag(p)
    
    #ridge parameters 
    #r1<- 1/(t(MPLE[u,])%*%MPLE[u,])
    r2<- min(1/(MPLE[u,]^2))
    
    #Kibria_lukman parameters
    bhat<- t(W)%*%MPLE[u,]
    KhatJ<- lam/(1+2*lam*(bhat^2))
    r3<- mean(KhatJ)
    r4<- max(KhatJ)
    r5<- sqrt(max(0, min(KhatJ)))
    
    #Liu parameters
    a1<- (bhat^2-1)/((lam+1)^2)
    a2<- (lam*(bhat^2)+1)/(lam*((lam+1)^2))
    dopt<- sum(a1)/ sum(a2)
    
    #estimators
    Liu1[u,]<- solve(DTD+I)%*%(DTD+as.numeric(dopt)*I)%*%MPLE[u,]
    #Liu1[u,]<- solve(DTD+I)%*%(DTD+as.numeric(dopt)*I)%*%MPLE[u,]
    
    
    #RE1[u,]<- solve(DTD+as.numeric(r1)*I)%*%DTD%*%MPLE[u,]
    RE2[u,]<- solve(DTD+as.numeric(r2)*I)%*%DTD%*%MPLE[u,]
    
    #KLE1[u,]<- solve(DTD+as.numeric(r1)*I)%*%(DTD-as.numeric(r1)*I)%*%MPLE[u,]
    KLE2[u,]<- solve(DTD+as.numeric(r2)*I)%*%(DTD-as.numeric(r2)*I)%*%MPLE[u,]
    KLE3[u,]<- solve(DTD+as.numeric(r3)*I)%*%(DTD-as.numeric(r3)*I)%*%MPLE[u,]
    KLE4[u,]<- solve(DTD+as.numeric(r4)*I)%*%(DTD-as.numeric(r4)*I)%*%MPLE[u,]
    KLE5[u,]<- solve(DTD+as.numeric(r5)*I)%*%(DTD-as.numeric(r5)*I)%*%MPLE[u,]
  }
  
  
  TMS1<- output(MPLE, Tbeta)
  TMS2<- output(Liu1, Tbeta)
  #TMS3<- output(RE1, Tbeta)
  TMS4<- output(RE2, Tbeta)
  #TMS5<- output(KLE1, Tbeta)
  TMS6<- output(KLE2, Tbeta)
  TMS7<- output(KLE3, Tbeta)
  TMS8<- output(KLE4, Tbeta)
  TMS9<- output(KLE5, Tbeta)
  #TMS10<- output(Liu2, Tbeta)
  
  
  
  TMSE<- c(TMS1$MSE, TMS2$MSE, TMS4$MSE, TMS6$MSE, TMS9$MSE, TMS7$MSE, TMS8$MSE)
  
  TMAE<- c(TMS1$MAE, TMS2$MAE, TMS4$MAE, TMS9$MAE, TMS6$MAE, TMS7$MAE, TMS8$MAE)
  
  TBias<- c(TMS1$Bias, TMS2$Bias, TMS4$Bias, TMS6$Bias, TMS9$Bias, TMS7$Bias, TMS8$Bias)
  
  OUTFINAL<- matrix(nrow = 3, ncol = 7)
  row.names(OUTFINAL)<- c("MSE", "MAE", "Bias")
  OUTFINAL[1,]<- TMSE
  OUTFINAL[2,]<- TMAE
  OUTFINAL[3,]<- TBias
  
  OUTFINAL
}


CoxKLp4n501<- CoxKL(50, 4, 0.2, 0.3, 0.90, 1000)
CoxKLp4n502<- CoxKL(50, 4, 0.2, 0.3, 0.95, 1000)  
CoxKLp4n503<- CoxKL(50, 4, 0.2, 0.3, 0.99, 1000)

CoxKLp4n1001<- CoxKL(100, 4, 0.2, 0.3, 0.90, 1000)
CoxKLp4n1002<- CoxKL(100, 4, 0.2, 0.3, 0.95, 1000)  
CoxKLp4n1003<- CoxKL(100, 4, 0.2, 0.3, 0.99, 1000)

CoxKLp4n2001<- CoxKL(200, 4, 0.2, 0.3, 0.90, 1000)
CoxKLp4n2002<- CoxKL(200, 4, 0.2, 0.3, 0.95, 1000)  
CoxKLp4n2003<- CoxKL(200, 4, 0.2, 0.3, 0.99, 1000)
