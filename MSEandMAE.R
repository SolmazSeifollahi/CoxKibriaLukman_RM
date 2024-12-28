output<- function(Me, Tb){
  dM<- matrix(0, nrow = nrow(Me), ncol=ncol(Me))
  SM<- matrix(0, nrow = nrow(Me), ncol=ncol(Me))
  Bias2<- matrix(0, nrow = nrow(Me), ncol=ncol(Me))
  for (i in 1:ncol(Me)) {
    dM[,i]<- abs(Me[,i]-Tb[i])
    SM[,i]<- (Me[,i]-Tb[i])^2
    Bias2[,i]<- abs((Me[,i]/Tb[i])-1)
  }
  MSE<- mean(colMeans(SM, na.rm = T))
  MAE<- mean(colMeans(dM, na.rm = T))
  B2<- mean(colMeans(Bias2, na.rm = T))
print(list(MSE=MSE, MAE=MAE, Bias=B2))
}

