output<- function(Me, Tb){
  dM<- matrix(0, nrow = nrow(Me), ncol=ncol(Me))
  SM<- matrix(0, nrow = nrow(Me), ncol=ncol(Me))
  
  for (i in 1:ncol(Me)) {
    dM[,i]<- abs(Me[,i]-Tb[i])
    SM[,i]<- (Me[,i]-Tb[i])^2
     }
  MSE<- mean(colMeans(SM, na.rm = T))
  MAE<- mean(colMeans(dM, na.rm = T))
  print(list(MSE=MSE, MAE=MAE))
}

