### Error Measure
# R. J. Hyndman and A. B. Koehler 2005
my_forecasting_measure<-function(y.hat, y, err.benchmk=y[,1])
{
  y<-as.matrix(y)
  y.hat<-as.matrix(y.hat)
  err.benchmk<-as.matrix(err.benchmk)
  n<-nrow(y)
  
  #Scale dependent: MAE RMSE MAE MdAE
  error<-y-y.hat
  error<-as.matrix(error)
  MSE<-apply(error^2,2,mean, na.rm=T)
  RMSE<-sqrt(MSE)
  MAE<-apply(abs(error),2,mean, na.rm=T)
  MdAE<-apply(abs(error),2,median, na.rm=T)
  
  #Percentage errors: MAPE MdAPE RMSPE RMdSPE
  per.err<-error/y*100
  MAPE<-apply(abs(per.err),2,mean, na.rm=T)
  MdAPE<-apply(abs(per.err),2,median, na.rm=T)
  RMSPE<-sqrt(apply(per.err^2,2,mean, na.rm=T))
  RMdSPE<-sqrt(apply(per.err^2,2,median, na.rm=T))
  
  #Symmetric
  s.per.err<-200*abs(y-y.hat)/(abs(y)+abs(y.hat))  # Changed by Ying to align with the corrected calculation
  s.per.err<-as.matrix(s.per.err)
  sMAPE<-apply(s.per.err,2,mean, na.rm=T)
  sMdAPE<-apply(s.per.err,2,median, na.rm=T)
  
  #Relative: MRAE MdRAE GMRAE
  rela.err<-error/rep(err.benchmk,ncol(error))
  MRAE<-apply(abs(rela.err),2,mean, na.rm=T)
  MdRAE<-apply(abs(rela.err),2,median, na.rm=T)
  GMRAE<- exp(apply(log(abs(rela.err)),2,mean, na.rm=T))
  
  #Scaled
  y<-as.matrix(y)
  n<-nrow(y)
  MASE<-apply(abs(error),2,sum,na.rm=T)/(n/(n-1)*apply(abs(diff(y)),2,sum,na.rm=T))
  
  error.measure<-rbind(MAE, MAPE, MASE,
                       MSE, RMSE,MdAE,
                       MdAPE, RMSPE, RMdSPE,
                       sMAPE, sMdAPE,
                       MRAE, MdRAE, GMRAE)
  rownames(error.measure)<-c('MAE', 'MAPE', 'MASE',
                             'MSE', 'RMSE','MdAE',
                             'MdAPE', 'RMSPE', 'RMdSPE',
                             'sMAPE', 'sMdAPE',
                             'MRAE', 'MdRAE', 'GMRAE')
  colnames(error.measure)<-1:ncol(y)
  return(error.measure)
}