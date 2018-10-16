#install.packages('EMD')
#install.packages('e1071')
#install.packages("hydroGOF")
#install.packages('caret')


library('EMD')
library(ggplot2)
library(grid)
library(dplyr)
library('forecast')
library('urca')
library('astsa')
library('tseries')
library('e1071')
library(hydroGOF)

### function required: my_vector_split,my_scale_ts, my_cv_svm (svm+grid search+roll/direct forecasting) ###
### function required: my_cv_partition, my_forecasting ###

MAX.IMF=6
HORIZON=1
step_ahead=6
tt2<-seq((len_train+1),length(Ottseries), by=1) # for test set time series
tt3<-seq((length(Ottseries)+1),(length(Ottseries)+step_ahead),by=1) # for h-steps forecasting
# extract EMD components
emd.xts<- emd(as.numeric(Ottseries), max.imf=MAX.IMF)
nimf<-emd.xts$nimf
# plot imf and residue
g0 <- df %>%
  select(DateTime, Growth_Pct) %>%
  ggplot() +
  geom_line(aes(x = DateTime, y = Growth_Pct), size=0.3) +
  theme_minimal()+
  theme(axis.title.x = element_blank())
g1 <- df %>%
  select(DateTime, imf1) %>%
  ggplot() +
  geom_line(aes(x = DateTime, y = imf1), size=0.3) +
  theme_minimal()+
  theme(axis.title.x = element_blank())
g2 <- df %>%
  select(DateTime, imf2) %>%
  ggplot() +
  geom_line(aes(x = DateTime, y = imf2), size=0.3) +
  theme_minimal()+
  theme(axis.title.x = element_blank())
g3 <- df %>%
  select(DateTime, imf3) %>%
  ggplot() +
  geom_line(aes(x = DateTime, y = imf3), size=0.3) +
  theme_minimal()+
  theme(axis.title.x = element_blank())
g4 <- df %>%
  select(DateTime, imf4) %>%
  ggplot() +
  geom_line(aes(x = DateTime, y = imf4), size=0.3) +
  theme_minimal()+
  theme(axis.title.x = element_blank())
g5 <- df %>%
  select(DateTime, imf5) %>%
  ggplot() +
  geom_line(aes(x = DateTime, y = imf5), size=0.3) +
  theme_minimal()+
  theme(axis.title.x = element_blank())
g6 <- df %>%
  select(DateTime, imf6) %>%
  ggplot() +
  geom_line(aes(x = DateTime, y = imf6), size=0.3) +
  theme_minimal()+
  theme(axis.title.x = element_blank())
g7 <- df %>%
  select(DateTime, res) %>%
  ggplot() +
  geom_line(aes(x = DateTime, y = res), size=0.3) +
  theme_minimal()+
  theme(axis.title.x = element_blank())

grid.newpage()
grid.draw(rbind(ggplotGrob(g0), ggplotGrob(g1), ggplotGrob(g2), ggplotGrob(g3),
                ggplotGrob(g4), ggplotGrob(g5), ggplotGrob(g6), ggplotGrob(g7), size = "last"))

# Order of the (ARIMA(p,d,q)) models fitted to each intrinsic mode function (IMF) and to the residue.
# The number of lagged values m = p + d was used to construct the input vectors
# for the empirical mode decomposition–support vector regression (EMD–SVR) models.

###############################################################
####################### Build EMD-AR #########################
##############################################################
AR_LAG <-list()
for (i in 1:(nimf+1))
    if(i>nimf) {
      fit<-auto.arima(emd.xts$residue,seasonal=FALSE) 
      p<-length(fit$model$phi)
      d<-length(fit$model$Delta)
      q<-length(fit$model$theta)
      AR_LAG[i]<-p+d
    }else {
      fit<-auto.arima(emd.xts$imf[,i],seasonal=FALSE) 
      p<-length(fit$model$phi)
      d<-length(fit$model$Delta)
      q<-length(fit$model$theta)
      AR_LAG[i]<-p+d
    }

AR_LAG<-as.numeric(unlist(AR_LAG))
#auto.arima(emd.xts$imf[,1],seasonal=FALSE)  #imf1 (1,0,4)  so m=1
#auto.arima(emd.xts$imf[,2],seasonal=FALSE)  #imf1 (2,0,3)  so m=2
#auto.arima(emd.xts$imf[,3],seasonal=FALSE)  #imf3 (4,0,2)  so m=4
#auto.arima(emd.xts$imf[,4],seasonal=FALSE)  #imf4 (4,0,4)  so m=4
#auto.arima(emd.xts$imf[,5],seasonal=FALSE)  #imf5 (3,0,5)  so m=3
#auto.arima(emd.xts$imf[,6],seasonal=FALSE)  #imf6 (0,1,0)  so m=1
#auto.arima(emd.xts$res,seasonal=FALSE)   #res (0,2,5)   so m=2

###############################################################
####################### Train SVR Model #######################
###############################################################
# Model Selection: RBF Kernel: C and rou + CV + Grid Search

# build SVR for each IMF
emd.trn<-emd.tst<-list()
emd.svr.predict<-array(NA,dim=c(nrow(tst$labels),ncol(tst$labels),nimf+1))
emd.svr.parameter<-list()
for (i in 1:(nimf+1)){ # for each IMF, +1 is for residue
  if(i>nimf){ # for residue
    tmp<-my_vector_split(emd.xts$residue,HORIZON=HORIZON,LAG=AR_LAG[i],RATIO=RATIO)
  }else{ # for imf
    tmp<-my_vector_split(emd.xts$imf[,i],HORIZON=HORIZON,LAG=AR_LAG[i],RATIO=RATIO)
  }
  emd.trn$data[[i]]<-tmp$trn$data
  emd.trn$labels[[i]]<-tmp$trn$labels
  emd.tst$data[[i]]<-tmp$tst$data
  emd.tst$labels[[i]]<-tmp$tst$labels
  
  # cv to determine best k and d
  emd.svr.p.model<-my_cv_svm(tmp$trn,tmp$tst,is.ts=TRUE,cv.criteria='MSE')
  # store model and parameter information
  emd.svr.parameter[i]<-emd.svr.p.model$model
  # store predicted value
  emd.svr.predict[,,i]<-emd.svr.p.model$predict
}

# take sum of IMFs
emd.svr.p.predict<-array(NA,dim=c(nrow(tst$labels),HORIZON))
for (h in 1:HORIZON)
  emd.svr.p.predict[,h]<-apply(emd.svr.predict[,h,],1,sum)

# calculate error
emd.svr.p.error<-my_forecasting_measure(emd.svr.p.predict,tst$labels)

# visulization comparison between actual and forecasting for IMF
plot(emd.xts$imf[,1],type='l', main="IMF", xlab='Year',ylab='IMF 1',col='black')
lines(tt2,emd.svr.predict[,,1],col='red',pch=16)
legend(1,95,legend=c("IMF 1", "Forecasting"),
       col=c("black", "red"), lty=1:2,cex=0.8)

plot(emd.xts$imf[,2],type='l', main="IMF", xlab='Year',ylab='IMF 2')
lines(tt2,emd.svr.predict[,,2],col='red',pch=16)

plot(emd.xts$imf[,3],type='l', main="IMF", xlab='Year',ylab='IMF 3')
lines(tt2,emd.svr.predict[,,3],col='red',pch=16)

plot(emd.xts$imf[,4],type='l', main="IMF", xlab='Year',ylab='IMF 4')
lines(tt2,emd.svr.predict[,,4],col='red',pch=16)

plot(emd.xts$imf[,5],type='l', main="IMF", xlab='Year',ylab='IMF 5')
lines(tt2,emd.svr.predict[,,5],col='red',pch=16)

plot(emd.xts$imf[,6],type='l', main="IMF", xlab='Year',ylab='IMF 6')
lines(tt2,emd.svr.predict[,,6],col='red',pch=16)

plot(emd.xts$residue,type='l', main="Residue", xlab='Year',ylab='Residue')
lines(tt2,emd.svr.predict[,,7],col='red',pch=16)

# visulizatio comparison between actual and forecasting
emd.svr.p.predict.ts<-ts(emd.svr.p.predict,frequency=12,start=c(2012,05))
plot.ts(Ottseries,type='l', main="Annaulized Monthly New House Price Index", xlab='Year',ylab='NHPI',axes=F)
lines(emd.svr.p.predict.ts,col='red',pch=16)
axis(1,at=seq(1998,2018,by=2))
axis(2)

###############################################################
####################### Forecasting ###########################
###############################################################
# use roll forecasting

# emd.trn.w<-emd.tst.w<-list()
# emd.svr.predict.whole<-array(NA,dim=c(nrow(tst$labels),ncol(tst$labels),nimf+1))
# emd.svr.parameter.whole<-list()
# for (i in 1:(nimf+1)){ # for each IMF, +1 is for residue
#   if(i>nimf){ # for residue
#     tmp<-my_vector_split(emd.xts$residue,HORIZON=HORIZON,LAG=AR_LAG[i],RATIO=1.0)
#   }else{ # for imf
#     tmp<-my_vector_split(emd.xts$imf[,i],HORIZON=HORIZON,LAG=AR_LAG[i],RATIO=1.0)
#   }
#   emd.trn.w$data[[i]]<-tmp$trn$data
#   emd.trn.w$labels[[i]]<-tmp$trn$labels
#   emd.tst.w$data[[i]]<-tmp$tst$data
#   emd.tst.w$labels[[i]]<-tmp$tst$labels
# }

# forecast  steps ahead
# Need to manually update based on previous AR result
model<-list()
newdata<-list()
forecast_value<-list()
for (i in 1:(nimf+1)){
  newdata[[i]]<-matrix(nrow=step_ahead,ncol=AR_LAG[i])
  if (i>nimf){
    model[[i]]<-svm(emd.trn$labels[[i]] ~ ., data=emd.trn$data[[i]], kernel="radial", cost=1000, gamma=0.5, epsilon=0.0001)
    newdata1<-emd.xts$residue[(length(emd.xts$residue)-AR_LAG[i]+1):length(emd.xts$residue)]
  }else if(i==1){
    model[[i]]<-svm(emd.trn$labels[[i]] ~ ., data=emd.trn$data[[i]], kernel="radial", cost=0.1, gamma=1, epsilon=0.0001)
    newdata1<-emd.xts$imf[,i][(length(emd.xts$imf[,i])-AR_LAG[i]+1):length(emd.xts$imf[,i])]
  }else if(i==2){
    model[[i]]<-svm(emd.trn$labels[[i]] ~ ., data=emd.trn$data[[i]], kernel="radial", cost=1, gamma=0.5, epsilon=0.1)
    newdata1<-emd.xts$imf[,i][(length(emd.xts$imf[,i])-AR_LAG[i]+1):length(emd.xts$imf[,i])]
  }else if(i==3){
    model[[i]]<-svm(emd.trn$labels[[i]] ~ ., data=emd.trn$data[[i]], kernel="radial", cost=100, gamma=0.25, epsilon=0.0001)
    newdata1<-emd.xts$imf[,i][(length(emd.xts$imf[,i])-AR_LAG[i]+1):length(emd.xts$imf[,i])]
  }else if(i==4){
    model[[i]]<-svm(emd.trn$labels[[i]] ~ ., data=emd.trn$data[[i]], kernel="radial", cost=1000, gamma=0.25, epsilon=0.001)
    newdata1<-emd.xts$imf[,i][(length(emd.xts$imf[,i])-AR_LAG[i]+1):length(emd.xts$imf[,i])]
  }else if(i==5){
    model[[i]]<-svm(emd.trn$labels[[i]] ~ ., data=emd.trn$data[[i]], kernel="radial", cost=100, gamma=0.33333, epsilon=0.001)
    newdata1<-emd.xts$imf[,i][(length(emd.xts$imf[,i])-AR_LAG[i]+1):length(emd.xts$imf[,i])]
  }else if(i==6){
    model[[i]]<-svm(emd.trn$labels[[i]] ~ ., data=emd.trn$data[[i]], kernel="radial", cost=1000, gamma=1, epsilon=0.1)
    newdata1<-emd.xts$imf[,i][(length(emd.xts$imf[,i])-AR_LAG[i]+1):length(emd.xts$imf[,i])]
  }
  for (h in 1:step_ahead){
    if (h==1) {
      newdata[[i]][h,]<-rbind(newdata1)
      forecast_value[[i]]<-predict(model[[i]],newdata[[i]])
    }else{
      # minus ? depends on the lag for each model
      newdata[[i]][h,]<-matrix(c(newdata[[i]][h-1,-AR_LAG[i]+1],forecast_value[[i]][h-1]))
      forecast_value[[i]]<-predict(model[[i]],newdata[[i]])
    }
  }
  
  imf_position=data.frame(i=i)
  if (i>nimf) {
    plot(emd.xts$residue,type='l', main="Residue", xlab='Year',axes=T, ylab='Residue')
    lines(tt3,forecast_value[[i]],col='blue',pch=16)
  }else{
    plot(emd.xts$imf[,i],type='l', main=bquote('IMF'~.(imf_position$i)), xlab='Year',axes=T, ylab=bquote('IMF'~.(imf_position$i)))
    lines(tt3,forecast_value[[i]],col='blue',pch=16)
  }
}



forecast.w<-array()
for (h in 1:step_ahead){
  for (i in 1:nimf+1){
    forecast.w[h]<-sum(forecast_value[[i]][h])
  }
}

# plot the whole part
plot.ts(Ottseries, main="Annaulized Monthly New House Price Index", xlab='Year', ylab='NHPI', axes=F, type='l')
forecast.w.ts<-ts(forecast.w,frequency=12,start=c(2018,06))
ott.test[-1]
lines(forecast.w.ts,col='blue')
axis(1,at=seq(1998,2018,by=2))
axis(2)

# plot the forecasted part
plot.ts(forecast.w,main="6 steps ahead", col='blue',ylab='NHPI',xlab='Month',type='l',ylim=c(0.033,0.035))
forecast.w
# compare accuracy bettwen arima and emd-svr
measure_compare<-cbind(arima_measure,emd.svr.p.error)
# MAE, MASE, MdAE, MdSPE, RMdSPE, sMAPE, sMdAP, MdRAE
Ottseries

