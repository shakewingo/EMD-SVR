# The data has already been seaonsal adjusted by a monthly frequeny index
#install.packages('dplyr')
#install.packages('stringr')
#install.packages('psych')
#install.packages('tidyverse')
library(dplyr)
library(stringr)
library('scales')
library(psych)
library(tidyverse)
NHPI<-read.csv(file='/Users/yingyao/Desktop/Projects/CQADS/House Price - TS/New Housing Price Monthly.csv',header=TRUE)
# REF_DATE: 30 years from 1988-01 to 2018-06, New.housing.price.indexes: Total (house and land), GEO, Status, Value
Dataset1<-filter(NHPI, str_detect(substr(REF_DATE,1,4), '1998|1999|2000|2001|2002|2003|2004|2005|2006|2007|2008|2009|2010|2011|2012|2013|2014|2015|2016|2017|2018'),
                 New.housing.price.indexes=="Total (house and land)",
                 str_detect(GEO,'Ottawa-Gatineau, Ontario part|Toronto|Vancouver'),
                 VALUE!='NA'
                 )
#Dataset2<-filter(Dataset1, REF_DATE != "2013-01", REF_DATE != "2013-02", REF_DATE != "2013-03", REF_DATE != "2013-04", REF_DATE != "2013-05", REF_DATE != "2013-06")
Dataset2<-select(Dataset1, REF_DATE, GEO, New.housing.price.indexes, VALUE)

# basic statistics for original dataset (NHPI)
summary<-do.call('rbind',describeBy(Dataset2$VALUE,Dataset2$GEO,digits=3))

# add annualized monthly rate
Dataset4<-Dataset2 %>%
  group_by(GEO) %>%
  mutate(lag = lag(VALUE)) %>%
  mutate(Growth_Pct = ((VALUE/ lag)^12)-1)
Dataset4$Growth_Pct[Dataset4$REF_DATE == "1998-01"] <- 0

# basic statistics for new dataset (Growth Rate)
summary2<-do.call('rbind',describeBy(Dataset4$Growth_Pct,Dataset4$GEO,digits=2))
statistics<-select(summary2,mean,sd,min,max,skew,kurtosis)
Dataset4<-Dataset4 %>% ungroup(GEO)

#############################################################################
# PLOTS
#############################################################################
# Use ggplot to plot from one dataframe, use aes(.., group=) to group them
library(reshape2)
library(ggplot2)

# Group by cities
Ottdata<-filter(Dataset4,str_detect(GEO,'Ottawa-Gatineau'))  #Should be change until 2018-06
Ottdata<-select(Ottdata,Growth_Pct,VALUE) #from 1998-01 to 2018-01
OttValue_Series<-ts(Ottdata$VALUE, frequency=12, start=c(1998,01))
Ottseries<-ts(Ottdata$Growth_Pct, frequency=12, start=c(1998,01))
Trtdata<-filter(Dataset4,str_detect(GEO,'Toronto')) 
Trtdata<-select(Trtdata,Growth_Pct,VALUE)
TrtValue_Series<-ts(Trtdata$VALUE, frequency=12, start=c(1998,01))
Trtseries<-ts(Trtdata$Growth_Pct, frequency=12, start=c(1998,01))
Vbcdata<-filter(Dataset4,str_detect(GEO,'Vancouver'))  #Should be change until 2018-06
Vbcdata<-select(Vbcdata,Growth_Pct,VALUE)
VbcValue_Series<-ts(Vbcdata$VALUE, frequency=12, start=c(1998,01))
Vbcseries<-ts(Vbcdata$Growth_Pct, frequency=12, start=c(1998,01))

date=unique(format(as.Date(paste(Dataset4$REF_DATE,sep='-','01'),'%Y-%m-%d'),'%Y-%m'))

# plot NHPI
plot.ts(OttValue_Series, main="Otawa NHPI Value", xlab='Year', ylab='NHPI', axes=F, type='l')
axis(1,at=seq(1998,2018,by=2))
axis(2, at=seq(55,105,by=10))


# plot Annualized Monthly Growth
plot.ts(Ottseries, main="Ottawa Annaulized Monthly NHPI", xlab='Year', ylab='Growth Percent', axes=F, type='l')
axis(1,at=seq(1998,2018,by=2))
axis(2, at=seq(-0.10,0.30, by=0.05))

plot.ts(Vbcseries, main="Vancouver Annaulized Monthly NHPI", xlab='Year', ylab='Growth Percent', axes=F, type='l')
axis(1,at=seq(1998,2018,by=2))
axis(2)

##############################################################
####################### Build ARIMA ##########################
##############################################################

#install.packages('forecast',dependencies=TRUE)
#install.packages('astsa')
library('forecast')
library('astsa')
library('tseries')
#install.packages('fUnitRoots')
#install.packages("https://cran.r-project.org/src/contrib/Archive/RcppArmadillo/RcppArmadillo_0.6.100.0.0.tar.gz", repos=NULL, type="source")
#install.packages('forecast')
#ott.test.start<-c(2012,05)

# split into in-sample and out-of-sample
RATIO=0.7    # train and test proportion
len_train=round(length(Ottseries)*RATIO)
tt <- seq(1, length(Ottseries), by=1) # for whole set time series

ott.test<-window(ts(Ottdata$Growth_Pct),start=(len_train+1))
ott.test.ts<-ts(ott.test, frequency=12, start=c(2012,05))
ott.train<-window(ts(Ottdata$Growth_Pct),end=(len_train))
plot.ts(Ottseries,axes=F)
axis(2)
axis(1, at = seq(from=1998, to=2018, by=4))
box()

# test stationary: unit root tests-ADF test: for general ARMA(p,q), determine if trending data has constant and time trend stationary
adf.test(ott.train,alternative = "stationary")
ott.train.diff1<-diff(ott.train,differences=1)
adf.test(ott.train.diff1,alternative = "stationary")
plot.ts(ott.train.diff1)
acf2(ott.train.diff1)
acf2(ott.train)
# ARIMA(0,1,1), \hat Y_{d_t} = 0.551 Y_{t-1} - 0.2496 e_{t-1} + E
# fit the model
arima.model<-auto.arima(ott.train,seasonal=FALSE)
(p<-length(arima.model$model$phi))
(d<-length(arima.model$model$Delta))
(q<-length(arima.model$model$theta))
arima.p<-data.frame(p=p,d=d,q=q)
ottseriesarima <- arima(ott.train, order=c(p,d,q))

# evaluate the residuals
tsdisplay(residuals(ottseriesarima), main=bquote((.(arima.p$p)~','~.(arima.p$d)~','~.(arima.p$q)) ~'Model Residuals'))

# validate test set
fit= arima(ott.train,order=c(p,d,q))
fcast <- forecast(fit,h=length(ott.test))
fcast.ts<-ts(fcast$mean,frequency=12,start=c(2012,05))
fcast.l.ts<-ts(fcast$lower[,2],frequency=12,start=c(2012,05))
fcast.h.ts<-ts(fcast$upper[,2],frequency=12,start=c(2012,05))
plot.ts(Ottseries, main="Annaulized Monthly New House Price Index", xlab='Year', ylab='NHPI', axes=F, type='l')
lines(fcast.ts,col='red')
lines(fcast.l.ts,col = 'grey')
lines(fcast.h.ts,col = 'grey')
axis(1,at=seq(1998,2018,by=2))
axis(2)

# calulate error
arima_measure<-my_forecasting_measure(fcast.ts,ott.test.ts)
#AIC = N ln(RMSE) + 2N(k + 1)ln(N)
#BIC = N ln(RMSE) + (k + 1)ln(N) − N ln(N)
