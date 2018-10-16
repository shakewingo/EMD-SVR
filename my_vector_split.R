my_vector_split<-function(ts,HORIZON=1,LAG=2,RATIO=0.7,is.point=F)
{
  if(is.point){#point forecasting, need to delect intermediate data
    org.HORIZON=HORIZON
    HORIZON=max(HORIZON)
  }
  # training data
  x.width<-LAG+HORIZON
  x.length<-length(ts)-x.width+1
  x<-array(0,dim=c(x.length, x.width))
  
  for (i in (1:x.width))
  {
    x[,i]=ts[i:(i+x.length-1)]
  }
  
  historical<-x[,(1:LAG),drop=F]
  future<-x[,LAG+(1:HORIZON),drop=F]
  
  # split according to ratio
  # trn.idx<-1:round(nrow(historical)*RATIO)
  tst.length<-length(ts)-round(length(ts)*RATIO)
  tst.idx<-(nrow(historical)-tst.length+1):nrow(historical)
  
  trn.data<-as.matrix(historical)[-tst.idx,,drop=F]
  trn.labels<-as.matrix(future)[-tst.idx,,drop=F]
  if(tst.length==0){# no tst
    tst.data=NULL
    tst.labels=NULL
  }else{
    tst.data<-as.matrix(historical)[tst.idx,,drop=F]
    tst.labels<-as.matrix(future)[tst.idx,,drop=F]
  }
  
  if(is.point){ # discard unused points
    trn.labels<-trn.labels[,org.HORIZON,drop=F]
    if(tst.length!=0)
      tst.labels<-tst.labels[,org.HORIZON,drop=F]
  }
  
  trn<-list(data=trn.data,labels=trn.labels)
  tst<-list(data=tst.data,labels=tst.labels)
  
  return(list(trn=trn,tst=tst))
}


my_scale_ts<-function(ts, scaled.range=c(0,1),ts.range=range(ts))
{
  ts.max<-max(ts.range)
  ts.min<-min(ts.range)
  d.max<-max(scaled.range)
  d.min<-min(scaled.range)
  
  scaled.ts<-(ts-ts.min)*(d.max-d.min)/(ts.max-ts.min)+d.min
  return(list(scaled=scaled.ts, range=c(ts.min,ts.max)))
}
