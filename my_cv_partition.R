my_cv_partition<-function(trn,is.ts=FALSE,cv=5){
  trn.length=nrow(trn$data)
  trn.width=ncol(trn$data)
  if(is.ts){ # for ts, sequential
    partition.number=cv*2-1 # rolling cv
    # make factor to partition
    f<-rep(1:partition.number,each=floor(trn.length/partition.number))
  }else{ # non cv partition
    # make factor to partition
    f<-ceiling(runif(trn.length,0,cv)) # cv=1 means no cv
  }
  # initialize 
  x<-y<-list()
  x.labels<-y.labels<-list()
  length(x)<-length(y)<-length(x.labels)<-length(y.labels)<-cv
  
  if(cv>1){# valid cv number
    for (c in 1:cv){
      y[[c]]<-trn$data[which(f==c),] # validation data
      x[[c]]<-trn$data[which(f!=c),] # training data
      if(is.null(nrow(trn$labels))){
        y.labels[[c]]<-trn$labels[which(f==c)] # validataion labels
        x.labels[[c]]<-trn$labels[which(f!=c)] # training labels
      }else{
        y.labels[[c]]<-trn$labels[which(f==c),,drop=F] # validataion labels
        x.labels[[c]]<-trn$labels[which(f!=c),,drop=F] # training labels
      }
    }
    return(list(trn.data=x, trn.labels=x.labels, val.data=y, val.labels=y.labels,k=cv))
  }else{#no cv
    return(NULL)
  }
}