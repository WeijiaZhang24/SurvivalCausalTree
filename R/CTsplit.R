# The 'evaluation' function.  Called once per node.
evaluation <- function(y, wt, parms) {
  # pre-process wt to treatment indicator and propensity score
  indicator<-floor(wt)-1
  pscore<-wt-floor(wt)

  treated.y<-y[indicator==1]
  control.y <- y[indicator==0]

  cate.node <- CATE(y,pscore,indicator,parms)
  # score <- sum((TOT-cate.node))^2/length(y)
  # list(label= cate.node, deviance = 1/result)
  # list(label= cate.node, deviance = 1/CATE(y,pscore,indicator,list(pscore=F)))
  list(label= cate.node, deviance = length(wt)-sqrt(length(wt)))
  # list(label= cate.node, deviance = score)
}

# The split function, where most of the work occurs.
#   Called once per split variable per node.
CATE<-function(y,pscore,indicator,parms){
  # estimate the CATE for a node, given treatment indicator and propensity score
  treated.y<-y[indicator==1]
  treated.pscore<-pscore[indicator==1]
  control.y <- y[indicator==0]
  control.pscore <- pscore[indicator==0]
  if (length(treated.y)>=1 & length(control.y)>=1)
    result<-sum(treated.y/treated.pscore) / sum(1 / treated.pscore) -
    sum(control.y/(1-control.pscore))/sum(1/(1-control.pscore))
  else
    result<-0
  return(result)
}
estimateCATE<-function(y,indicator,pscore,position,parms){
  # this function estimate the Causal Tree goodness of fit, given the position of split
  # y: outcome variable
  # indicator: treatment indicator
  # position: the position of the split, the fixed feature is ordered, so the split point is completely
  #           determined by position

    y.left<-y[1:position];  y.right<-y[(position+1):length(y)]
    indicator.left<-indicator[1:position];  indicator.right<-indicator[(position+1):length(y)]
    pscore.left<-pscore[1:position];  pscore.right<-pscore[(position+1):length(y)]

    left<-CATE(y.left,pscore.left,indicator.left,parms)
    right<-CATE(y.right,pscore.right,indicator.right,parms)
    # all<-CATE(y,pscore,indicator,parms)
    # result<-abs((result.left^2*length(y.left)+result.right^2*length(y.right))/length(y))
    # result<-left^2+right^2-all^2
    result<-left^2 + right^2 + 2*left*right
    # result<-abs(left)*length(y.left)+abs(right)*length(y.right)
    # result<-exp(abs(left)+abs(right)-abs(all))
}

splitting <- function(y, wt, x, parms, continuous) {
  # indicator<-parms$indicator
  indicator<-floor(wt)-1
  pscore<-wt-floor(wt)

  n <- length(y)
  if (continuous) {
    # continuous x variable
    goodness<-vector(,n-1)
    for (i in 1:(length(y)-1)){
      goodness[i]<-estimateCATE(y,indicator,pscore,i,parms)
    }
    dir<-rep(-1,length(y)-1)
    list(goodness= goodness, direction=dir)
  }
  else{
    # categorical variables
    # only consider two level factors, factor with multiple level should be converted

    x.unique<-sort(unique(x))
    y.split<-split(y,x)
    indicator.split<-split(indicator,x)
    pscore.split<-split(pscore,x)
    CATE.all <- CATE(y,pscore,indicator,parms)
    # CATE.split<-vector()
    fit<-vector(,length(x.unique)-1)

    for (i in 1:(length(x.unique)-1)){
      y.left <- unlist(y.split[1:i])
      indicator.left <- unlist(indicator.split[1:i])
      pscore.left <- unlist(pscore.split[1:i])

      y.right <- unlist(y.split[(i+1):length(y.split)])
      indicator.right <- unlist (indicator.split[(i+1):length(indicator.split)])
      pscore.right <- unlist (pscore.split[(i+1):length(pscore.split)])

      # cate.split <- CATE(y.left,pscore.left,indicator.left,parms)^2+CATE(y.right,pscore.right,indicator.right,parms)^2-CATE.all^2
      cate.split <- abs(CATE(y.left,pscore.left,indicator.left,parms))+abs(CATE(y.right,pscore.right,indicator.right,parms))-abs(CATE.all)
      fit[i]<-cate.split
    }

#     for (i in 1:length(x.unique)){
#       CATE.split[i]<-CATE(y.split[[i]],pscore.split[[i]],indicator.split[[i]])
#     }
#     ord<-order(abs(CATE.split))
#     list(goodness= CATE.split[1]^2+CATE.split[2]^2-CATE.all^2,
#          direction = x.unique[ord])
    list(goodness= fit, direction = x.unique)
  }
}

# The init function:
init <- function(y, offset, parms, wt) {
  if (!is.null(offset)) y <- y-offset
  list(y=y, parms=0, numresp=1, numy=1, wt=wt,
       summary= function(yval, dev, wt, ylevel, digits ) {
         paste("  mean=", format(signif(yval, digits)),
               ", MSE=" , format(signif(dev/wt, digits)),
               sep='')
       })
}


splitlist <- list(eval=evaluation, split=splitting, init=init)


