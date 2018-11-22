

censoringProbability <- function(Y,treatment,completeCase){
  treated <- Y[treatment==1];
  treated.cc <- 1 - completeCase[treatment==1];
  treated.data <- as.data.frame(cbind(treated, treated.cc))
  control <- Y[treatment==0];
  control.cc <- 1 - completeCase[treatment==0];
  control.data <- as.data.frame(cbind(control, control.cc))
  colnames(treated.data) <- c("time","status"); colnames(control.data) <- c("time","status")
  treated.surv <- survfit(Surv(time,status)~1, treated.data)
  control.surv <- survfit(Surv(time,status)~1, control.data)
  treated.survest <- stepfun(treated.surv$time, c(1, treated.surv$surv))
  control.survest <- stepfun(control.surv$time, c(1, control.surv$surv))
  censoringProb <- vector(length = length(Y))
  #reconstruct the censoring probability vector
  for (i in 1:length(Y)){
    if (i %in% which(treatment == 1))
      censoringProb[i] <- treated.survest(Y[i])
    else
      censoringProb[i] <- control.survest(Y[i])
  }
  return(censoringProb)
}
