susan1 <- function(n, p=2, outcome = T, pt = 1, py = 1, treatsize = .5){
  # Generate data
  # parameters for data generating
  # n is the  total size of the dataset
  asym <- .5 # whether treatment effects are distributed asymmetrically across treated and control
  propens <- .5 #treatment probability
  levsize <- 1
  sig <- .01
  # draw W
  w <- rbinom(n, 1, propens)
  w_ <- 1-w

  # draw X
  X <- matrix(rnorm(n * p), nrow = n, ncol = p)

  # generate treatment effects as function of X
  if (p<pt+py) print("error: p>=pt+py required")
  tau <- 0
  for (iii in 1:pt) {
    tau <- tau + treatsize*pmax(X[,iii],array(0,n))*(2*X[,iii])
    # tau <- tau + treatsize*pmax(X[,iii],array(0,n))
  }

  # generate average value of outcomes
  mu <- 0.5 * X[,1] + X[,2]

  # generate outcomes as function of treatment status, mu, tau, and noise
  y <- mu + 0.5 * (2 * w - 1) * tau + (asym-1)*(1-w)*tau + rnorm(n, 0, sig)
  y_ <- mu + 0.5*(2*w_ - 1) * tau + (asym-1)*(1-w_)*tau + rnorm(n, 0, sig)
  # y <- mu + asym*w*tau + (asym-1)*(1-w)*tau + rnorm(n,0,sig)
  # y_ <- mu + asym*w_*tau + (asym-1)*(1-w_)*tau + rnorm(n,0,sig)

  # create formulas for estimation
  # if modifying code for another dataset, need to simply name outcome variable y, treatment variable w, and
  # create a formula such as f with the list of x variables separated by "+"
  f <- ""
  nextx <- ""
  if (p>1) {
    for (ii in 1:(p-1)) {
      nextx <- paste("x",ii, sep="")
      if (ii==1) {name <- nextx}
      if (ii>1) {name <- c(name, nextx)}
      f <- paste(f, nextx, "+", sep="")
    }
    f <- paste(f, "x", ii+1, sep="")
  } else if (p==1) {
    f <- "x1"
  }

  for (ii in 1:p) {
    nextx <- paste("x",ii, sep="")
    if (ii==1) {name <- nextx}
    if (ii>1) {name <- c(name, nextx)}
  }

  name <- c( name,  "y", "w", "tau_true")

  tau_true <- (1-2*w)*(y_ - y)

  ntr <- round(.5*n)
  # nest <- round(.333*n)
  ntest <- n - ntr


  X<-data.frame(X)

  data <- data.frame(X[1:n,], y[1:n], w[1:n], tau_true[1:n])

  names(data)=name
  temp<-list(data,f)
  return(temp)
}
susan11 <- function(n, p=2, outcome = T, pt = 1, py = 1, treatsize = .5){
  # Generate data
  # parameters for data generating
  # n is the  total size of the dataset
  asym <- .5 # whether treatment effects are distributed asymmetrically across treated and control
  propens <- .5 #treatment probability
  levsize <- 1
  sig <- .01
  # draw W
  w <- rbinom(n, 1, propens)
  w_ <- 1-w

  # draw X
  X <- matrix(rnorm(n * p), nrow = n, ncol = p)

  # generate treatment effects as function of X
  if (p<pt+py) print("error: p>=pt+py required")
  tau <- 0
  for (iii in 1:pt) {
    # tau <- tau + treatsize*pmax(X[,iii],array(0,n))*(2*X[,iii])
    tau <- tau + treatsize*pmax(X[,iii],array(0,n))
    # tau <- tau + treatsize*X[,iii]
  }

  # generate average value of outcomes
  mu <- 0.5 * X[,1] + X[,2]

  # generate outcomes as function of treatment status, mu, tau, and noise
  y <- mu + 0.5 * (2 * w - 1) * tau + (asym-1)*(1-w)*tau + rnorm(n, 0, sig)
  y_ <- mu + 0.5*(2*w_ - 1) * tau + (asym-1)*(1-w_)*tau + rnorm(n, 0, sig)
  # y <- mu + asym*w*tau + (asym-1)*(1-w)*tau + rnorm(n,0,sig)
  # y_ <- mu + asym*w_*tau + (asym-1)*(1-w_)*tau + rnorm(n,0,sig)

  # create formulas for estimation
  # if modifying code for another dataset, need to simply name outcome variable y, treatment variable w, and
  # create a formula such as f with the list of x variables separated by "+"
  f <- ""
  nextx <- ""
  if (p>1) {
    for (ii in 1:(p-1)) {
      nextx <- paste("x",ii, sep="")
      if (ii==1) {name <- nextx}
      if (ii>1) {name <- c(name, nextx)}
      f <- paste(f, nextx, "+", sep="")
    }
    f <- paste(f, "x", ii+1, sep="")
  } else if (p==1) {
    f <- "x1"
  }

  for (ii in 1:p) {
    nextx <- paste("x",ii, sep="")
    if (ii==1) {name <- nextx}
    if (ii>1) {name <- c(name, nextx)}
  }

  name <- c( name,  "y", "w", "tau_true")

  tau_true <- (1-2*w)*(y_ - y)

  ntr <- round(.5*n)
  # nest <- round(.333*n)
  ntest <- n - ntr


  X<-data.frame(X)

  data <- data.frame(X[1:n,], y[1:n], w[1:n], tau_true[1:n])

  names(data)=name
  temp<-list(data,f)
  return(temp)
}


susan2 <- function(n, p = 10, outcome = T, pt = 4, py = 4, treatsize = .5){
  # Setting No.2 in Susan's PNAS paper
  # Generate data
  # parameters for data generating

  asym <- .5 # whether treatment effects are distributed asymmetrically across treated and control
  propens <- .5 #treatment probability
  sig = .01
  levsize <- 1

  # draw W
  w <- rbinom(n, 1, propens)
  w_ <- 1-w

  # draw X
  X <- matrix(rnorm(n * p), nrow = n, ncol = p)

  # generate treatment effects as function of X
  if (p<pt+py) print("error: p>=pt+py required")
  tau <- 0
  for (iii in 1:pt) {
    # tau <- tau + treatsize*pmax(X[,iii],array(0,n))
    tau <- tau + treatsize*pmax(X[,iii],array(0,n))*(2*X[,iii])
  }

  # generate average value of outcomes
  if (outcome == T){
    mu <- treatsize * rowSums(X[,1:pt]) + rowSums(X[,(pt+1):(pt+py)])
    # mu <- rowSums(X[,(pt+1):(pt+py)])

  }
  else
    mu <- 0


  # generate outcomes as function of treatment status, mu, tau, and noise
  y <- mu + 0.5 * (2 * w - 1) * tau + (asym-1)*(1-w)*tau + rnorm(n, 0, sig)
  y_ <- mu + 0.5*(2*w_ - 1) * tau + (asym-1)*(1-w_)*tau + rnorm(n, 0, sig)

  # y <- mu + asym*w*tau + (asym-1)*(1-w)*tau + rnorm(n,0,sig)
  # y_ <- mu + asym*w_*tau + (asym-1)*(1-w_)*tau + rnorm(n,0,sig)

  # create formulas for estimation
  # if modifying code for another dataset, need to simply name outcome variable y, treatment variable w, and
  # create a formula such as f with the list of x variables separated by "+"
  f <- ""
  nextx <- ""
  if (p>1) {
    for (ii in 1:(p-1)) {
      nextx <- paste("x",ii, sep="")
      if (ii==1) {name <- nextx}
      if (ii>1) {name <- c(name, nextx)}
      f <- paste(f, nextx, "+", sep="")
    }
    f <- paste(f, "x", ii+1, sep="")
  } else if (p==1) {
    f <- "x1"
  }

  for (ii in 1:p) {
    nextx <- paste("x",ii, sep="")
    if (ii==1) {name <- nextx}
    if (ii>1) {name <- c(name, nextx)}
  }

  name <- c( name,  "y", "w", "tau_true")

  tau_true <- (1-2*w)*(y_ - y)

  X<-data.frame(X)

  data <- data.frame(X[1:n,], y[1:n], w[1:n], tau_true[1:n])

  names(data)=name
  temp<-list(data,f)
  return(temp)
}
susan3 <- function(n, p = 20, outcome = T, pt = 4, py = 4, treatsize = .5){
  # Setting No.3 in Susan's PNAS paper
  #     Generate data
  # parameters for data generating
  asym <- .5 # whether treatment effects are distributed asymmetrically across treated and control
  propens <- .5 #treatment probability
  sig = .01
  levsize <- 1

  # draw W
  w <- rbinom(n, 1, propens)
  w_ <- 1-w

  # draw X
  X <- matrix(rnorm(n * p), nrow = n, ncol = p)

  # generate treatment effects as function of X
  if (p<pt+py) print("error: p>=pt+py required")
  tau <- 0
  for (iii in 1:pt) {
    tau <- tau + treatsize*pmax(X[,iii],array(0,n))*(2*X[,iii])
  }

  # generate average value of outcomes
  if (outcome == T){
    mu <- treatsize * rowSums(X[,1:pt]) + rowSums(X[,(pt+1):(pt+py)])
    # mu <- rowSums(X[,(pt+1):(pt+py)])

  }
  else
    mu <- 0

  # generate outcomes as function of treatment status, mu, tau, and noise
  y <- mu + 0.5 * (2 * w - 1) * tau + (asym-1)*(1-w)*tau + rnorm(n, 0, sig)
  y_ <- mu + 0.5*(2*w_ - 1) * tau + (asym-1)*(1-w_)*tau + rnorm(n, 0, sig)

  # y <- mu + asym*w*tau + (asym-1)*(1-w)*tau + rnorm(n,0,sig)
  # y_ <- mu + asym*w_*tau + (asym-1)*(1-w_)*tau + rnorm(n,0,sig)

  # create formulas for estimation
  # if modifying code for another dataset, need to simply name outcome variable y, treatment variable w, and
  # create a formula such as f with the list of x variables separated by "+"
  f <- ""
  nextx <- ""
  if (p>1) {
    for (ii in 1:(p-1)) {
      nextx <- paste("x",ii, sep="")
      if (ii==1) {name <- nextx}
      if (ii>1) {name <- c(name, nextx)}
      f <- paste(f, nextx, "+", sep="")
    }
    f <- paste(f, "x", ii+1, sep="")
  } else if (p==1) {
    f <- "x1"
  }

  for (ii in 1:p) {
    nextx <- paste("x",ii, sep="")
    if (ii==1) {name <- nextx}
    if (ii>1) {name <- c(name, nextx)}
  }

  name <- c( name,  "y", "w", "tau_true")

  tau_true <- (1 - 2 * w) * (y_ - y)


  data <- data.frame(X[1:n,], y[1:n], w[1:n], tau_true[1:n])

  names(data)=name
  temp<-list(data,f)
  return(temp)
}

nonlinear1 <- function(n, p = 10, outcome = T, pt = 4, py = 4, treatsize = .5){
  # Setting No.3 in Susan's PNAS paper
  #     Generate data
  # parameters for data generating
  asym <- .5 # whether treatment effects are distributed asymmetrically across treated and control
  propens <- .5 #treatment probability
  sig = .01
  levsize <- 1

  # draw W
  w <- rbinom(n, 1, propens)
  w_ <- 1-w

  # draw X
  X <- matrix(rnorm(n * p), nrow = n, ncol = p)

  # generate treatment effects as function of X
  if (p<pt+py) print("error: p>=pt+py required")
  tau <- 0
  for (iii in 1:pt) {
    if (iii <2){
      tau <- tau + treatsize * sin(X[,iii] - .5)
    }
    else
      tau <- tau + treatsize * X[,iii]^iii
  }

  # generate average value of outcomes
  if (outcome == T){
    mu <- treatsize * rowSums(X[,1:pt]) + rowSums(X[,(pt+1):(pt+py)])
    # mu <- rowSums(X[,(pt+1):(pt+py)])

  }
  else
    mu <- 0

  # generate outcomes as function of treatment status, mu, tau, and noise
  y <- mu + 0.5 * (2 * w - 1) * tau + (asym-1)*(1-w)*tau + rnorm(n, 0, sig)
  y_ <- mu + 0.5*(2*w_ - 1) * tau + (asym-1)*(1-w_)*tau + rnorm(n, 0, sig)

  # y <- mu + asym*w*tau + (asym-1)*(1-w)*tau + rnorm(n,0,sig)
  # y_ <- mu + asym*w_*tau + (asym-1)*(1-w_)*tau + rnorm(n,0,sig)

  # create formulas for estimation
  # if modifying code for another dataset, need to simply name outcome variable y, treatment variable w, and
  # create a formula such as f with the list of x variables separated by "+"
  f <- ""
  nextx <- ""
  if (p>1) {
    for (ii in 1:(p-1)) {
      nextx <- paste("x",ii, sep="")
      if (ii==1) {name <- nextx}
      if (ii>1) {name <- c(name, nextx)}
      f <- paste(f, nextx, "+", sep="")
    }
    f <- paste(f, "x", ii+1, sep="")
  } else if (p==1) {
    f <- "x1"
  }

  for (ii in 1:p) {
    nextx <- paste("x",ii, sep="")
    if (ii==1) {name <- nextx}
    if (ii>1) {name <- c(name, nextx)}
  }

  name <- c( name,  "y", "w", "tau_true")

  tau_true <- (1 - 2 * w) * (y_ - y)


  data <- data.frame(X[1:n,], y[1:n], w[1:n], tau_true[1:n])

  names(data)=name
  temp<-list(data,f)
  return(temp)
}

constant <- function(n, p = 10, outcome = T, pt = 4, py = 4, treatsize = .5){
  # Setting No.3 in Susan's PNAS paper
  #     Generate data
  # parameters for data generating
  asym <- .5 # whether treatment effects are distributed asymmetrically across treated and control
  propens <- .5 #treatment probability
  sig = .01
  levsize <- 1

  # draw W
  w <- rbinom(n, 1, propens)
  w_ <- 1-w

  # draw X
  X <- matrix(rnorm(n * p), nrow = n, ncol = p)

  # generate treatment effects as function of X
  if (p<pt+py) print("error: p>=pt+py required")
  tau <- treatsize * pt

  # generate average value of outcomes
  if (outcome == T){
    mu <- treatsize * rowSums(X[,1:pt]) + rowSums(X[,(pt+1):(pt+py)])
    # mu <- rowSums(X[,(pt+1):(pt+py)])

  }
  else
    mu <- 0

  # generate outcomes as function of treatment status, mu, tau, and noise
  y <- mu + 0.5 * (2 * w - 1) * tau + (asym-1)*(1-w)*tau + rnorm(n, 0, sig)
  y_ <- mu + 0.5*(2*w_ - 1) * tau + (asym-1)*(1-w_)*tau + rnorm(n, 0, sig)

  # y <- mu + asym*w*tau + (asym-1)*(1-w)*tau + rnorm(n,0,sig)
  # y_ <- mu + asym*w_*tau + (asym-1)*(1-w_)*tau + rnorm(n,0,sig)

  # create formulas for estimation
  # if modifying code for another dataset, need to simply name outcome variable y, treatment variable w, and
  # create a formula such as f with the list of x variables separated by "+"
  f <- ""
  nextx <- ""
  if (p>1) {
    for (ii in 1:(p-1)) {
      nextx <- paste("x",ii, sep="")
      if (ii==1) {name <- nextx}
      if (ii>1) {name <- c(name, nextx)}
      f <- paste(f, nextx, "+", sep="")
    }
    f <- paste(f, "x", ii+1, sep="")
  } else if (p==1) {
    f <- "x1"
  }

  for (ii in 1:p) {
    nextx <- paste("x",ii, sep="")
    if (ii==1) {name <- nextx}
    if (ii>1) {name <- c(name, nextx)}
  }

  name <- c( name,  "y", "w", "tau_true")

  tau_true <- (1 - 2 * w) * (y_ - y)


  data <- data.frame(X[1:n,], y[1:n], w[1:n], tau_true[1:n])

  names(data)=name
  temp<-list(data,f)
  return(temp)
}

# wagner1 is not implemented
wagner1 <- function(n){
  # Setting No.2 in Susan's PNAS paper

  # Generate data
  # parameters for data generating
  # total size of the dataset

  p <- 8 # number of total covariates
  pt <- 2 # number of covariates affecting treatment effects
  py <- 0 # number of covariates affecting outcomes but not treatment effects
  asym <- .5 # whether treatment effects are distributed asymmetrically across treated and control
  propens <- .5 #treatment probability
  sig = .01
  treatsize <- .5 # treatment effect size
  levsize <- 1

  # draw W
  w <- rbinom(n, 1, propens)
  w_ <- 1-w

  # draw X
  X <- matrix(rnorm(n * p), nrow = n, ncol = p)

  # generate treatment effects as function of X
  if (p<pt+py) print("error: p>=pt+py required")
  tau <- 0
  for (iii in 1:pt) {
    tau <- tau * (1 + 1 / (1+exp(-20*(X[,iii]-1/3))))
  }

  # generate outcomes as function of treatment status, mu, tau, and noise
  y <- 0.5 * (2 * w - 1) * tau + (asym-1)*(1-w)*tau + rnorm(n, 0, sig)
  y_ <- 0.5*(2*w_ - 1) * tau + (asym-1)*(1-w_)*tau + rnorm(n, 0, sig)


  # create formulas for estimation
  # if modifying code for another dataset, need to simply name outcome variable y, treatment variable w, and
  # create a formula such as f with the list of x variables separated by "+"
  f <- ""
  nextx <- ""
  if (p>1) {
    for (ii in 1:(p-1)) {
      nextx <- paste("x",ii, sep="")
      if (ii==1) {name <- nextx}
      if (ii>1) {name <- c(name, nextx)}
      f <- paste(f, nextx, "+", sep="")
    }
    f <- paste(f, "x", ii+1, sep="")
  } else if (p==1) {
    f <- "x1"
  }

  for (ii in 1:p) {
    nextx <- paste("x",ii, sep="")
    if (ii==1) {name <- nextx}
    if (ii>1) {name <- c(name, nextx)}
  }

  name <- c( name,  "y", "w", "tau_true")

  tau_true <- (1-2*w)*(y_ - y)

  ntr <- round(.5*n)
  # nest <- round(.333*n)
  ntest <- n - ntr


  X<-data.frame(X)

  data <- data.frame(X[1:n,], y[1:n], w[1:n], tau_true[1:n])

  names(data)=name
  temp<-list(data,f)
  return(temp)

}
wagner2 <- function(n, p=8){
  # Setting No.2 in Susan's PNAS paper

  # Generate data
  # parameters for data generating
  # n <- 5000 # total size of the dataset

  p <- 8 # number of total covariates
  pt <- 2 # number of covariates affecting treatment effects
  py <- 0 # number of covariates affecting outcomes but not treatment effects
  asym <- .5 # whether treatment effects are distributed asymmetrically across treated and control
  propens <- .5 #treatment probability
  sig = .01
  treatsize <- .5 # treatment effect size
  levsize <- 1

  # draw W
  w <- rbinom(n, 1, propens)
  w_ <- 1-w

  # draw X
  X <- matrix(rnorm(n * p), nrow = n, ncol = p)

  # generate treatment effects as function of X
  if (p<pt+py) print("error: p>=pt+py required")
  tau <- 1
  for (iii in 1:pt) {
    tau <- tau * (1 + 1 / (1+exp(-20*(X[,iii]-1/3))))
  }

  # generate outcomes as function of treatment status, mu, tau, and noise
  y <- 0.5 * (2 * w - 1) * tau + (asym-1)*(1-w)*tau + rnorm(n, 0, sig)
  y_ <- 0.5*(2 * w_ - 1) * tau + (asym-1)*(1-w_)*tau + rnorm(n, 0, sig)


  # create formulas for estimation
  # if modifying code for another dataset, need to simply name outcome variable y, treatment variable w, and
  # create a formula such as f with the list of x variables separated by "+"
  f <- ""
  nextx <- ""
  if (p>1) {
    for (ii in 1:(p-1)) {
      nextx <- paste("x",ii, sep="")
      if (ii==1) {name <- nextx}
      if (ii>1) {name <- c(name, nextx)}
      f <- paste(f, nextx, "+", sep="")
    }
    f <- paste(f, "x", ii+1, sep="")
  } else if (p==1) {
    f <- "x1"
  }

  for (ii in 1:p) {
    nextx <- paste("x",ii, sep="")
    if (ii==1) {name <- nextx}
    if (ii>1) {name <- c(name, nextx)}
  }

  name <- c( name,  "y", "w", "tau_true")

  tau_true <- (1 -2 * w)*(y_ - y)

  ntr <- round(.5*n)
  # nest <- round(.333*n)
  ntest <- n - ntr


  X<-data.frame(X)

  data <- data.frame(X[1:n,], y[1:n], w[1:n], tau_true[1:n])

  names(data)=name
  temp<-list(data,f)
  return(temp)

}
wagner3 <- function(n, p=8, outcome = F, pt = 2, treatsize = 5){
  # Setting No.2 in Susan's PNAS paper

  # Generate data
  # parameters for data generating
  # n <- 5000 # total size of the dataset

  # p <- 8 # number of total covariates
  # pt <- 2 # number of covariates affecting treatment effects
  py <- 4 # number of covariates affecting outcomes but not treatment effects
  asym <- .5 # whether treatment effects are distributed asymmetrically across treated and control
  propens <- .5 #treatment probability
  sig = .01
  levsize <- 1

  # draw W
  w <- rbinom(n, 1, propens)
  w_ <- 1-w

  # draw X
  X <- matrix(rnorm(n * p), nrow = n, ncol = p)

  # generate treatment effects as function of X
  if (p<pt+py) print("error: p>=pt+py required")
  tau <- 1
  for (iii in 1:pt) {
    # tau <- tau * (1 + 2 / (1+exp(-12*(X[,iii]-1/2))))
    tau <- tau * (1 + 2 / (1+exp(-12*(pmax(X[,iii],array(0,n))-1/2))))
  }
  if (outcome == T){
    # mu <- treatsize * rowSums(X[,1:pt]) + rowSums(X[,(pt+1):(pt+py)])
    mu <- rowSums(X[,(pt+1):(pt+py)])
  }
  else
    mu <- 0

  # generate outcomes as function of treatment status, mu, tau, and noise
  y <- mu + 0.5 * (2 * w - 1) * tau + (asym-1)*(1-w)*tau + rnorm(n, 0, sig)
  y_ <- mu + 0.5*(2*w_ - 1) * tau + (asym-1)*(1-w_)*tau + rnorm(n, 0, sig)
  # y <-  0.5 * (2 * w - 1) * tau + (asym-1)*(1-w)*tau + rnorm(n, 0, sig)
  # y_ <-  0.5*(2*w_ - 1) * tau + (asym-1)*(1-w_)*tau + rnorm(n, 0, sig)


  # create formulas for estimation
  # if modifying code for another dataset, need to simply name outcome variable y, treatment variable w, and
  # create a formula such as f with the list of x variables separated by "+"
  f <- ""
  nextx <- ""
  if (p>1) {
    for (ii in 1:(p-1)) {
      nextx <- paste("x",ii, sep="")
      if (ii==1) {name <- nextx}
      if (ii>1) {name <- c(name, nextx)}
      f <- paste(f, nextx, "+", sep="")
    }
    f <- paste(f, "x", ii+1, sep="")
  } else if (p==1) {
    f <- "x1"
  }

  for (ii in 1:p) {
    nextx <- paste("x",ii, sep="")
    if (ii==1) {name <- nextx}
    if (ii>1) {name <- c(name, nextx)}
  }

  name <- c( name,  "y", "w", "tau_true")

  tau_true <- (1-2*w)*(y_ - y)

  ntr <- round(.5*n)
  # nest <- round(.333*n)
  ntest <- n - ntr


  X<-data.frame(X)

  data <- data.frame(X[1:n,], y[1:n], w[1:n], tau_true[1:n])

  names(data)=name
  temp<-list(data,f)
  return(temp)

}

min1 <- function(n){
  # Setting No.1 in Min's preprint
  # Generate data
  # parameters for data generating
  p <- 20 # number of total covariates
  pt <- 2 # number of covariates affecting treatment effects
  py <- 4 # number of covariates affecting outcomes but not treatment effects
  asym <- .5 # whether treatment effects are distributed asymmetrically across treated and control
  propens <- .5 #treatment probability
  sig = .01
  treatsize <- .5 # treatment effect size
  levsize <- 1

  # draw W
  w <- rbinom(n, 1, propens)
  w_ <- 1-w

  # draw X
  X.1 <- matrix(rnorm(n * 11), nrow = n, ncol = 11)
  X.2 <- matrix(rbinom(n * 9, 1, 0.5), nrow = n, ncol = 9)
  X <- cbind(X.1, X.2)

  # generate treatment effects as function of X
  if (p<pt+py) print("error: p>=pt+py required")
  tau <- 0
  for (iii in 1:pt) {
    tau <- tau + treatsize*pmax(X[,iii],array(0,n))
  }

  # generate average value of outcomes
  mu <- 0.5 * rowSums(X[,1:pt]) + rowSums(X[,(pt+1):(pt+py)])

  # generate outcomes as function of treatment status, mu, tau, and noise
  y <- mu + 0.5 * (2 * w - 1) * tau + (asym - 1) * (1 - w) * tau + rnorm(n, 0, sig)
  y_ <- mu + 0.5*(2*w_ - 1) * tau + (asym - 1) * (1 - w_) * tau + rnorm(n, 0, sig)

  # y <- mu + asym*w*tau + (asym-1)*(1-w)*tau + rnorm(n,0,sig)
  # y_ <- mu + asym*w_*tau + (asym-1)*(1-w_)*tau + rnorm(n,0,sig)

  # create formulas for estimation
  # if modifying code for another dataset, need to simply name outcome variable y, treatment variable w, and
  # create a formula such as f with the list of x variables separated by "+"
  f <- ""
  nextx <- ""
  if (p>1) {
    for (ii in 1:(p-1)) {
      nextx <- paste("x",ii, sep="")
      if (ii==1) {name <- nextx}
      if (ii>1) {name <- c(name, nextx)}
      f <- paste(f, nextx, "+", sep="")
    }
    f <- paste(f, "x", ii+1, sep="")
  } else if (p==1) {
    f <- "x1"
  }

  for (ii in 1:p) {
    nextx <- paste("x",ii, sep="")
    if (ii==1) {name <- nextx}
    if (ii>1) {name <- c(name, nextx)}
  }

  name <- c( name,  "y", "w", "tau_true")

  tau_true <- (1-2*w)*(y_ - y)

  data <- data.frame(X[1:n,], y[1:n], w[1:n], tau_true[1:n])
  names(data)=name
  temp<-list(data,f)
  return(temp)
}
