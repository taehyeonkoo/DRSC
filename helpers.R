# Generate AR series
generate.AR.series <- function(T01,rho){
  u           <- rep(NA,T01)
  epsl        <- rnorm(T01)*sqrt(1-rho^2) 
  startvalue  <- rnorm(1)
  u[1]        <- rho*startvalue + epsl[1]
  for (t in 2:T01){
    u[t] <- rho*u[(t-1)]+epsl[t]
  }
  return(u)
}
A1gen <- function(rho,p){
  A1=matrix(0,p,p)
  for(i in 1:p){
    for(j in 1:p){
      A1[i,j]<-rho^(abs(i-j))
    }
  }
  return(A1)
}

vecl <- function(A){
  A[upper.tri(A)] <- NA
  A.vec <- as.vector(A)[!is.na(A)]
  return(A.vec)
}

sc.repro <- function(Y1,Y0,T0,cond = F,
                     M = 100,eps =1e-3,method = "qreg",alpha = 0.05,reg = T){
  Y1.pre <- Y1[1:T0]
  Y0.pre <- Y0[1:T0,,drop = F]
  Y1.post <- Y1[-(1:T0)]
  Y0.post <- Y0[-(1:T0),,drop = F]
  # T0 <- nrow(Y0.pre)
  J <- ncol(Y0.pre)
  T1 <- nrow(Y0.post)

  Y0.vec <- vecl(Y0.pre)

  # J <- dim(Y0)[2]
  e <- matrix(1,1,J)
  f <- 1
  g <- diag(x=1,J,J)
  h <- matrix(0,J,1)
  w.hat <- limSolve::lsei(A=Y0,B=Y1,E=e,F=f,G=g,H=h,type=2)$X
  sc.resid <- Y1-Y0%*%w.hat
  # sc.resid <- sc.fit$u.hat
  # sc.beta <- sc.fit$w.hat
  if (reg) {
    sigma.u <- sqrt(mean(sc.resid^2))
    sigma.x <- apply(Y0.pre,2,sd)
    # T0 <- nrow(Y0.pre)
    c <- sigma.u/min(sigma.x)
    thres <- c*sqrt(log(T0))/sqrt(T0)

  } else{
    thres <- 0
  }


  # Y0.vec <- matrix(0,nrow = T0,ncol = J*(J+1)/2)
  # for (t in 1:T0) {
  #   temp <- t(Y0.pre[t,,drop = F])%*%Y0.pre[t,,drop = F]
  #   temp[upper.tri(temp)] <- NA
  #   Y0.vec[t,] <-  as.vector(temp)[!is.na(temp)]
  #
  # }
  #
  Y0.vec <- vecl(Y0.pre)
  tau0 <- 0.2
  Y0.vec.cov <- cov(Y0.vec) # covariance of vecl(cov(Y0))
  d0 <- max(max(apply(abs(Y0.vec.cov),1,sum))*tau0,1)
  Y0.cov <- cov(Y0.pre)
  Y0.cov[upper.tri(Y0.cov)] <- NA
  Y0.cov <- as.vector(Y0.cov)[!is.na(Y0.cov)]

  Y1Y0.vec <- (Y1.pre*Y0.pre)
  Y1Y0.mean <- apply(Y1Y0.vec,2,mean)
  Y1Y0.cov <- cov(Y1Y0.vec)
  tau.mat <- matrix(0,nrow = M, ncol = T1)
  CI.tau <-  matrix(0,nrow = M, ncol = 2*T1)
  CI.in.sample = matrix(0,nrow = M, ncol = 2*T1)

  resid.mat <- matrix(0,nrow = M,ncol = T0)
  screen.idx <- vector()
  Sigma.cov <- Y0.vec.cov/T0+d0/T0*diag(nrow(Y0.vec.cov))
  alpha0 <- 0.01
  value <- vector()

  for (m in 1:M) {
    if (cond) {
      XX.repro <- Y0.cov
    } else{
      XX.repro <- MASS::mvrnorm(mu = Y0.cov,Sigma = Sigma.cov)
    }
    ### Check m to screen ###
    value[m] <- max(abs(XX.repro-Y0.cov)/sqrt(diag(Sigma.cov)))
    if (max(abs(XX.repro-Y0.cov)/sqrt(diag(Sigma.cov)))<1.1*qnorm(1-alpha0/(J*(J+1)))) {
      screen.idx <- c(screen.idx,m)
    }
    XY.repro <- MASS::mvrnorm(mu = Y1Y0.mean,Sigma = Y1Y0.cov/T0)
    tmp <- matrix(0,nrow = J,ncol = J)
    tmp[lower.tri(tmp,diag = T)] <- XX.repro
    # Symmetric
    for(l in 2:J) {
      for(k in 1:(l-1)) {
        tmp[k,l] = tmp[l,k]
      }
    }

    Diag.XX<-diag(eigen(tmp)$values)
    for(ind in 1:J){
      Diag.XX[ind,ind]<-max(tmp[ind,ind],eps)
    }
    # eigen(XX.positive)
    XX.positive<-eigen(tmp)$vectors%*%Diag.XX%*%t(eigen(tmp)$vectors)
    Amat <- t(rbind(rep(1,J),diag(nrow =J)))
    bvec <- c(1,rep(0,J))
    w.repro <- quadprog::solve.QP(Dmat = XX.positive, dvec = XY.repro, Amat = Amat, bvec = bvec, meq = 1)$solution

    ### In-sample Uncertainty ###
    v.repro <- w.repro>thres
    Y0.repro <- Y0.pre[,v.repro,drop = F]
    if (length(Y0.repro)==0) {
      next

    }
    Y0.post.repro <- Y0.post[,v.repro,drop = F]
    fit.repro <- lm(Y1.pre~0+Y0.repro) # no intercept
    var.repro <- sandwich::vcovHC(fit.repro)/T0
    sd.repro <- sqrt(diag(var.repro)) # Heteroscedastic error
    beta.repro <- fit.repro$coefficients
    # construct confidence interval for each beta.j by the distribution of u.repro
    # beta.mat[m,] <- sc.repro
    var.in.sample <- diag(Y0.post.repro%*%var.repro%*%t(Y0.post.repro))
    CI.in.sample[m,] <- as.vector(t(cbind( Y0.post.repro%*%beta.repro-qnorm(1-alpha/2)*var.in.sample,
                                           Y0.post.repro%*%beta.repro+qnorm(1-alpha/2)*var.in.sample)))

    tau.t.repro <- Y1.post-Y0.post%*%w.repro
    resid.pre <- Y1.pre - Y0.pre%*%w.repro
    resid.mat[m,] <- resid.pre
    bound <-  scpi.out(resid.pre,Y0.pre,Y0.post,method = method,alpha = alpha,weight = w.repro)
    lb <- bound$lb; ub <- bound$ub
    tau.t.repro <- Y1.post-Y0.post%*%w.repro

    # construct confidence interval for tau

    tau.mat[m,] <- tau.t.repro
    CI.tau[m,]<- as.vector(t(cbind(tau.t.repro-ub,tau.t.repro-lb)))
  }
    CI.lower <- apply(CI.tau[, c(TRUE,FALSE),drop = FALSE],2,min)
    CI.upper <- apply(CI.tau[,!c(TRUE,FALSE),drop = FALSE],2,max)
    CI_tau <- cbind(CI.lower,CI.upper)

    CI.in.lower <- apply(CI.in.sample[, c(TRUE,FALSE),drop = FALSE],2,min)
    CI.in.upper <- apply(CI.in.sample[,!c(TRUE,FALSE),drop = FALSE],2,max)
    CI_in_sample <- cbind(CI.in.lower,CI.in.upper)


    CI.lower.screen <- apply(CI.tau[screen.idx, c(TRUE,FALSE),drop = FALSE],2,min)
    CI.upper.screen <- apply(CI.tau[screen.idx,!c(TRUE,FALSE),drop = FALSE],2,max)
    CI_tau.screen <- cbind(CI.lower.screen,CI.upper.screen)

  return(list(CI_tau = CI_tau,CI_in_sample = CI_in_sample,CI.tau = CI.tau,screen.idx = screen.idx,value = value))
}

# # no intercept
sc <- function(Y0,X0, intercept = FALSE, u = NULL){
  if (length(u)==0) {
    if (intercept) {
      J <- dim(X0)[2]
      X0 <- cbind(1,X0)
      e <- cbind(0,matrix(1,1,J))
      f <- 1
      g <- diag(x=1,J+1,J+1)
      g <- g[-1,]
      h <- matrix(0,J,1)
      w.hat <- limSolve::lsei(A=X0,B=Y0,E=e,F=f,G=g,H=h,type=2)$X
      u.hat <- Y0-X0%*%w.hat
      return(list(u.hat=u.hat,w.hat=w.hat))
    } else {
      J <- dim(X0)[2]
      e <- matrix(1,1,J)
      f <- 1
      g <- diag(x=1,J,J)
      h <- matrix(0,J,1)
      w.hat <- limSolve::lsei(A=X0,B=Y0,E=e,F=f,G=g,H=h,type=2)$X
      u.hat <- Y0-X0%*%w.hat
      return(list(u.hat=u.hat,w.hat=w.hat))
    }
  } else {
    if (intercept) {
      J <- dim(X0)[2]
      X0 <- cbind(1,X0)
      e <- cbind(0,matrix(1,1,J),0)
      f <- 1
      g <- cbind(diag(x=1,J+1,J+1),rep(0,J+1))
      g <- g[-1,]
      h <- matrix(0,J,1)
      w.hat <- limSolve::lsei(A=cbind(X0,u),B=Y0,E=e,F=f,G=g,H=h,type=2)$X
      u.hat <- Y0-X0%*%w.hat
      return(list(u.hat=u.hat,w.hat=w.hat))
    } else{
      J <- dim(X0)[2]
      e <- cbind(matrix(1,1,J),0)
      f <- 1
      g <- cbind(diag(x=1,J,J),rep(0,J))
      h <- rbind(matrix(0,J,1))
      temp <- limSolve::lsei(A=cbind(X0,u),B=Y0,E=e,F=f,G=g,H=h,type=2)$X
      w.hat <- c(temp[1:J])
      u.hat <- Y0-X0%*%w.hat
      sigma <- c(temp[J+1])
      return(list(u.hat=u.hat,w.hat=w.hat,sigma = sigma))
    }
  }
}

DRSC <- function(Y0,Y1,X0,X1,lambda = NULL, intercept = F, alpha = 0.01, 
                 gamma = 0.2,M = 100, CVXR = F, step = 0,
                 cond = F,Inference = F){
  
  T0 <- length(Y0)
  T1 <- length(Y1)
  N <- ncol(X0)
  
 
  # lambda0 <- max(abs(1/T0*t(X0)%*%(Y0-X0%*%w.hat)))
  if (length(lambda)==0) {
    lambda <- 0
  }
  
  X0.scale <- t(apply(X0,1,function(x){x/sqrt(apply(X0^2,2,mean))}))
  
  SC.fit <- sc(Y0,X0,intercept = intercept)
  w.hat <- SC.fit$w.hat
  SC.tau <- mean(Y1-X1%*%w.hat)
  sig <- sqrt(mean((SC.fit$u.hat)^2))
  
  ## Define sample mean ###
  ## For Constraint Set ##
  # t(X0)%*%X0
  X0.cov <- t(X0)%*%X0/T0
  
  Y0X0.vec <- (Y0*X0)
  Y0X0.mean <- apply(Y0X0.vec,2,mean)
  Y0X0.scale.mean <- apply((Y0*X0.scale),2,mean)
  ## Post-treatment period
  EY1 <- mean(Y1)
  EX1 <- apply(X1,2,mean)
  VarY1 <- var(Y1)
  CovX1 <- cov(X1)
  
  # X1.scale <- scale(X1, center=TRUE, scale=FALSE)
  # Y1.scale <- scale(Y1,center = T,scale = F)
  ### Sum-to-one ###
  e <- cbind(matrix(1,1,N))
  f <- 1
  
  # Normalized Dantzig #
  lambda.norm <- lambda+sqrt(2*log(N)/T0)*sig
  g <- rbind(diag(x=1,N,N),t(X0.scale)%*%(X0)/T0,-t(X0.scale)%*%(X0)/T0)
  h <- rbind(matrix(0,N,1),cbind(c(Y0X0.scale.mean-lambda.norm,-(Y0X0.scale.mean+lambda.norm))))
  tryCatch(betaHat.norm <-limSolve::lsei(A=EX1,B=EY1,E=e,F=f,G=g,H=h,type=2)$X,
           error = function(e) { betaHat.norm <<- NA})
  tryCatch(tauHat.norm <- EY1-EX1%*%betaHat.norm,
           error = function(e) { tauHat.norm <<-NA})
  
  # Multiplier Bootstrap #
  tuning.vec <- vector()
  for (m in 1:M) {
    u <- rnorm(T0,mean = 0, sd = 1)
    tuning.vec[m] <- sig*max(abs(apply(X0.scale*u,2,mean)))
  }
  c0 <- quantile(tuning.vec,c(0.9,0.95,0.99))
  lambda90 <- lambda+c0[1]
  lambda95 <- lambda+c0[2]
  lambda99 <- lambda+c0[3]
  tauHat90 <- tauHat95 <- tauHat99 <- NA
  g90 <- rbind(diag(x=1,N,N),t(X0.scale)%*%(X0)/T0,-t(X0.scale)%*%(X0)/T0)
  h90 <- rbind(matrix(0,N,1),cbind(c(Y0X0.scale.mean-lambda90,-(Y0X0.scale.mean+lambda90))))
  tryCatch(betaHat90 <-limSolve::lsei(A=EX1,B=EY1,E=e,F=f,G=g90,H=h90,type=2)$X,
           error = function(e) { betaHat90 <<- NA})
  tryCatch(tauHat90<- EY1-EX1%*%betaHat90,
           error = function(e) { tauHat90 <<-NA})
  g95 <- rbind(diag(x=1,N,N),t(X0.scale)%*%(X0)/T0,-t(X0.scale)%*%(X0)/T0)
  h95 <- rbind(matrix(0,N,1),cbind(c(Y0X0.scale.mean-lambda95,-(Y0X0.scale.mean+lambda95))))
  tryCatch(betaHat95 <-limSolve::lsei(A=EX1,B=EY1,E=e,F=f,G=g95,H=h95,type=2)$X,
           error = function(e) { betaHat95 <<- NA})
  tryCatch(tauHat95<- EY1-EX1%*%betaHat95,
           error = function(e) { tauHat95 <<-NA})
  g99 <- rbind(diag(x=1,N,N),t(X0.scale)%*%(X0)/T0,-t(X0.scale)%*%(X0)/T0)
  h99 <- rbind(matrix(0,N,1),cbind(c(Y0X0.scale.mean-lambda99,-(Y0X0.scale.mean+lambda99))))
  tryCatch(betaHat99 <-limSolve::lsei(A=EX1,B=EY1,E=e,F=f,G=g99,H=h99,type=2)$X,
           error = function(e) { betaHat99 <<- NA})
  tryCatch(tauHat99<- EY1-EX1%*%betaHat99,
           error = function(e) { tauHat99 <<-NA})
  # c0
  
  # if (tuning =="Normalize") {
  #   
  # } else if (tuning == "Sample") {
  #   
  #   
  # }
  # 
  # Increasing C #
  tauHat.c <- NA
  betaHat.c <- NA
  if (step>0) {
    C <- 0
    
    while (is.na(tauHat.c)){
      lambda <- lambda+C*sig*sqrt(log(N)/T0) # Dantzig
      ### Point Estimator ###
      g <- rbind(diag(x=1,N,N),X0.cov,-X0.cov)
      h <- rbind(matrix(0,N,1),cbind(c(Y0X0.mean-lambda,-(Y0X0.mean+lambda))))
      tryCatch(betaHat.c <-limSolve::lsei(A=EX1,B=EY1,E=e,F=f,G=g,H=h,type=2)$X,
               error = function(e) { betaHat.c <<- NA})
      tryCatch(tauHat.c <- EY1-EX1%*%betaHat.c,
               error = function(e) { tauHat.c <<-NA})
      C <- C+step
    }
  }
  
  
  
  
  ###################
  #### Inference ####
  ###################
  X0.cov.vec <- vecl(X0.cov) # vecl(\Sigma)
  vecl.X0t <- matrix(0,nrow = T0,ncol = length(X0.cov.vec))
  for (t in 1:T0) {
    vecl.X0t[t,] <- vecl(t(X0[t,,drop = F])%*%X0[t,,drop = F])
  }
  cov.vecl.X0 <- cov(vecl.X0t)
  # sqrt(ncol(cov.vecl.X0)/T0)
  dW <- max(gamma*max(cov.vecl.X0),1)*sqrt(log(ncol(cov.vecl.X0))/T0) # Some problem
  Sigma.cov <- cov.vecl.X0/T0+dW/T0*diag(ncol(cov.vecl.X0))
  
  
  Y0X0.cov <- cov(Y0X0.vec)
  dZ <- max(gamma*max(Y0X0.cov),1)*sqrt(log(ncol(Y0X0.cov))/T0)
  
 
    
  ### non-zero and Omega ###
  beta.mat <- matrix(NA,nrow = M,ncol = N)
  tau.vec <- vector()
  Int.mat <- matrix(NA,nrow = M,ncol = 2)
  # dY <- max(gamma*norm(VarY1,"I"),1)
  dX <- max(gamma*max(CovX1),1)*sqrt(log(ncol(CovX1))/T1)
  if (Inference) {
    for (m in 1:M) {
      skip_to_next <<- FALSE
      if (!cond) { # unconditional
        XX.repro <- MASS::mvrnorm(mu = X0.cov.vec,Sigma = Sigma.cov)
        # drop outlier
        if (max(abs(XX.repro-X0.cov.vec)/sqrt(diag(Sigma.cov)))>1.1*qnorm(1-0.01/length(XX.repro))) {
          next
        }
        tmp <- matrix(0,nrow = N,ncol = N)
        tmp[lower.tri(tmp,diag = T)] <- XX.repro
        # Symmetric
        for(l in 2:N) {
          for(k in 1:(l-1)) {
            tmp[k,l] = tmp[l,k]
          }
        }
        
        
        Diag.XX<-diag(eigen(tmp)$values)
        for(ind in 1:N){
          Diag.XX[ind,ind]<-max(Diag.XX[ind,ind],0.001)
        }
        ev <- eigen(tmp)
        XX.positive<-ev$vectors%*%Diag.XX%*%t(ev$vectors)
      } else{
        XX.positive <- X0.cov
      }
      XY.repro <- MASS::mvrnorm(mu = Y0X0.mean,
                                Sigma = Y0X0.cov/T0+dZ/T0*diag(ncol(Y0X0.cov)))
      if (max(abs(XY.repro-Y0X0.mean)/sqrt(diag(Y0X0.cov/T0)+dZ/T0*diag(ncol(Y0X0.cov))))>1.1*qnorm(1-0.01/length(XY.repro))) {
        next
      }
      EY1.m <- rnorm(1,mean = EY1, sd = sqrt((VarY1)/T1))
      if (max(abs(EY1.m-EY1)/sqrt((VarY1)/T1))>1.1*qnorm(1-0.01/length(EY1.m))) {
        next
      }
      if (cond) {
        EX1.m <- EX1
      } else{
        EX1.m <- MASS::mvrnorm(n=1,mu = EX1,Sigma = CovX1/T1+dX/T1*diag(N))
        if (max(abs(EX1.m-EX1)/sqrt(diag(CovX1/T1+dX/T1*diag(N))))>1.1*qnorm(1-0.01/length(EX1.m))) {
          next
        }
      }
      
      A <- (EX1.m)%*%t(EX1.m)
      if (CVXR) {
        
        betaHat <- CVXR::Variable(N)
        objective   <- CVXR::Minimize( CVXR::quad_form(betaHat,A) -2*sum(EY1.m*EX1.m*betaHat)) 
        constraints <- list(betaHat>=0,sum(betaHat)==1,
                            max(abs(XY.repro-XX.positive%*%betaHat))<=lambda)
        problem     <- CVXR::Problem(objective, constraints)
        suppressWarnings({
          solution <- CVXR::solve(problem) 
          tryCatch( beta.m <- as.vector(solution$getValue(betaHat)),
                    error = function(e) { skip_to_next <<- TRUE})
          if(skip_to_next) { next }
          if(any(is.na(beta.m))) {next}
        })
        
        
      } else{
        g <- rbind(diag(x=1,N,N),XX.positive,-XX.positive)
        h <- rbind(matrix(0,N,1),cbind(c(XY.repro-lambda,-(XY.repro+lambda))))
        tryCatch(beta.m <-limSolve::lsei(A=EX1.m,B=EY1.m,E=e,F=f,G=g,H=h,type=2)$X,
                 error = function(e) { skip_to_next <<- TRUE})
        
        if(skip_to_next) { next }
        
      }
      beta.mat[m,] <- beta.m
      tau.m <- EY1 - EX1%*%beta.m
      tau.vec[m] <- tau.m
      SE.m <- sd(Y1-X1%*%beta.m)
      Int.mat[m,] <- c(tau.m-qnorm(1-alpha/2)*SE.m/sqrt(T1),tau.m+qnorm(1-alpha/2)*SE.m/sqrt(T1))
      
    }
  }
 
  if (all(is.na(Int.mat))) {
    CI.tau = NA
  } else{
    suppressWarnings({
      uni = intervals::Intervals(Int.mat)
      CI.tau = as.matrix(intervals::interval_union(uni))
    })
  }
  return(list(betaHat.norm = betaHat.norm,tauHat.norm = tauHat.norm,
              betaHat90 = betaHat90,tauHat90 = tauHat90,
              betaHat95 = betaHat95,tauHat95 = tauHat95,
              betaHat99 = betaHat99,tauHat99 = tauHat99,
              betaHat.c = betaHat.c,tauHat.c = tauHat.c,
              CI.tau = CI.tau))
  
}


# 
# scsens <- function(Y0,Y1,X0,X1,lambda = NULL, intercept = F,tau.grid = NULL){
#   T0 <- length(Y0)
#   T1 <- length(Y1)
#   N <- ncol(X0)
#   ind.vec <- vector()
#   for (i in 1:length(tau.grid)) {
#     Y10 <- Y1-tau.grid[i]
#     beta.hat <- sc(Y10,X1,intercept = intercept)$w.hat
#     if (max(abs(1/T0*t(X0)%*%(Y0-X0%*%beta.hat)))<=lambda) {
#       ind.vec <- c(ind.vec,i)
#     }
#   }
#   tau.vec <- tau.grid[ind.vec]
#   return(tau.vec)
# }


# sc.repro <- function(Y1,Y0,T0,intercept = F, M = 1000, center = F,alpha.beta = 0.05,
#                      alpha2 = 0.05,method = "qreg",true.vec = NULL,sample_est = T,IQR = T,thres = NULL){
# 
#   Y1 <- cbind(Y1)
#   Y0 <- cbind(Y0)
#   T01 <- length(Y1)
#   T1 <- T01-T0
#   J <- dim(Y0)[2]
#   Y1.pre <- Y1[1:T0,,drop = FALSE]; Y1.post <- Y1[(T0+1):T01,drop = FALSE]
#   Y0.pre <- Y0[1:T0,,drop = FALSE];  Y0.post <- Y0[(T0+1):T01,,drop = FALSE]
#   if (center) {
#     Y1.pre <- Y1.pre - mean(Y1.pre)
#     Y0.pre <-  t(apply(Y0.pre,1,function(x){x - colMeans(Y0.pre)}))
#   }
# 
# 
#   CI.in.sample = matrix(nrow = M, ncol = 2*T1)
#   CI.tau = matrix(nrow = M, ncol = 2*T1)
#   beta.mat = beta.left  = beta.right = matrix(0,nrow = M,ncol = J)
#   CI_beta <- matrix(0,nrow = J, ncol = 2)
#   tau.mat <- matrix(0,nrow = M, ncol = T1)
#   beta.zero = rep(0,J)
#   sc.fit <- sc(Y1.pre, Y0.pre)
#   sc.resid <- sc.fit$u.hat
#   sc.beta <- sc.fit$w.hat
#   if (length(thres)==0) {
#     sigma.u <- sqrt(mean(sc.resid^2))
#     sigma.x <- apply(Y0.pre,2,sd)
#     # T0 <- nrow(Y0.pre)
#     c <- sigma.u/min(sigma.x)
#     thres <- c*sqrt(log(T0))/sqrt(T0)
#   }
# 
#   # idx <- which(sc.beta>0)
#   # Y0.pre.idx <- Y0.pre[,idx,drop = F]
#   # Y0.post.idx <- Y0.post[,idx,drop = F]
#   # res.fit <- scpi.out(sc.resid,Y0.pre.idx,Y0.post.idx,weight = sc.beta,
#   #                     method = method,alpha = alpha2,thres = thres)
#   # lb <- res.fit$lb
#   # ub <- res.fit$ub
#   # sig.vec <- vector()
#   res.mat <- matrix(0,nrow = M,ncol = T0)
#    lb.mat <- ub.mat <- matrix(0,nrow = M, ncol = T1)
#   beta.norm <- vector()
#   for (m in 1:M) {
#     u.repro <- rnorm(T0,sd = 1)
#     sc.fit.repro <- sc(Y1.pre, Y0.pre,F, u.repro)
#     sc.repro <- sc.fit.repro$w.hat
#     # check two norm between true beta and beta^{[m]}.
#     if (length(true.vec)>0) {
#       beta.norm[m] <- sum((sc.repro-true.vec)^2)
#     }
#     v.repro <- sc.repro>thres
#     Y0.repro <- Y0.pre[,v.repro,drop = F]
#     if (length(Y0.repro)==0) {
#       next
#     }
#     Y0.post.repro <- Y0.post[,v.repro,drop = F]
#     fit.repro <- lm(Y1.pre~0+Y0.repro) # no intercept
#     var.repro <- sandwich::vcovHC(fit.repro)/T0
#     sd.repro <- sqrt(diag(var.repro)) # Heteroscedastic error
#     beta.repro <- fit.repro$coefficients
#     # construct confidence interval for each beta.j by the distribution of u.repro
#     beta.mat[m,] <- sc.repro
#     var.in.sample <- diag(Y0.post.repro%*%var.repro%*%t(Y0.post.repro))
#     CI.in.sample[m,] <- as.vector(t(cbind( Y0.post.repro%*%beta.repro-qnorm(1-alpha.beta/2)*var.in.sample,
#                                            Y0.post.repro%*%beta.repro+qnorm(1-alpha.beta/2)*var.in.sample)))
# 
#     beta.left[m,which(v.repro)] <- beta.repro-qnorm(1-alpha.beta/2)*sd.repro
#     beta.right[m,which(v.repro)] <- beta.repro+qnorm(1-alpha.beta/2)*sd.repro
#     # sig.vec[m] <- sc.fit.repro$sigma
#     sc.resid <- sc.fit.repro$u.hat
#     res.mat[m,] <- sc.resid
#     if (sample_est) {
#       res.fit <- scpi.out(sc.resid,Y0.pre,Y0.post,weight = sc.repro,
#                           method = method,alpha = alpha2,IQR = IQR,thres = thres)
#       lb <- res.fit$lb
#       ub <- res.fit$ub
#     }
#     lb.mat[m,] <- lb
#     ub.mat[m,] <- ub
# 
#     # # out-of-sample bound by sigma^2
#     # lb <- -qnorm(1-alpha2/2)*sqrt(sig.repro)
#     # ub <- -lb
#     tau.t.repro <- Y1.post-Y0.post%*%sc.repro
#     # cat("sigma hat = ",sig.temp,", tau hat = ",tau.t.repro[1],",",tau.t.repro[2],"\n")
# 
#     # construct confidence interval for tau
# 
#     tau.mat[m,] <- tau.t.repro
#     CI.tau[m,]<- as.vector(t(cbind(tau.t.repro-ub,tau.t.repro-lb)))
#   }
# 
#   CI.in.lower <- apply(CI.in.sample[, c(TRUE,FALSE),drop = FALSE],2,function(x) {min(x,na.rm = T)})
#   CI.in.upper <- apply(CI.in.sample[,!c(TRUE,FALSE),drop = FALSE],2,function(x) {max(x,na.rm = T)})
#   CI_in_sample <- cbind(CI.in.lower,CI.in.upper)
# 
#   CI.lower <- apply(CI.tau[, c(TRUE,FALSE),drop = FALSE],2,function(x) {min(x,na.rm = T)})
#   CI.upper <- apply(CI.tau[,!c(TRUE,FALSE),drop = FALSE],2,function(x) {max(x,na.rm = T)})
#   CI_tau <- cbind(CI.lower,CI.upper)
# 
#   # check whether beta hat is close to true beta
#   beta_min <- min(beta.norm)
# 
# 
#   for (j in 1:J) {
#     idx.j <- which(beta.right[,j]>0)
#     if (length(idx.j)<length(beta.right[,j])) {
#       beta.zero[j] <- 1
#     }
#     if (length(idx.j)==0) {
#       next
#     }
#     CI.j <- cbind(beta.left[idx.j,j],beta.right[idx.j,j])
#     CI_beta[j,] <- c(max(min(CI.j[,1]),0),min(max(CI.j[,2]),1))
# 
#   }
#   return(list(CI_tau = CI_tau, CI_in_sample = CI_in_sample,
#               CI.in.sample=CI.in.sample,
#               CI_beta = CI_beta,CI.tau = CI.tau,tau.mat = tau.mat,
#               lb.mat = lb.mat,ub.mat = ub.mat,beta.mat = beta.mat,res.mat = res.mat,
#               beta.left = beta.left,beta.right = beta.right,
#               beta_min = beta_min))
#   # return(list(CI_tau = CI_tau, CI_beta = CI_beta,beta.zero = beta.zero))
# }



scpi.out <- function(res,X.pre,X.post,weight = NULL,thres = NULL,
                     method = "qreg",alpha = 0.05,IQR = T,intercept = T){
  if (length(thres)==0) {
    sigma.u <- sqrt(mean(res^2))
    sigma.x <- apply(X.pre,2,sd)
    T0 <- nrow(X.pre)
    c <- sigma.u/min(sigma.x)
    thres <- c*sqrt(log(T0))/sqrt(T0)
  }
  
  if (length(weight)>0) {
    idx <- which(weight>thres)
    # cat("idx = ",idx,"\n")
    # print(X.pre)
    X.pre <- X.pre[,idx ,drop = FALSE]
    X.post <- X.post[,idx ,drop = FALSE]
  }
  if (intercept) {
    X.pre <- cbind(X.pre,rep(1,nrow(X.pre)))
    X.post <- cbind(X.post,rep(1,nrow(X.post)))
  }
  if (method == "gaussian") {
    fit      <- predict(y = res, x = X.pre,
                        eval = rbind(X.post,X.pre), type = "lm")
    e.mean   <- fit[1:nrow(X.post)]
    res.fit  <- fit[-(1:nrow(X.post))]
    
    var.pred <- predict(y = log((res-res.fit)^2), x = X.pre, 
                        eval = rbind(X.post,X.pre), type = "lm")
    var.pred <- var.pred[1:nrow(X.post)]
    
    e.sig2   <- exp(var.pred)
    e.sig <- sqrt(e.sig2)
    if (IQR) {
      q.pred <- predict(y = res- res.fit, x = X.pre, 
                        eval = rbind(X.post,X.pre), type = "qreg", tau = c(0.25, 0.75))
      q3.pred <- q.pred[1:nrow(X.post), 2]
      q1.pred <- q.pred[1:nrow(X.post), 1]
      IQ.pred <- q3.pred - q1.pred
      IQ.pred <- abs(IQ.pred)
      e.sig <- apply(cbind(sqrt(e.sig2), IQ.pred/1.34), 1, min)
    }
    eps <- sqrt(-log(alpha/2)*2)*e.sig
    lb <- e.mean - eps
    ub <- e.mean + eps
  } else if (method == "qreg"){
    e.pred <- predict(y = res, x = X.pre, 
                      eval = X.post, type = "qreg", tau = c(alpha/2, 1-alpha/2))
    e.pred.lb <- e.pred[,1]
    e.pred.ub <- e.pred[,2]
    lb <- e.pred.lb
    ub <- e.pred.ub
  }
  return(list(lb = lb,ub = ub))
}

predict <- function(y, x, eval, type="lm", tau=NULL, verbose = FALSE) {
  
  if (type == "lm") {
    betahat <-.lm.fit(x, y)$coeff
    pred <- eval %*% betahat
    
  } else if (type == "qreg") {
    
    tryCatch(
      {
        betahat <- Qtools::rrq(y~x-1, tau=tau)$coefficients
      },
      
      warning = function(war) {
        message("Warning produced when estimating moments of the out-of-sample residuals with quantile regressions.")
        war$call <- NULL
        if (verbose) warning(war)
      }
    )
    # quantreg::rq
    betahat <- suppressWarnings(Qtools::rrq(y~x-1, tau=tau)$coefficients)
    pred <- eval %*% betahat
  }
  
  
  return(pred)
}



