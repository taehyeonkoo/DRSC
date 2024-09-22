####################
# M3 does not hold #
####################
rm(list=ls())
library(ggplot2)
source("~/Dropbox/Taehyeon/Synthetic Control/R-code/helpers.R")
T0 <- 50
T1 <- 50
N <- 16
Sigma0 <- diag(1,N)#0.05*diag(1,N)+0.95*matrix(1,N,N)#A1gen(0.1,N)

# Sigma0[1,N] <- Sigma0[N,1] <- 0.99
# Sigma0[2,(N-1)] <- Sigma0[(N-1),2] <- 0.99
# Sigma0[3,(N-2)] <- Sigma0[(N-2),3] <- 0.99


eigen(Sigma0)$values

Sigma1 <- diag(1,N)
# Sigma1[1,N] <- Sigma1[N,1] <- -0.99
# Sigma1[2,(N-1)] <- Sigma1[(N-1),2] <- -0.99
# Sigma1[3,(N-2)] <- Sigma1[(N-2),3] <- -0.99

eigen(Sigma1)$values


mu_X <- rep(10,N)+rnorm(N,mean = 0, sd = 0.05) 
# mu_X[N] <- mu_X[N]
mu1 <- 1:N+rnorm(N,mean = 0, sd = 0.2)

SC.true <- c(1/3,1/3,1/3,rep(0,N-3))


# SC.true <- c(rep(0,N/2-2),0.25,0.25,0.25,0.25,rep(0,N/2-2))
# SC.true <- rep(1/N,N)
# SC.true <- c((1-c0)/2*((N/2-1):1)/sum((N/2-1):1),c0/2,c0/2,(1-c0)/2*((N/2-1):1)/sum((N/2-1):1))
c <- 1
SC.post <- SC.true


sum(SC.post)
# Subset of simplex
A <- mu_X%*%t(mu_X)+matrix(1,N,N)+diag(N)
lambda <-0


gg.list <- list()
beta.list <- beta.list2 <- beta.list3 <- beta.list4 <- beta.list5 <-beta.list6 <- list()
bias.mat <- sd.mat <- mse.mat <- rmse.mat <- matrix(nrow = 21,ncol = 6)

colnames(mse.mat) <- colnames(rmse.mat) <- c(
  "SC","DRSC-norm","DRSC-90","DRSC-95","DRSC-99","DRSC-C"
)
nsims <- 500
for (i in 1:21) {
  # i <- 1
  tau <- 5*(-1+(i-1)/10)
  # tau <- -2
  beta.star <- SC.post
  tau.star <- c(tau+mu1%*%(SC.post-beta.star))
  sc.vec <- tau.vec <- tau.vec2 <-tau.vec3 <-tau.vec4 <-vector()
  drsc.mat <-drsc.mat2 <-drsc.mat3 <-drsc.mat4 <- drsc.mat5 <-sc.mat <- matrix(nrow = nsims, ncol = N)
  dt <- data.frame()
  lambda <- 0
  for (m in 1:nsims) {
    F1 <- rnorm(T0+T1)
    F2 <- rnorm(T0+T1)
    eps_y <- rnorm(T0+T1,sd = sqrt(1))
    eps_X <- MASS::mvrnorm(T0,mu = mu_X,Sigma = Sigma0)
    # X <- cbind(rep(1,T0+T1))%*%rbind(mu_X)+cbind(F1)%*%rbind(rep(c,N))+
    #   cbind(F2)%*%rbind(mu_X)+eps_X
    X0 <- eps_X
    X1 <- MASS::mvrnorm(T1,mu = mu1,Sigma = Sigma1)
    # X <- matrix(rnorm((T0+T1)*N,sd = 1),nrow = T0+T1,ncol = N)
    # X0 <- X[1:T0,];X1 <- X[-(1:T0),]
    Y0 <- as.vector(X0%*%SC.true+eps_y[1:T0])
    Y1 <- as.vector(X1%*%SC.true+eps_y[-(1:T0)])
    # X <- sqrt(3)*matrix(runif((T0+T1)*N),nrow = T0+T1,ncol = N)
    # X <- eps_X
    # Y <- X%*%SC.true+eps_y
    # Y0 <- as.vector(X[1:T0,]%*%SC.true+eps_y[1:T0])
    # Y1 <- as.vector(X[-(1:T0),]%*%SC.post+eps_y[-(1:T0)])
    # Y0 <- Y[1:T0];Y1 <- Y[-(1:T0)]
    # tau.t <- rnorm(T1,tau,sd = 0.05)
    Y1 <- Y1+tau
    SC.beta <- sc(Y0,X0)$w.hat
    tau.SC <- mean(Y1-X1%*%SC.beta)
    sc.vec[m] <- tau.SC
    sc.mat[m,] <- SC.beta
    DRSC.fit <- DRSC(Y0,Y1,X0,X1,lambda = lambda, M = 1000,step = 1e-7)
    drsc.mat[m,] <- DRSC.fit$betaHat.norm
    drsc.mat2[m,] <- DRSC.fit$betaHat90
    drsc.mat3[m,] <- DRSC.fit$betaHat95
    drsc.mat4[m,] <- DRSC.fit$betaHat99
    drsc.mat5[m,] <- DRSC.fit$betaHat.c
    
    # tau.DRSC <- DRSC.fit$tauHat
    # tau.vec[m]<- tau.DRSC
    dt <- rbind(dt,c(lambda,tau.star,T0,T1,tau.SC,"SC"),
                c(lambda,tau.star,T0,T1,DRSC.fit$tauHat.norm,"DRSC-norm"),
                c(lambda,tau.star,T0,T1,DRSC.fit$tauHat90,"DRSC-90"),
                c(lambda,tau.star,T0,T1,DRSC.fit$tauHat95,"DRSC-95"),
                c(lambda,tau.star,T0,T1,DRSC.fit$tauHat99,"DRSC-99"),
                c(lambda,tau.star,T0,T1,DRSC.fit$tauHat.c,"DRSC-C")
                )
  }
  beta.list[[i]] <- sc.mat
  beta.list2[[i]] <- drsc.mat
  beta.list3[[i]] <- drsc.mat2
  beta.list4[[i]] <- drsc.mat3
  beta.list5[[i]] <- drsc.mat4
  beta.list6[[i]] <- drsc.mat5
  
  
  
  colnames(dt) <- c(
    "lambda","tau","T0","T1","tauHat","Method"
  )
  head(dt)
  # dt <- na.omit(dt)
  dt$tauHat <- as.numeric(dt$tauHat)
  dt$Method <- as.factor(dt$Method)
  gg.list[[i]] <- ggplot(dt,aes(x=Method, y=tauHat, fill = Method)) +
    geom_violin(alpha = 0.5)+
    geom_hline(yintercept = tau,colour = 'blue',linetype="dashed")+
    geom_hline(yintercept = tau.star,colour ='red',linetype="dashed")+
    ylab(latex2exp::TeX("$\\hat{\\tau}$"))
  mse.mat[i,] <- c(mean((dt[dt$Method == "SC",]$tauHat-tau)^2),
  mean((dt[dt$Method == "DRSC-norm",]$tauHat-tau.star)^2,na.rm = T),
  mean((dt[dt$Method == "DRSC-90",]$tauHat-tau)^2,na.rm = T),
  mean((dt[dt$Method == "DRSC-95",]$tauHat-tau)^2,na.rm = T),
  mean((dt[dt$Method == "DRSC-99",]$tauHat-tau)^2,na.rm = T),
  mean((dt[dt$Method == "DRSC-C",]$tauHat-tau)^2,na.rm = T))
  bias.mat[i,] <- c(mean(dt[dt$Method == "SC",]$tauHat,na.rm = T)-tau,
                    mean(dt[dt$Method == "DRSC-norm",]$tauHat,na.rm = T)-tau,
                    mean(dt[dt$Method == "DRSC-90",]$tauHat,na.rm = T)-tau,
                    mean(dt[dt$Method == "DRSC-95",]$tauHat,na.rm = T)-tau,
                    mean(dt[dt$Method == "DRSC-99",]$tauHat,na.rm = T)-tau,
                    mean(dt[dt$Method == "DRSC-C",]$tauHat,na.rm = T)-tau)
  sd.mat[i,] <- c(sd(dt[dt$Method == "SC",]$tauHat,na.rm = T),
                    sd(dt[dt$Method == "DRSC-norm",]$tauHat,na.rm = T),
                    sd(dt[dt$Method == "DRSC-90",]$tauHat,na.rm = T),
                    sd(dt[dt$Method == "DRSC-95",]$tauHat,na.rm = T),
                    sd(dt[dt$Method == "DRSC-99",]$tauHat,na.rm = T),
                    sd(dt[dt$Method == "DRSC-C",]$tauHat,na.rm = T))
}

df <- data.frame(beta2=beta.list[[1]][,1],beta3 = beta.list[[1]][,2],
                 beta4=beta.list[[1]][,3])
df2 <- data.frame(beta2=beta.list2[[1]][,1],beta3 = beta.list2[[1]][,2],
                 beta4=beta.list2[[1]][,3])
df3 <- data.frame(beta2=beta.list3[[1]][,1],beta3 = beta.list3[[1]][,2],
                 beta4=beta.list3[[1]][,3])

ggplot(df, aes(x=beta2)) + 
  geom_histogram(aes(y=..density..), colour="black", fill="white")+
  geom_density(alpha=.2, fill="#FF6666") +
  geom_vline(aes(xintercept=1/3), color="blue",
             linetype="dashed")

beta.list[[1]][,c(1,2,3,N-2,N-1,N)]
beta.list2[[1]][,c(1,2,3,N-2,N-1,N)]
head(dt)

gg.list[[21]]

mse.mat
head(mse.dt)
rmse.dt <- data.frame(sqrt(mse.mat),tau = 5*(-1+(1:21-1)/10))

ggplot(rmse.dt,aes(x = tau))+
  geom_point(aes(y = SC,colour = 'SC'))+
  geom_line(aes(y = SC,colour = 'SC'))+
  geom_point(aes(y = DRSC.norm,colour = 'DRSC.norm'))+
  geom_line(aes(y = DRSC.norm,colour = 'DRSC.norm'))+
  geom_point(aes(y = DRSC.90,colour = 'DRSC.90'))+
  geom_line(aes(y = DRSC.90,colour = 'DRSC.90'))+
  geom_point(aes(y = DRSC.95,colour = 'DRSC.95'))+
  geom_line(aes(y = DRSC.95,colour = 'DRSC.95'))+
  geom_point(aes(y = DRSC.99,colour = 'DRSC.99'))+
  geom_line(aes(y = DRSC.99,colour = 'DRSC.99'))+
  geom_point(aes(y = DRSC.C,colour = 'DRSC.C'))+
  geom_line(aes(y = DRSC.C,colour = 'DRSC.C'))+
  xlab(latex2exp::TeX("$\\tau$")) +
  ylab("RMSE")+
  labs(colour='Methods')+
  ggtitle(latex2exp::TeX(paste0( "$\\lambda=$",lambda)))+
  theme(plot.title = element_text(hjust = 0.5))

bias.dt <- data.frame(bias.mat,tau = 5*(-1+(1:21-1)/10))
abs.bias.dt <- data.frame(abs(bias.mat),tau = 5*(-1+(1:21-1)/10))
colnames(bias.dt) <- c("SC","DRSC.norm","DRSC.90","DRSC.95","DRSC.99","DRSC.C","tau")
colnames(abs.bias.dt) <- c("SC","DRSC.norm","DRSC.90","DRSC.95","DRSC.99","DRSC.C","tau")

ggplot(bias.dt,aes(x = tau))+
  geom_point(aes(y = SC,colour = 'SC'))+
  geom_line(aes(y = SC,colour = 'SC'))+
  geom_point(aes(y = DRSC.norm,colour = 'DRSC.norm'))+
  geom_line(aes(y = DRSC.norm,colour = 'DRSC.norm'))+
  geom_point(aes(y = DRSC.90,colour = 'DRSC.90'))+
  geom_line(aes(y = DRSC.90,colour = 'DRSC.90'))+
  geom_point(aes(y = DRSC.95,colour = 'DRSC.95'))+
  geom_line(aes(y = DRSC.95,colour = 'DRSC.95'))+
  geom_point(aes(y = DRSC.99,colour = 'DRSC.99'))+
  geom_line(aes(y = DRSC.99,colour = 'DRSC.99'))+
  geom_point(aes(y = DRSC.C,colour = 'DRSC.C'))+
  geom_line(aes(y = DRSC.C,colour = 'DRSC.C'))+
  xlab(latex2exp::TeX("$\\tau$")) +
  ylab("Bias")+
  labs(colour='Methods')+
  ggtitle("Bias")+
  theme(plot.title = element_text(hjust = 0.5))

ggplot(abs.bias.dt,aes(x = tau))+
  geom_point(aes(y = SC,colour = 'SC'))+
  geom_line(aes(y = SC,colour = 'SC'))+
  geom_point(aes(y = DRSC.norm,colour = 'DRSC.norm'))+
  geom_line(aes(y = DRSC.norm,colour = 'DRSC.norm'))+
  geom_point(aes(y = DRSC.90,colour = 'DRSC.90'))+
  geom_line(aes(y = DRSC.90,colour = 'DRSC.90'))+
  geom_point(aes(y = DRSC.95,colour = 'DRSC.95'))+
  geom_line(aes(y = DRSC.95,colour = 'DRSC.95'))+
  geom_point(aes(y = DRSC.99,colour = 'DRSC.99'))+
  geom_line(aes(y = DRSC.99,colour = 'DRSC.99'))+
  geom_point(aes(y = DRSC.C,colour = 'DRSC.C'))+
  geom_line(aes(y = DRSC.C,colour = 'DRSC.C'))+
  xlab(latex2exp::TeX("$\\tau$")) +
  ylab("Absolute value of Bias")+
  labs(colour='Methods')+
  ggtitle("Absolute value of Bias")+
  theme(plot.title = element_text(hjust = 0.5))

sd.dt <- data.frame(sd.mat,tau = 5*(-1+(1:21-1)/10))
colnames(sd.dt) <- c("SC","DRSC.norm","DRSC.90","DRSC.95","DRSC.99","DRSC.C","tau")

ggplot(sd.dt,aes(x = tau))+
  geom_point(aes(y = SC,colour = 'SC'))+
  geom_line(aes(y = SC,colour = 'SC'))+
  geom_point(aes(y = DRSC.norm,colour = 'DRSC.norm'))+
  geom_line(aes(y = DRSC.norm,colour = 'DRSC.norm'))+
  geom_point(aes(y = DRSC.90,colour = 'DRSC.90'))+
  geom_line(aes(y = DRSC.90,colour = 'DRSC.90'))+
  geom_point(aes(y = DRSC.95,colour = 'DRSC.95'))+
  geom_line(aes(y = DRSC.95,colour = 'DRSC.95'))+
  geom_point(aes(y = DRSC.99,colour = 'DRSC.99'))+
  geom_line(aes(y = DRSC.99,colour = 'DRSC.99'))+
  geom_point(aes(y = DRSC.C,colour = 'DRSC.C'))+
  geom_line(aes(y = DRSC.C,colour = 'DRSC.C'))+
  xlab(latex2exp::TeX("$\\tau$")) +
  ylab("Standard error")+
  labs(colour='Methods')+
  ggtitle("Standard error")+
  theme(plot.title = element_text(hjust = 0.5))

gg.list[[19]]
gg.list[[21]]
gg.list[[10]]

save.image("C2.Rdata")

####################
# M2 does not hold #
####################
library(ggplot2)
source("~/Dropbox/Taehyeon/Synthetic Control/R-code/helpers.R")
T0 <- 31
T1 <- 13
N <- 16
# sig <- 10
Sigma0 <- A1gen(0.9999,N)
Sigma1 <- A1gen(0,N)*4
mu_X <- (1:N)/N
A <- mu_X%*%t(mu_X)+Sigma0


SC.true <- c(rep(0,N-4),0.1,0.3,0.3,0.3)
# SC.true <- c(rep(0,N/2-2),0.25,0.25,0.25,0.25,rep(0,N/2-2))
# SC.true <- rep(1/N,N)+rep(c(0.02,-0.02),N/2)
# SC.true <- c((1-c0)/2*((N/2-1):1)/sum((N/2-1):1),c0/2,c0/2,(1-c0)/2*((N/2-1):1)/sum((N/2-1):1))
c <- 0.44*2*2
SC.post <- SC.true + c*c(rep(0,N-6),0.1,0.1,0.1,-0.1,-0.1,-0.1)
sum(SC.post)
min(SC.post)
round(max(abs(A%*%(SC.true-SC.post))),3)

# Subset of simplex
lambda <-0.1
gg.list <- list()
mse.mat <- rmse.mat <- matrix(nrow = 21,ncol =3)

colnames(mse.mat) <- colnames(rmse.mat) <- c(
  "SC","DRSC","DRSC_star"
)
c <- 1

for (i in 1:21) {
  tau <- -1+(i-1)/10
  # tau <- -1
  
  g <- rbind(diag(x=1,N,N),A,-A)
  
  
  h <- rbind(matrix(0,N,1),cbind(c(A%*%SC.true-lambda,-(A%*%SC.true+lambda))))
  if (lambda!=0) {
    beta.star <-limSolve::lsei(A=mu_X,B=tau+mu_X%*%SC.post,E=cbind(matrix(1,1,N)),F=1,
                               G=g,H=h,type=2)$X
  } else{
    beta.star <- SC.post
  }
  
  tau.star <- c(tau+mu_X%*%(SC.post-beta.star))
  sc.vec <- tau.vec <- vector()
  sc.mat <- matrix(nrow = nsims, ncol = N)
  dt <- data.frame()
  nsims <- 500
  for (m in 1:nsims) {
    F1 <- rnorm(T0+T1)
    F2 <- rnorm(T0+T1)
    eps_y <- rnorm(T0+T1,sd = sqrt(0.5))
    eps_X <- MASS::mvrnorm(T0,mu = mu_X,Sigma = Sigma0)
    # X <- cbind(rep(1,T0+T1))%*%rbind(mu_X)+cbind(F1)%*%rbind(rep(c,N))+
    #   cbind(F2)%*%rbind(mu_X)+eps_X
    X0 <- eps_X
    X1 <- MASS::mvrnorm(T1,mu = mu_X,Sigma = Sigma1)
    # X <- matrix(rnorm((T0+T1)*N,sd = 1),nrow = T0+T1,ncol = N)
    # X0 <- X[1:T0,];X1 <- X[-(1:T0),]
    Y0 <- as.vector(X0%*%SC.true+eps_y[1:T0])
    Y1 <- as.vector(X1%*%SC.true+eps_y[-(1:T0)])
    # X <- sqrt(3)*matrix(runif((T0+T1)*N),nrow = T0+T1,ncol = N)
    # X <- eps_X
    # Y <- X%*%SC.true+eps_y
    # Y0 <- as.vector(X[1:T0,]%*%SC.true+eps_y[1:T0])
    # Y1 <- as.vector(X[-(1:T0),]%*%SC.post+eps_y[-(1:T0)])
    # Y0 <- Y[1:T0];Y1 <- Y[-(1:T0)]
    # tau.t <- rnorm(T1,tau,sd = 0.05)
    Y1 <- Y1+tau
    SC.beta <- sc(Y0,X0)$w.hat
    tau.SC <- mean(Y1-X1%*%SC.beta)
    sc.vec[m] <- tau.SC
    sc.mat[m,] <- SC.beta
    tau.DRSC <- DRSC(Y0,Y1,X0,X1,lambda = lambda,step=1e-7)$tauHat
    tau.vec[m]<- tau.DRSC
    dt <- rbind(dt,c(lambda,tau,tau.star,T0,T1,tau.SC,"SC"),c(lambda,tau,tau.star,T0,T1,tau.DRSC,"DRSC"))
  }
  colnames(dt) <- c(
    "lambda","tau","tau.star","T0","T1","tauHat","Method"
  )
  dt$tauHat <- as.numeric(dt$tauHat)
  dt$Method <- as.factor(dt$Method)
  gg.list[[i]] <- ggplot(dt,aes(x=Method, y=tauHat, fill = Method)) +
    geom_violin(alpha = 0.5)+
    geom_hline(yintercept = tau,colour = 'blue',linetype="dashed")+
    geom_hline(yintercept = tau.star,colour ='red',linetype="dashed")+
    ggtitle(paste0("tau = ",round(tau,3),", lambda = ",round(lambda,3),", T0 = ",T0,", T1 = ",T1))
  
  mse.mat[i,1] <- mean((dt[dt$Method =='SC' ,'tauHat']-tau)^2)
  mse.mat[i,2] <-mean((dt[dt$Method =='DRSC' ,'tauHat']-tau)^2)
  mse.mat[i,3] <-mean((dt[dt$Method =='DRSC' ,'tauHat']-tau.star)^2)
  rmse.mat[i,1] <- sqrt(mse.mat[i,1])
  rmse.mat[i,2] <- sqrt(rmse.mat[i,2])
}



mse.dt <- data.frame(mse.mat,tau = -1+(1:21-1)/10)

ggplot(mse.dt,aes(x = tau))+
  geom_point(aes(y = SC,colour = 'SC'))+
  geom_line(aes(y = SC,colour = 'SC'))+
  geom_point(aes(y = DRSC,colour = 'DRSC'))+
  geom_line(aes(y = DRSC,colour = 'DRSC'))+
  geom_point(aes(y = DRSC_star,colour = 'DRSC_star'))+
  geom_line(aes(y = DRSC_star,colour = 'DRSC_star'))+
  xlab(latex2exp::TeX("$\\tau$")) +
  ylab("MSE")+
  labs(colour='Methods')+
  ggtitle(latex2exp::TeX(paste0("$\\lambda = $",lambda)))+
  theme(plot.title = element_text(hjust = 0.5))

gg.list[[20]]
gg.list[[11]]
gg.list[[5]]

save.image("M3.Rdata")
