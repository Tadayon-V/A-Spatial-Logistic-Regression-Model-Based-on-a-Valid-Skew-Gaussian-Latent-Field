rm(list = ls())

# Packages
library(mvtnorm)
library(truncnorm)
library(TruncatedNormal)
library(geoR)
library(sp)
#install.packages("optimx")
library(optimx)

set.seed(1400)

n=100; m=1
L<-100 # for E7t
psi=3.5; delta=0.9; beta0=1.5; beta1=0.8; gamma=0.5; tau2=0.4
teta<-rnorm(n,0,1)      
X<-rnorm(n,0,1); epsilon<-rnorm(n,0,sqrt(tau2))

coor1 <- c(runif(n,0,10))
coor2 <- c(runif(n,0,10))
coordinates <- cbind(coor1,coor2)
class(coordinates); head(coordinates); dim(coordinates)
plot(coordinates)

NORMMH <- function(coordinate){
  NMH <- as.matrix(dist(coordinates))
  return(NMH)
}
normmatrix <- NORMMH(coordinates)
class(normmatrix); dim(normmatrix)

expocorr <- function(normmatrix, psi){
  M <- exp(-(normmatrix/psi))
  c <- c(min(M[lower.tri(M)][order(-row(M)[lower.tri(row(M))])]),
         max(M[lower.tri(M)][order(-row(M)[lower.tri(row(M))])]))
  return(c) 
}

z.s <- function(delta, psi, beta0, beta1, gamma, tau2){
  V.tranc <- rtruncnorm(n=n, a=0, b=Inf, mean=0, sd=1)
  #mean(V.tranc)
  #hist(V.tranc)
  mu.v <- delta*(V.tranc - sqrt(2/pi)*c(rep(1,n)))
  #mean(mu.v)
  #hist(mu.v)
  W.given.v <- grf(n=n, cov.model="exponential", cov.pars=c(1,psi),
                   grid=coordinates,
                   mean=mu.v
  )
  #mean(W.given.v$data)
  #hist(W.given.v$data)
  y.s <- beta0 + beta1*X + 
    gamma*W.given.v$data + epsilon 
  #mean(y.s)
  #hist(y.s)
  p.s <- exp(y.s)/(1 + exp(y.s))
  #mean(p.s)
  #hist(p.s)
  z.generates <- rbinom(n=n, size=1, prob=p.s)
  return(list(z.generates, W.given.v$data))
}
resu <- z.s(delta, psi, beta0, beta1, gamma, tau2)
Z.s <- resu[[1]]
W.skew <- resu[[2]]

par(mfrow=c(1,2))
hist(Z.s, prob=TRUE, col="grey")
lines(density(Z.s), col="blue", lwd=2)
hist(W.skew, prob=TRUE, col="grey")
lines(density(W.skew), col="blue", lwd=2)

K<-function(x) 1/(1+exp(-x))
Landa<-function(teta) (K(teta)-0.5)/(2*teta)

mu.v<- function(delta){
  V.tranc <- rtruncnorm(n=n, a=0, b=Inf, mean=0, sd=1)
  mv<- delta*(V.tranc - sqrt(2/pi)*c(rep(1,n)))
  return(mv)
}
# mu.v(delta)

sigmaepsilon.wz<-function(teta,tau2){
  d1<-(1/tau2)+2*Landa(teta) 
  sigma<-diag(1/d1)
  return(sigma) 
}
# sigmaepsilon.wz(teta,tau2)

miuepsilon.wz<-function(beta0,beta1,gamma,teta,tau2){
  z<-Z.s
  d1<-z-2*Landa(teta)*(beta0+beta1*X+gamma*W.skew)-0.5
  d2<-(1/tau2)+2*Landa(teta) 
  miu<-d1/d2
  return(miu)
  
}
# miuepsilon.wz(beta0,beta1,gamma,teta,tau2)

repsilon.wz<-function(m,beta0,beta1,gamma,teta,tau2){
  MIU<-miuepsilon.wz(beta0,beta1,gamma,teta,tau2)
  SIGMA<-sigmaepsilon.wz(teta,tau2)
  Results<-rmvnorm(m,MIU,SIGMA)
  return(Results)
}
# repsilon.wz(m=1,beta0,beta1,gamma,teta,tau2)

sigmaepsiloni.wz<-function(teta,tau2,i){
  tetai<-teta[i]
  sigma<-(1/tau2)+2*Landa(tetai) 
  return(1/sigma)
}
# sigmaepsiloni.wz(teta,tau2,i=2)

miuepsiloni.wz<-function(beta0,beta1,gamma,teta,tau2,i){
  z<-Z.s
  zi<-z[i];tetai<-teta[i];wi<-W.skew[i];Xi<-X[i]
  d1<-zi-2*Landa(tetai)*(beta0+beta1*Xi+gamma*wi)-0.5
  d2<-sigmaepsiloni.wz(teta,tau2,i) 
  miu<-d1/d2
  return(miu)
}
# miuepsiloni.wz(beta0,beta1,gamma,teta,tau2,i=2)

repsiloni.wz<-function(n,beta0,beta1,gamma,teta,tau2,i){
  MIU<-miuepsiloni.wz(beta0,beta1,gamma,teta,tau2,i)
  SIGMA<-sigmaepsiloni.wz(teta,tau2,i)
  Results<-rnorm(n,MIU,SIGMA)
  return(Results)
}
# repsiloni.wz(n,beta0,beta1,gamma,teta,tau2,i=2)

sigmaw.epsilonzold<-function(psi,gamma,teta){
  H <- exp(-(normmatrix/psi))
  landateta <- Landa(teta)
  D <- 2*(gamma^2)*diag(landateta)
  sigma <- solve(solve(H)+D)
  return(sigma)
}
# sigmaw.epsilonzold(psi,gamma,teta)

sigmaw.epsilonz<-function(psi,gamma,teta){
  H <- exp(-(normmatrix/psi))
  landateta <- Landa(teta)
  D <- 2*(gamma^2)*diag(landateta)
  invH<-chol2inv(chol(H))
  sigma <- chol2inv(chol(invH+D))
  return(sigma)
}
# sigmaw.epsilonz(psi,gamma,teta)

miuw.epsilonzold<-function(psi,beta0,beta1,gamma,teta){
  z<-Z.s
  H <- exp(-(normmatrix/psi))
  landateta<-Landa(teta)
  sigmawz<-sigmaw.epsilonz(psi,gamma,teta)
  miuv<-mu.v(delta)
  c1<-rep(0,n)
  for(i in 1:n){
    c1[i]<-gamma*(1+4*landateta[i]*(beta0+beta1*X[i]+epsilon[i])-2*z[i])}
  d1<- solve(H)%*%miuv-c1/2 
  miu<-sigmawz%*%d1
  return(miu)
}
# miuw.epsilonzold(psi,beta0,beta1,gamma,teta)

miuw.epsilonz<-function(psi,beta0,beta1,gamma,teta){
  z<-Z.s
  H <- exp(-(normmatrix/psi))
  landateta<-Landa(teta)
  sigmawz<-sigmaw.epsilonz(psi,gamma,teta)
  miuv<-mu.v(delta)
  c1<-rep(0,n)
  for(i in 1:n){
    c1[i]<-gamma*(1+4*landateta[i]*(beta0+beta1*X[i]+epsilon[i])-2*z[i])}
  invH<-chol2inv(chol(H))
  d1<- invH%*%miuv-c1/2 
  miu<-sigmawz%*%d1
  return(miu)
}
# miuw.epsilonz(psi,beta0,beta1,gamma,teta)

rw.epsilonzold<-function(m,psi,beta0,beta1,gamma,teta){
  MIU<-miuw.epsilonzold(psi,beta0,beta1,gamma,teta)
  SIGMA<-sigmaw.epsilonzold(psi,gamma,teta)
  Results<-rmvnorm(m,mean=MIU,sigma=SIGMA)
  return(Results) 
}
# rw.epsilonzold(m=1,psi,beta0,beta1,gamma,teta)

rw.epsilonz<-function(m,psi,beta0,beta1,gamma,teta){
  MIU<-miuw.epsilonz(psi,beta0,beta1,gamma,teta)
  SIGMA<-sigmaw.epsilonz(psi,gamma,teta)
  Results<-rmvnorm(m,mean=MIU,sigma=SIGMA)
  return(Results) 
}
# rw.epsilonz(m=1,psi,beta0,beta1,gamma,teta)

sigma.v.givenold <- function(n,delta,psi){
  H <- exp(-(normmatrix/psi))
  d1 <- (delta^2)*solve(H)
  d2 <- d1 + diag(x=1,nrow=n,ncol=n)
  return(solve(H))
}
# sigma.v.givenold(n,delta,psi)

sigma.v.given <- function(n,delta,psi){
  H <- exp(-(normmatrix/psi))
  d1 <- (delta^2)*chol2inv(chol(H))
  d2 <- d1 + diag(n)
  return(chol2inv(chol(d2)))
}
# sigma.v.given(n,delta,psi)

mu.v.givenold <- function(n,delta,psi,W.skew){
  H <- exp(-(normmatrix/psi))
  d1 <- delta*solve(H)%*%W.skew
  d2 <- sqrt(2/pi)*solve(H)%*%c(rep(1,n))
  d3 <- sigma.v.given(n,delta,psi)%*%(d1-d2)
  return(d3)
}
# mu.v.givenold(n,delta,psi,W.skew)

mu.v.given <- function(n,delta,psi,W.skew){
  H <- exp(-(normmatrix/psi))
  d1 <- delta*chol2inv(chol(H))%*%W.skew
  d2 <- sqrt(2/pi)*chol2inv(chol(H))%*%c(rep(1,n))
  d3 <- sigma.v.given(n,delta,psi)%*%(d1-d2)
  return(d3)
}
# mu.v.given(n,delta,psi,W.skew)

rv.givenold <- function(m,delta,psi,W.skew){
  d1 <- sigma.v.givenold(n,delta,psi)
  diag(d1) <- sqrt(diag(d1))
  # library(truncnorm)
  # d2 <- rtmvnorm(n=m,
  #                lower=rep(0, length=n),
  #                upper=rep(Inf, length=n),
  #                mean=c(mu.v.givenold(n,delta,psi,W.skew)),
  #                sigma=d1)
  d2 <- rtmvnorm(n=m,
                 mu=c(mu.v.givenold(n,delta,psi,W.skew)),
                 sigma=d1,
                 lb=rep(0, length=n),
                 ub=rep(Inf, length=n))
  return(d2)
}
# rv.givenold(m=1,delta,psi,W.skew)

rv.given <- function(m,delta,psi,W.skew){
  d1 <- sigma.v.given(n,delta,psi)
  diag(d1) <- sqrt(diag(d1))
  # library(truncnorm)
  # d2 <- rtruncnorm(n=m, a=0, b=Inf,
  #                mean=mu.v.given(n,delta,psi,W.skew),
  #                sd=d1)
  d2 <- rtmvnorm(n=m,
                 mu=c(mu.v.given(n,delta,psi,W.skew)),
                 sigma=d1,
                 lb=rep(0, length=n),
                 ub=rep(Inf, length=n))
  return(d2)
}
# rv.given(m=1,delta,psi,W.skew)

E2ti<-function(psi,beta0,beta1,gamma,teta,i){ # scaler
  d<-miuw.epsilonz(psi,beta0,beta1,gamma,teta)
  return(d[i])
}
# E2ti(psi,beta0,beta1,gamma,teta,2)

E3ti<-function(beta0,beta1,gamma,teta,tau2,i){  # scaler
  d<-miuepsiloni.wz(beta0,beta1,gamma,teta,tau2,i)
  return(d)
}
# E3ti(beta0,beta1,gamma,teta,tau2,3)

E4t<-function(beta0,beta1,gamma,teta,tau2){  # scaler
  d1<-miuepsilon.wz(beta0,beta1,gamma,teta,tau2)
  d2<-sigmaepsilon.wz(teta,tau2)
  d<-sum(d1^2+d2)
  return(d)
}
# E4t(beta0,beta1,gamma,teta,tau2)

E5t<-function(psi,beta0,beta1,gamma,teta){ # matrix
  d1<-sigmaw.epsilonz(psi,gamma,teta)
  d2<-miuw.epsilonz(psi,beta0,beta1,gamma,teta)
  d<-d1+d2%*%t(d2)
  return(d)
}
# E5t(psi,beta0,beta1,gamma,teta) 

E6t<-function(psi,beta0,beta1,gamma,teta){ # vector
  d<-miuw.epsilonz(psi,beta0,beta1,gamma,teta)
  return(d)
}
# E6t(psi,beta0,beta1,gamma,teta) 

E7t<-function(psi,beta0,beta1,gamma,teta,delta){ # scaler
  e1<-rw.epsilonz(L,psi,beta0,beta1,gamma,teta) 
  d1<-e1+sqrt(2/pi)*delta
  H <- exp(-(normmatrix/psi))
  e2<-H+diag(delta^2,n); e3<-chol2inv(chol(e2))
  d2<-delta*e3
  d3<-d2%*%t(d1)
  DELTA<-diag(1,n)-delta*d2
  result<-rep(0,L)
  for(i in 1:L){
    result[i]<-log(mvtnorm::pmvnorm(lower=rep(-Inf,n),upper=d1[i,],
                                    mean=as.vector(rep(0,n)),sigma=DELTA)[1] )
  }  
  Result<-mean(result)
  return(Result)
}
# E7t(psi,beta0,beta1,gamma,teta,delta) 

E8t<-function(psi,beta0,beta1,gamma,teta,i){ # scaler
  d1 <- miuw.epsilonz(psi,beta0,beta1,gamma,teta)[i,i]
  d2 <- sigmaw.epsilonz(psi,gamma,teta)[i,i]
  d <- ((d1^2)+d2)
  return(mean(d))}
# E8t(psi,beta0,beta1,gamma,teta) 

E9t<-function(beta0,beta1,gamma,teta,tau2){  # scaler
  d1<-miuepsiloni.wz(beta0,beta1,gamma,teta,tau2)
  d2<-sigmaepsilon.wz(teta,tau2)
  d<-d1^2+diag(d2)
  return(d)
}
# E9t(beta0,beta1,gamma,teta,tau2)

tau2.t1 <- function(E4.t0){
  d1 <- (1/n)*E4.t0
  return(d1)
}
# tau2.t1(E4.t0=4)

beta.t1 <- function(teta.t0, gamma.t0, E2.t0, E3.t0){
  d1 <- Landa(teta.t0)
  X.disg <- cbind(rep(1,n),X)
  first.term <- matrix(c(rep(0,4)),nrow=2)
  for (i in 1:n){
    first.term <- first.term + (d1[i] * (t(t(X.disg[i,])) %*% X.disg[i,]))
    # print(i)
  }
  first.term <- -chol2inv(chol(first.term))
  
  d2 <- gamma.t0*E2.t0 + E3.t0
  d3 <- matrix(c(rep(0,2)),nrow=2)
  for (i in 1:n){
    d3 <-  d3 + (d1[i]*d2[i])*t(t(X.disg[i,]))
  }
  d3 <- d3 + (n/4) + sum(Z.s/2)
  return(first.term %*% d3)
}
# beta.t1(teta,gamma, E2.t0=c(1:n), E3.t0=c(1:n))

gamma.t1 <- function(teta.t0, beta0.t0, beta1.t0, E2.t0, E3.t0, E8.t0){
  d1 <- Landa(teta.t0)
  X.disg <- cbind(rep(1,n),X)
  d2 <- c()
  for (i in 1:n){
    d2[i] <- (Z.s[i] -0.5 - 
                2*(d1[i])* (X.disg[i,]%*%t(t(c(beta0.t0,beta1.t0))) + E3.t0[i]))*E2.t0[i]
  }
  numer <- sum(d2)
  denomin <- 2*sum(d1*E8.t0)
  return(numer/denomin)
}
# gamma.t1(teta, beta0, beta1, E2.t0=c(1:n), E3.t0=c(1:n), E8.t0=c(1:n))

teta.t1 <- function(beta0.t1, beta1.t1, gamma.t1,
                    E8.t1, E9.t1, E2.t1, E3.t1){
  d1 <- c()
  d2 <- c()
  d3 <- c()
  X.disg <- cbind(rep(1,n),X)
  for (i in 1:n){
    d1[i] <- (X.disg[i,]%*%t(t(c(beta0.t1,beta1.t1))))^2 +
      gamma.t1^2*E8.t1[i] + E9.t1[i]
    d2[i] <- 2*gamma.t1*(X.disg[i,]%*%t(t(c(beta0.t1,beta1.t1))))*E2.t1[i]
    d3[i] <- 2*(X.disg[i,]%*%t(t(c(beta0.t1,beta1.t1))))*E3.t1[i] + 
      2*gamma.t1*E2.t1[i]*E3.t1[i]
  }
  return(d1+d2+d3)
}
# teta.t1(beta0,beta1,gamma,E8.t1=c(1:n),E9.t1=c(1:n),E2.t1=c(1:n),E3.t1=c(1:n))

delta.t1.f <- function(delta, psi.t0, E5.t0, E6.t0, E7.t0){
  H.t0 <- exp(-(normmatrix/psi.t0))
  d0 <- H.t0 + diag(delta^2, nrow=n)
  #
  d1 <- -log(det(d0))
  d2 <- -(1/2)*sum(diag( chol2inv(chol(d0))%*%E5.t0 ))
  d3 <- -sqrt(2/pi)*c(rep(1,n)) %*% chol2inv(chol(d0)) %*% t(t(E6.t0))
  d4 <- -(1/pi)*(delta^2)*c(rep(1,n)) %*% 
    chol2inv(chol(d0)) %*% t(t(c(rep(1,n))))
  d5 <- c(d1+d2+d3+d4+E7.t0)
  return(d5)
}
# delta.t1.f(delta, psi, E5.t0=exp(-normmatrix/psi), E6.t0=c(rep(1,n)), E7.t0=2)

Opti1 <- function(psi.t0, E5.t0, E6.t0, E7.t0, LB, UB){
  methods1 <- c("Nelder-Mead","BFGS","CG","L-BFGS-B","nlm",
                "nlminb","spg","ucminf","newuoa","bobyqa","nmkb",
                "hjkb","Rcgmin","Rvmmin")
  Naiv <- optimx(1,delta.t1.f,
                 psi.t0=psi.t0, E5.t0=E5.t0, E6.t0=E6.t0, E7.t0=E7.t0,
                 lower=LB, upper=UB,
                 method=methods1, control=list(maximize=TRUE))$value
  result <- c()
  acc <- Naiv[LB < Naiv & Naiv < UB]
  if(is.null(acc)==FALSE){
    if(length(acc)!=0){
      result <- mean(acc)
    }
  }
  if(is.null(result)==TRUE){
    result <- rnorm(1,delta,0.1)
  }
  return(result)
}

# Opti1(psi, E5.t0=exp(-normmatrix/psi),
#       E6.t0=c(rep(1,n)), E7.t0=2,
#       LB=0, UB=20)
######################## Run
# Naiv <- optimx(1,delta.t1,
#        psi.t0=psi, E5.t0=exp(-normmatrix/psi),
#        E6.t0=c(rep(1,n)), E7.t0=2,
#        lower=0, upper=20,
#        method=methods1, control=list(maximize=TRUE))$value
# acc <- Naiv[0 < Naiv & Naiv < 20]

psi.t1.f <- function(psi, delta.t0, E5.t0, E6.t0, E7.t0){
  H.t0 <- exp(-(normmatrix/psi))
  d0 <- H.t0 + diag(delta.t0^2, nrow=n)
  #
  d1 <- -log(det(d0))
  d2 <- -(1/2)*sum(diag( chol2inv(chol(d0))%*%E5.t0 ))
  d3 <- -sqrt(2/pi)*c(rep(1,n)) %*% chol2inv(chol(d0)) %*% t(t(E6.t0))
  d4 <- -(1/pi)*(delta.t0^2)*c(rep(1,n)) %*% 
    chol2inv(chol(d0)) %*% t(t(c(rep(1,n))))
  d5 <- c(d1+d2+d3+d4+E7.t0)
  return(d5)
}
# psi.t1.f(psi, delta, E5.t0=exp(-normmatrix/psi), E6.t0=c(rep(1,n)), E7.t0=2)

Opti2 <- function(delta.t0, E5.t0, E6.t0, E7.t0, LB, UB){
  methods1 <- c("Nelder-Mead","BFGS","CG","L-BFGS-B","nlm",
                "nlminb","spg","ucminf","newuoa","bobyqa","nmkb",
                "hjkb","Rcgmin","Rvmmin")
  Naiv <- optimx(1,psi.t1.f,
                 delta.t0=delta.t0, E5.t0=E5.t0, E6.t0=E6.t0, E7.t0=E7.t0,
                 lower=LB, upper=UB,
                 method=methods1, control=list(maximize=TRUE))$value
  result <- c()
  acc <- Naiv[LB < Naiv & Naiv < UB]
  if(is.null(acc)==FALSE){
    if(length(acc)!=0){
      result <- mean(acc)
    }
  }
  if(is.null(result)==TRUE){
    result <- rnorm(1,psi,0.1)
  }
  return(result)
}

# Opti2(delta, E5.t0=exp(-normmatrix/psi),
#       E6.t0=c(rep(1,n)), E7.t0=2,
#       LB=0, UB=20)













