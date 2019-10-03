rm(list=ls())

require(Rcpp)
require(RcppArmadillo)
require(BayesMallows)
require(fields)

#########################pseudo denominator calculator############################
pseudoDenom<-function(alpha,data,rho,ordering){
  n<-dim(data)[2]
  logdenom_marginal<-rep(NA, n)
  log_num_joint<- (-alpha/n)* sum(abs(t(data)-rho))
  log_denom<-0
  for(i in 1:n){
    item<-which(ordering==i)
    supports<-rho[ordering>=i]
    factor_i<-0
    for(support in supports){
      factor_i<-factor_i + exp(-(alpha/n)*sum(abs(data[,item] - support)))
    }
    logdenom_marginal[i] <-log(factor_i)
  }
  
  return(logdenom_marginal)
}

################generate all permutations###################
source("./allFunctions.R")
n<-20
load("./Cdfootrule.RData")
#################generate some data###################
alpha0<-c(2)
alpha0_orderings<-c(10)
N<-200
rho0<-1:n
sourceCpp('MCMC_old.cpp')

data<-sample_mallows(rho0 = rho0,alpha0 = alpha0, n_samples = N)
###############calculating Mallows posterior probability, ground truth####################
Rtrue<-rho0
nmc <- 500000 # Number of iterations of MCMC algorithm
thin <- 100    # Number of thinning iterations for saving purposes
alphaJump <- 10
distance<-"footrule"
L<-2
sdAlpha <- 0.5
lambda <- 1/10
fitvec<-0
alphaInit<-1
pars <- list("nmc"=nmc, "thin"=thin, "alphaJump"=alphaJump, "lambda"=lambda, "L"=L, "sdAlpha"=sdAlpha,
             "dist"=distance, "fit"=fitvec, "Cd"=seq2, aug=FALSE, Rmiss=sample(n), rho0=sample(n), alpha0 = alphaInit)
result<- MCMCfunction(data, pars)

rhoMat_ML<-result$rho[1001:5000,]
heatMat_ML<-heatMap(rhoMat_ML,rho0)

################using Pseudo Likelihood#####################
n_samples<-500
perturbMethod<-'V'
start <- proc.time()
tmpRank<-rank(apply(data,2,sum),ties.method='first')
rhoMat<-matrix(data=NA, nrow = n_samples,ncol = n)

for(sample_i in 1:n_samples){
  print(sample_i)
  support<-1:n
  rho<-rep(NA,n)
  In<-generateVOrderings(tmpRank)
  for(i in 1:n){
    i_curr<-which(In==i)
    dist<-sapply(support, oneDimfootrule, Rs=data[,i_curr])
    log_num<-(-alpha0/(n)*(dist))
    log_denom<- log(sum(exp(log_num)))
    probs<-exp((log_num-log_denom))
    rand<-runif(1)
    indOfCdf<-length(support)+1-sum(rand<=cdfForSamples(probs))
    rho[i_curr]<-support[indOfCdf]
    support<-setdiff(support,rho[i_curr])
  }
  rhoMat[sample_i,]<-rho
}




