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
alpha0<-c(1.6)
alpha0_orderings<-c(10)
N<-100
rho0<-1:n
sourceCpp('MCMC_old.cpp')
sds<-c(0.1,0.3,0.5,1,3,5,7,9,15,20)
resultTable<-matrix(data=NA, nrow = 1,ncol = length(sds)+2)
colnames(resultTable)<-c("alpha data", "sd data", sds)
for(alpha0 in c(0.5,0.7,1,1.2,1.4,1.5,1.7,2,3,4,5)){
  KLs<-vector()
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
n_samples<-4000
tmpRank<-rank(apply(data,2,sum),ties.method='first')
rhoMat<-matrix(data=NA, nrow = n_samples,ncol = n)

for(sdNorm in sds){
  print(paste("SD = ",sdNorm))
  orderings <- vector()
  for(i in 1:n_samples){
    orderings <- rbind(orderings,rank(rnorm(n, mean =generateVOrderings(tmpRank),sd = sdNorm)))
  }
  
  for(sample_i in 1:n_samples){
    if(sample_i %% 100 ==0){
      print(sample_i)
    }
    support<-1:n
    rho<-rep(NA,n)
    #In<-generateVOrderings(tmpRank)
    In<-orderings[sample_i,]
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
  heatMat_ML<-heatMap(rhoMat_ML,rho0)
  heatMat_pseudo<-heatMap(rhoMat,rho0)
  KL<-0
  for(k in 1:n){
    KL<-KL+KL_margin(heatMat_ML[1,],heatMat_pseudo[1,],margin=0.000001)
  }
  
  KLs<-c(KLs,KL)

}


par(mai=c(1,1,0.65,1))
letters <- c(1:n)
image(heatMat_ML,col=tim.colors(64*10),zlim=c(0,1),axes=F,cex.lab=2,main = "Mallows")
par(mai=c(1,1,0.65,1))
image.plot(heatMat_ML, zlim=c(0,1),legend.only=T,horizontal = F)

par(mai=c(1,1,0.65,1))
letters <- c(1:n)
image(heatMat_pseudo,col=tim.colors(64*10),zlim=c(0,1),axes=F,cex.lab=2,main = "Pseudo")
par(mai=c(1,1,0.65,1))
image.plot(heatMat_ML, zlim=c(0,1),legend.only=T,horizontal = F)



resultTable<-rbind(resultTable,c(alpha0,mean(apply(data,2,sd)),KLs))
}

resultTable<-na.exclude(resultTable)
resultTable<-resultTable[order(resultTable[,2],decreasing = TRUE),]
save(resultTable,file=paste("./tmpResultTableN",N,"n",n,".RData",sep=""))

rowsToPlot<-1:30

plot(sds,resultTable[min(rowsToPlot),3:dim(resultTable)[2]],type='b',ylim=c(min(resultTable[rowsToPlot,3:dim(resultTable)[2]]),max(resultTable[rowsToPlot,3:dim(resultTable)[2]])))
for(i in min(rowsToPlot):max(rowsToPlot)){
  lines(sds,resultTable[i,3:dim(resultTable)[2]],type='b',ylim=c(min(resultTable),max(resultTable)),col=i)
}

whichMin<-function(vec){
  return(which(vec == min(vec)))
}

minIndex<-apply(resultTable[,3:dim(resultTable)[2]],1, whichMin)

plot(resultTable[,2],sds[minIndex])
