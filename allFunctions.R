
##########run this first##############
vecEquate<-function(vec1,vec2){
  n<-length(vec1)
  return(sum(vec1==vec2) == n)
}

entropy<-function(dist){
  nonZeros<-dist[dist>0]
  return(sum(nonZeros * log(nonZeros)))
}
generateAllPerm<-function(n){
  numSamples<-50*factorial(n)
  samples<-matrix(data=NA, nrow=numSamples,ncol=n)
  for(i in 1:numSamples){
    print(i)
    set.seed(i)
    samples[i,]<-sample(n)
  }

  allNum<-apply(samples,1,convertToNumber)
  length(unique(allNum))==factorial(n)

  allCombinations<-unique(allNum)
  allPerm<-t(sapply(allCombinations,backTovec))

  return(allPerm)
}

KL_margin<-function(distP,distQ, margin){
  distP_mod<-distP
  distQ_mod<-distQ
  zeroP<-which(distP == 0)
  nonzeroP<-which(distP > 0)
  zeroQ<-which(distQ == 0)
  nonzeroQ<-which(distQ > 0)

  if(length(zeroP)>0){
    distP_mod[zeroP]<-margin
    distP_mod[nonzeroP]<-distP[nonzeroP]-margin*length(zeroP)/length(nonzeroP)
  }
  if(length(zeroQ)>0){
    distQ_mod[zeroQ]<-margin
    distQ_mod[nonzeroQ]<-distQ[nonzeroQ]-margin*length(zeroQ)/length(nonzeroQ)
  }
  if((sum(distP_mod<0)+sum(distQ_mod<0))>0){
    print("margin too large")
  }

  return(sum(log(distP_mod / distQ_mod) *distP_mod))
}

oneDimfootrule<-function(rho_i,Rs){
  return(sum(abs(Rs-rho_i)))
}

convertToNumber<-function(vec){
  n<-length(vec)
  num<-0
  for(j in 1:n){
    num<-num+vec[n-j+1]*10^(j-1)
  }
  return(num)
}

backTovec<-function(num){
  n<-ceiling(log10(num))
  vec<-rep(NA,n)
  for(i in 1:n){
    vec[n-i+1]<-num %% 10
    num<-num %/%10
  }
  return(vec)
}


MallowsZ_footRule_log<-function(intSeq,alpha){
  n<-length(intSeq)
  if((n %% 2) == 0){
    m<-n/2
    Di=2*seq(0, m^2, 1)
  }else{
    m <-(n-1)/2
    Di=2*seq(0, m^2+m, 1)
  }
  Cdn <-intSeq
  res = log(sum((Cdn)*exp(-alpha*Di/n)))
  return(res)
}


pseudoJoint<-function(data, alpha, rho, ordering){
  n<-length(ordering)
  log_num_joint<- (-alpha/n)* sum(abs(t(data)-rho))
  log_denom<-0
  for(i in 1:n){
    item<-which(ordering==i)
    supports<-rho[ordering>=i]
    factor_i<-0
    for(support in supports){
      factor_i<-factor_i + exp(-(alpha/n)*sum(abs(data[,item] - support)))
    }
    log_denom<-log_denom+ log(factor_i)
  }
  prob<- exp(log_num_joint - log_denom)
  return(prob)
}

KL_divergence<-function(distP, distQ){
  return(sum(log(distP / distQ) *distP))
}


freqMat<-function(allPerm,probs){
  n<-dim(allPerm)[2]
  freqMap<-matrix(data=NA, nrow=n,ncol=n)
  for(i in 1:n){
    for(j in 1:n){
      freqMap[i,j]<-sum(probs[which(allPerm[,i]==j)])
    }
  }
  return(freqMap)
}


distToData<-function(data, vec){
  return(  sum(abs(t(data)-vec)))
}
cdfForSamples<-function(probs){
  cdf_result<-vector()
  numOfInterval<-length(probs)
  for(int in 1:numOfInterval){
    cdf_result[int]<-sum(probs[1:int])
  }
  return(cdf_result)
}


PL_SylviaDefined<-function(score){
  n<-length(score)
  support<-1:n
  ordering<-rep(NA,n)
  for(i in 1:n){
    probs<-score[support]/sum(score[support])
    cdf<-cdfForSamples(probs)
    rand<-runif(1)
    ordering[i]<-support[length(support)+1-sum(rand<=cdf)]
    support<-setdiff(support,ordering[i])
  }
  return(ordering)
}

####################visualization########################
heatMap<-function(rhoMat,Rtrue){
  n<-dim(rhoMat)[2]
  freqMat<-matrix(data=NA, nrow=n,ncol=n)
  for(j in 1:n){
    ind <- which(Rtrue == j)
    freqMat[j, ] <- sapply(1:n, function(x) sum(rhoMat[,ind] == x))
  }
  freqMat <- freqMat/dim(rhoMat)[1]
  return(freqMat)
}

####################
generateVOrderings<-function(centre){
  n<-length(centre)
  ordering<-rep(NA,n)
  distToCentre<-centre - (max(centre)+1)/2
  if(n%%2==1){
    ordering[which(distToCentre==0)]<-1
  }
  for(i in distToCentre[distToCentre>0]){
    if(runif(1)>=0.5){
      ordering[which(distToCentre==i)]<-abs(i)*2+1
      ordering[which(distToCentre==(-1)*i)]<-abs(i)*2
    }else{
      ordering[which(distToCentre==i)]<-abs(i)*2
      ordering[which(distToCentre==(-1)*i)]<-abs(i)*2+1
    }
  }
  return(ordering)
}


isV<-function(permutation,centre){
  n<-length(permutation)
  trueFalse<-rep(NA, n)
  if(n %%2 ==1){
    middle<-(1+n)/2
    trueFalse[middle]<-(permutation[which(centre == middle)]==1)
    for(i in 1:((n-1)/2)){
      trueFalse[i]<-as.numeric(permutation[which(centre==i)]==(n-2*i+1))+as.numeric(permutation[which(centre==i)]==(n-2*i+2))
      trueFalse[(n-i+1)]<-as.numeric(permutation[which(centre==(n-i+1))]==(n-2*i+1))+as.numeric(permutation[which(centre==(n-i+1))]==(n-2*i+2))
    }
  }else{
    for(i in 1:((n)/2)){
      trueFalse[i]<-as.numeric(permutation[which(centre==i)]==(n-2*i+1))+as.numeric(permutation[which(centre==i)]==(n-2*i+2))
      trueFalse[(n-i+1)]<-as.numeric(permutation[which(centre==(n-i+1))]==(n-2*i+1))+as.numeric(permutation[which(centre==(n-i+1))]==(n-2*i+2))
    }
  }
  return(sum(trueFalse)==n)
}
