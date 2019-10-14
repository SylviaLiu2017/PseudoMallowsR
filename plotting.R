#####does varying number of users matter?
n<-25

load(paste("./results/ResultTableN50n",n,".RData",sep=""))
result1<-resultTable
load(paste("./results/ResultTableN100n",n,".RData",sep=""))
result2<-resultTable
load(paste("./results/ResultTableN200n",n,".RData",sep=""))
result3<-resultTable

whichMin<-function(vec){
  return(which(vec == min(vec))[1])
}

minInd1<-apply(result1[,3:13],1,whichMin)
minInd2<-apply(result2[,3:13],1,whichMin)
minInd3<-apply(result3[,3:13],1,whichMin)

par(mfrow = c(1,3))
plot(result1[,2],minInd1,main='50 users',xlab = 'dataset sd', ylab = 'optimal sd')
plot(result2[,2],minInd2,col='red',main = '100 users',xlab = 'dataset sd', ylab = 'optimal sd')
plot(result3[,2],minInd3,col='blue',main='200 users',xlab = 'dataset sd', ylab = 'optimal sd')

par(mfrow = c(1,3))
plot(result1[,1],minInd1,xlim=rev(range(result1[,1])),main='50 users',xlab = 'dataset alpha', ylab = 'optimal sd')
plot(result2[,1],minInd2,xlim=rev(range(result2[,1])),col='red',main = '100 users',xlab = 'dataset alpha', ylab = 'optimal sd')
plot(result3[,1],minInd3,xlim=rev(range(result3[,1])),col='blue',main='200 users',xlab = 'dataset alpha', ylab = 'optimal sd')



##############same n , N, alpha, just different dataset############
par(mfrow = c(1,1))
line1<-result1[1,3:13]
line2<-result1[2,3:13]
line3<-result1[3,3:13]

plot(line1,type='b',ylim=c(min(line1,line2,line3),max(line1,line2,line3)))
lines(line2,col='blue',type='b')
lines(line3,col='red',type='b')

line4<-result1[16,3:13]
line5<-result1[18,3:13]
line6<-result1[19,3:13]
plot(line4,type='b',ylim=c(min(line4,line5,line6),max(line4,line5,line6)))
lines(line5,col='blue',type='b')
lines(line6,col='red',type='b')

###############different n, same N, same sd###################
load("./results/ResultTableN100n10.RData")
result1<-resultTable
load("./results/ResultTableN100n15.RData")
result2<-resultTable
load("./results/ResultTableN100n20.RData")
result3<-resultTable
load("./results/ResultTableN100n25.RData")
result4<-resultTable

minInd1<-apply(result1[,3:13],1,whichMin)
minInd2<-apply(result2[,3:13],1,whichMin)
minInd3<-apply(result3[,3:13],1,whichMin)
minInd4<-apply(result4[,3:13],1,whichMin)

par(mfrow = c(2,2))
plot(result1[,2],minInd1,main='10 items',xlab = 'dataset sd', ylab = 'optimal sd')
abline(v=2.0,lty =4)
plot(result2[,2],minInd2,col='red',main = '15 items',xlab = 'dataset sd', ylab = 'optimal sd')
abline(v=3.0,lty =4,col='red')
plot(result3[,2],minInd3,col='blue',main='20 items',xlab = 'dataset sd', ylab = 'optimal sd')
abline(v=4,lty =4,col='blue')
plot(result4[,2],minInd4,col='green',main='25 items',xlab = 'dataset sd', ylab = 'optimal sd')
abline(v=5.0,lty =4,col='green')

##################once again#########################
load("./results/ResultTableN100n10.RData")
result1<-resultTable
load("./results/ResultTableN100n15.RData")
result2<-resultTable
load("./results/ResultTableN100n20.RData")
result3<-resultTable
load("./results/ResultTableN100n50.RData")
result4<-resultTable

minInd1<-apply(result1[,3:13],1,whichMin)
minInd2<-apply(result2[,3:12],1,whichMin)
minInd3<-apply(result3[,3:12],1,whichMin)
minInd4<-apply(result4[,3:12],1,whichMin)

par(mar=c(2,2,2,2))
plot(result1[,2],minInd1,main='10 items',xlab = 'dataset sd', ylab = 'optimal sd')
plot(result2[,2],minInd2,col='red',main = '15 items',xlab = 'dataset sd', ylab = 'optimal sd')
plot(result3[,2],minInd3,col='blue',main='20 items',xlab = 'dataset sd', ylab = 'optimal sd')
plot(result4[,2],minInd4,col='green',main='50 items',xlab = 'dataset sd', ylab = 'optimal sd')

dev.off()
par(mfrow=c(2,2))
plot(result2[,2],minInd2,main='10 items',xlab = 'dataset sd', ylab = 'optimal sd')
plot(result1[,2],minInd1,main='15 items',xlab = 'dataset sd', ylab = 'optimal sd')

