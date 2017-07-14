# ISyE 6740 - Homework 1
# Javier Recasens

# Clear memory, include libraries
rm(list=ls())
require("ggplot2")
require("dplyr")
require("mvtnorm")
require("RColorBrewer")

#/////////////////Load Data///////////////////
# Read in data
myData=read.csv("semeion.csv",header=FALSE)
# Covert Data Frame to Matrix
myX=data.matrix(myData[,1:256])
# Get true Labels
myLabel=apply(myData[,257:266],1,function(xx){
  return(which(xx=="1")-1)
})

# Set initial parameters
k=10
d=dim(myX)[2]
n=dim(myX)[1]
iterations=25
PCs=c(0,2,4,6)

#/////////////////Create Functions///////////////////

# FUNCTION: E Step
EStep <- function(myX, mu, pi, VCV){
  px <- matrix(0,nrow=n, ncol=k)
  for(c in 1:k) {px[,c] <- pi[c]*dmvnorm(myX, mu[,c], VCV[c,,], log = FALSE)}
  return(px)
}

# FUNCTION: M Step
MStep <- function(myX, Gamma, pc){
  #Calculate mean vector and pi for all clusters
  SumM <- t(myX)%*%Gamma
  GammaSums <- colSums(Gamma)
  mu = SumM%*% diag(1/GammaSums)
  pi = GammaSums / n
  #Calculate the Variance Covariance Matrix
  VCV <- array(dim=c(k,d,d))
  for(c in 1:k) {
    CovTemp <- matrix(0,d,d)
    for(i in 1:n) {
      xbar <- myX[i,] - mu[,c]
      CovTemp <- CovTemp + Gamma[i,c]*((xbar)%*%t(xbar))
    }
    CovTemp <- CovTemp / GammaSums[c]
    # Do the Spectral Decomposition on CovTemp
    Decomposition <- eigen(CovTemp, symmetric=TRUE)
    PCvectors <- Decomposition$vectors[,1:pc]
    sigma2 <- sum(Decomposition$values[pc+1:d], na.rm = TRUE)/(d-pc)
    if (pc != 0) {
      Wq <- PCvectors%*%diag(sqrt(Decomposition$values[1:pc] - sigma2))
      VCV[c,,] <-  Wq%*%t(Wq) + sigma2*diag(d)
    } else {
      VCV[c,,] <-  sigma2*diag(d)
    }
  }
  return(list(mu,pi,VCV))
}


#/////////////////EM Algorithm///////////////////

set.seed(99)
loglike <- matrix(0,length(PCs),iterations)
AIC <- rep(0,length(PCs))
Xline <- seq(1:iterations)
Yline <- c()

# EM algorithm
for(q in PCs) {
    PCi <- which(PCs == q)
    PCplot <- toString(q)
    # k-means for initialization
    myCluster <- kmeans(myX, k, iter.max = 20, nstart = 10)
    # Initialization of Gamma and parameters
    Gamma <- matrix(0, nrow=n, ncol=k)
    
    for(i in 1:n) {Gamma[i, myCluster$cluster[i]] <- 1}
    phi <- MStep(myX, Gamma, q)
    
    #Get initialization of the 3 parameters: Mean of clusters, probability per cluster
    #and the variance covariance matrix.
    mu <- phi[[1]]
    pi <- phi[[2]]
    VCV <- phi[[3]]
    
    oldloglike <- 0
    
    for(it in 1:iterations) {
        #E Step
        px <- EStep(myX, mu, pi, VCV)
        Gamma <- px / rowSums(px)
        
        #M Step
        phi <- MStep(myX, Gamma, q)
        mu <- phi[[1]]
        pi <- phi[[2]]
        VCV <- phi[[3]]
        
        #/////////////////Log-Likelihood Convergence///////////////////
        loglike[PCi,it] <- sum(log(rowSums(px)))
        if(abs(loglike[PCi,it] - oldloglike) < 1){
            loglike[PCi,it:iterations] = oldloglike
            break
        }
        oldloglike <- loglike[PCi,it]
    }
    #AIC
    AIC[PCi] <- -2*loglike[PCi, iterations] + 2*(d*q + 1 - q*(q-1)/2)
    #log-likelihood Plot per Principal Component
    Yline <-loglike[PCi,]
    loglikePlot <- ggplot(
      data = as.data.frame(cbind(Xline,Yline)), aes(x = Xline, y = Yline)) 
    + geom_point() + xlab("Iteration Number") + ylab("Log-likelihood") 
    + ggtitle(bquote(list("Log Likelihood for q", .(PCplot))))
    print(loglikePlot)
    ggsave(loglikePlot,filename=paste("loglikePlot",PCplot,".png",sep=""), 
           width = 3, height = 3)
}

#/////////////////Choice of Number of Principle Components///////////////////
#Number of PC that minimize AIC
OptimalQ <- PCs[which(AIC == min(AIC))]
OptimalQ

#/////////////////Visualization of Clusters///////////////////

pal <- colorRampPalette(brewer.pal(9,"BuPu"))(100)

dev.new(width=7,height=3.5)
par(mai=c(0.05,0.05,0.05,0.05),mfrow=c(10,6))
for(c in 1:k){
  image(t(matrix(mu[,c], byrow=TRUE,16,16)[16:1,]),col=pal,axes=FALSE)
  box()
  for(j in 1:5){
    Normaldist <- rmvnorm(n=1, mean=mu[,c], sigma=VCV[c,,])
    image(t(matrix(Normaldist, byrow=TRUE,16,16)[16:1,]),col=pal,axes=FALSE)
    box()
  }
}

#/////////////////Accuracy Assessment///////////////////

#Get maximum probability from Gamma
myLabelEM <- apply(Gamma, 1, which.max)

# Group Data and select most common
Groups <- split(myLabel, myLabelEM)
PointsPerCluster <- unlist(lapply(Groups, length))
MostCommon <- lapply(Groups, function(g){return(sort(table(g), decreasing=TRUE)[1])})

# Accuracy
AccuracyMatrix <- matrix(0,10,4)

for(c in 1:k) {
  #The most common number
  AccuracyMatrix[c,1] <- as.integer(names(MostCommon[[c]]))
  #The most common number repetitions
  AccuracyMatrix[c,2] <- as.integer(MostCommon[[c]][[1]])
  #PointsPerCluster
  AccuracyMatrix[c,3] <- PointsPerCluster[c]
  #Accuracy per cluster
  AccuracyMatrix[c,4] <- 1 - (AccuracyMatrix[c,2] / AccuracyMatrix[c,3])
  #Print result per cluster
  print(c('Cluster ',c, 'for number:', AccuracyMatrix[c,1],
            'has',round(AccuracyMatrix[c,4], 2)*100, "% mis-cat. rate" ))
}

OverallAccuracy <- sum(AccuracyMatrix[,2]) / sum(AccuracyMatrix[,3])

print(c('The overall mis-categorization rate is:', round(1-OverallAccuracy,2)*100,'%'))

write.csv(AccuracyMatrix, "AccuracyMatrix.csv")







