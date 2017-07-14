# Problem 2

# Clear memory, include libraries
rm(list=ls())
require("ggplot2")

# Read Data
data = read.csv(file="MLR.csv", head=FALSE, sep =",")

#Prepare Data:
# number of dimensions and observations
n = dim(data)[1]
p = dim(data)[2]-1

# Define Design Matrix X and Response
X = as.matrix(data[,1:p])
Y = as.matrix(data[,p+1])
data = as.matrix(data)

# implement feature scaling
#X.scaled <- scale(X, center = TRUE, scale = TRUE)
#Y.scaled <- scale(Y, center = TRUE, scale = TRUE)
#data.scaled <- scale(data, center = TRUE, scale = TRUE)
X.scaled <- X
Y.scaled <- Y

##### Stochastic Gradient descent  ##### 

#----ALGORITHM INPUTS----- 
#differentiable function g()
g = function (X.scaled, Y.scaled, beta){
  return (sum((Y.scaled-X.scaled %*% beta)^2)/(2*n))
}
#Lipschitz
L = ((norm(X.scaled,type="2"))^2)/ n
# gradient function
avg_gradient = function (X.scaled, Y.scaled, beta , b){
  minibatch = sample (1:n, b, replace = TRUE)
  minigrad = matrix(0, nrow = p, ncol = 1)
  for(k in minibatch){
    minigrad = minigrad - (X.scaled[k,]%*%beta - Y.scaled[k])*X.scaled[k,]
  }
  return (- 1/b * minigrad )
}

#----BEGIN ALGORITHM per b-----
for(bb in c(10,25,100,1)){
      #fixed step length
      alpha = bb/(n*L)
      
      #initial point beta_0
      beta = matrix(0, nrow = p, ncol = 1)
      #----BEGIN ALGORITHM-----
      threshold = 0.0000001
      max_i = 1000
      FunctionValues = rep(0,max_i)
      FunctionValues[1] = g(X.scaled,Y.scaled,beta)
      plotseq = seq(1:max_i)
      best_beta = 0
      i = 2
      while (i <= max_i){
        grad = avg_gradient(X.scaled,Y.scaled,beta, bb)
        beta = beta - alpha * grad
        FunctionValues[i] = g(X.scaled,Y.scaled,beta)
        #Convergence Check
        if((norm(grad, type ="2") < threshold) 
           && (FunctionValues[i] - FunctionValues[i-1] < threshold) ){
          best_beta = beta
          print("break")
          break 
        }
        i=i+1
        if(i <= max_i) {best_beta = beta}
      }

      #----CREATE PLOT-----
      Plot <- ggplot(data = as.data.frame(cbind(plotseq,FunctionValues)), 
                     aes(x = plotseq, y = FunctionValues)) + 
        geom_point(size = 0.7, alpha = 1/3) + 
        xlab("Iteration Number") + 
        ylab("g(beta)") + 
        ggtitle(bquote(list("g(beta) versus number of iterations. b", .(bb))))
      print(Plot)
      ggsave(Plot,filename=paste("P2_Plot",bb,".png",sep=""),width = 5, height = 2.5)
      
      #----CHECK MSE-----
      true_beta = as.matrix(read.csv(file ="True_Beta.csv", head=FALSE ,sep =","))
      MSE = sum((best_beta - true_beta) ^ 2) / n
      print(MSE)

}