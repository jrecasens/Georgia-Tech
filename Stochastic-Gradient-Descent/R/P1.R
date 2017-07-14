# Problem 1

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

##### Gradient descent (with fixed step length) ##### 

#----ALGORITHM INPUTS----- 
# differentiable function f()
f = function (X.scaled, Y.scaled, beta ){
  return (sum((Y.scaled-X.scaled %*% beta)^2)/(2*n))
}
# fixed step length
L = ((norm(X.scaled,type="2"))^2)/ n
alpha = 1/L
# gradient function
gradient = function (X.scaled, Y.scaled, beta ){
  return ( 1/n * t(X.scaled) %*% (X.scaled %*% beta - Y.scaled))
}
#initial point beta_0
beta = matrix(0, nrow = p, ncol = 1)
#----BEGIN ALGORITHM-----
threshold = 5.053251e-20
max_i = 100
FunctionValues = rep(0,max_i)
FunctionValues[1] = f(X.scaled,Y.scaled,beta)
plotseq = seq(1:max_i)
best_beta = 0
i = 2
while (i <= max_i){
  g = gradient(X.scaled,Y.scaled,beta)
  beta = beta - alpha * g
  FunctionValues[i] = f(X.scaled,Y.scaled,beta)
  #Convergence Check
  if((norm(g, type ="2") < threshold) 
     && (FunctionValues[i] - FunctionValues[i-1] < threshold) ){
      best_beta = beta
      print("break")
      break 
  }
  i=i+1
  if(i <= max_i) {best_beta = beta}
}
best_beta

#----CREATE PLOT-----
P1d_Plot <- ggplot(data = as.data.frame(cbind(plotseq,FunctionValues)), 
                   aes(x = plotseq, y = FunctionValues)) + 
  geom_point(size = 0.8) + 
  xlab("Iteration Number") + 
  ylab("f(beta)") + 
  ggtitle(bquote(list("f(beta) versus number of iterations")))
print(P1d_Plot)
ggsave(P1d_Plot,filename=paste("P1d_Plot.png",sep=""),width = 5, height = 2.5)

#----CHECK MSE-----
true_beta = as.matrix(read.csv(file ="True_Beta.csv", head=FALSE ,sep =","))
MSE = sum((best_beta - true_beta) ^ 2) / n
MSE