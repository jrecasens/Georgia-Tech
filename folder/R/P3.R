# Problem 3

# Clear memory, include libraries
rm(list=ls())
require("ggplot2")

# Read Data
data = read.csv(file="OPCA.csv", head=FALSE, sep =",")
totalrows = dim(data)[1]
d = dim(data)[2]
n = totalrows/d

# Create list of A[i]
A = lapply(1:n, function (x) as.matrix(data[((x-1)*d+1):( x*d),]))

# Get True Eigenvector
true_eigenvector = as.matrix(read.csv(file="True_eigvector.csv", head=FALSE, sep =","))

#Define Step Function
step = function (i, q){
  return (1/(100 + i*q))
}

#Visualization stuff
nlist = c(0:n)
part = c("a)","b)")
part[1]

# Oja's algorithm
Oja = function (q){
  # Initialization: w0
  w = matrix(1/sqrt(d), nrow=d, ncol = 1)
  #Simmilarity for Plot
  similarity = rep(0,n+1)
  similarity[1] = 1 - (t(w)%*%true_eigenvector)^2
  #Iterate through n
  for(i in 1:n){
    w = w + step(i,q)* A[[i]]%*%w
    w = w/norm(w, type = c("2"))
    similarity[i+1] = 1 - (t(w)%*%true_eigenvector)^2
  }
  #----CREATE PLOT-----
  p = part[q+1]
  Plot <- ggplot(data = as.data.frame(cbind(nlist,similarity)), 
                 aes(x = nlist, y = similarity)) + 
    geom_point(size = 0.7, alpha = 1/3) + 
    xlab("Iteration Number") + 
    ylab("dist(wi; v)") + 
    ggtitle(bquote(list("Measure of similarity v/s iteration. Problem", .(p))))
  print(Plot)
  ggsave(Plot,filename=paste("P3_Plot",p,".png",sep=""),width = 5, height = 2.5)
  print(similarity[length(similarity)])
}

Oja(0)
Oja(1)

