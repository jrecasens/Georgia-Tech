# Clean console and variables , time computation

rm(list=ls())

data = read.csv( file ="MLR.csv", head =FALSE ,sep =",")
set.seed (456)
# read dimensions and values
n = nrow ( data )
p = ncol ( data )-1
X = as.matrix ( data [ ,1:p])
Y = as.matrix ( data [,p +1])
#Set L size
L= norm(t(X )%*%X, type = "2")/ n
# define g(x) function :
g = function (X, Y, Beta ){
  return (sum ((Y-X %*% Beta )^2 ) / (2*n))
}
# Define gradient function
gradient = function (X, Y, Beta , b){
  rands = sample (1:n, b, replace = TRUE )
  grad_vector = matrix (0, nrow = p, ncol = 1)
  for(r in rands ){
    grad_vector = grad_vector + (Y[r]-X[r,]%*% Beta )*X[r ,]
  }
  return (- 1/b * grad_vector )
}

#run the whole stochastic gradient with mini - batch
stoch_gradient = function (b, pngname){
  #Set step size
  step =b/n/L
  # Initialize
  tolerance = 0.001
  max_iter = 8000
  Beta = matrix (0, nrow = p, ncol = 1)
  iter = 1
  GG= list ()
  GG[ iter ] = g(X,Y, Beta )
  # Count number of iterations not improving
  not_improving_max = 10
  not_improving_count = 0
  best_g = GG[ iter ][[1]] + 1
  best_Beta = Beta
  # Iterate until convergence is found
  while (iter <= max_iter ){
    # update Beta values
    grad = gradient (X,Y,Beta ,b)
    Beta = Beta - step * grad
    iter = iter +1
    # calculate g
    new_g = g(X,Y, Beta )
    GG[ iter ] = new_g
    #if we found better solution
    if(new_g < best_g ){
      best_Beta = Beta
      best_g = new_g
      not_improving_count = 0
    }
    else
      not_improving_count = not_improving_count + 1
    # stop if norm ( grad ) is smaller than tolerance or if we are not improving
    if( norm (grad , type ="2") < tolerance || not_improving_count == not_improving_max ){
      Beta = best_Beta
      break
    }
  }
  # Plot f over the iterations
  png( filename = pngname , width = 1280*3/4 , height = 720/2)
  plot(c(1:length(GG)), GG , xlab =" Iterations ", ylab = "g(\ u03B2 )", cex.lab =1.5 , cex.axis =1.5, title(main = paste ("g(\ u03B2 ) vs Iterations , b = ",b)))
        dev.off()
        plot(c(1:length(GG)), GG , xlab ="Iterations ", ylab = "g(\ u03B2 )")
        # Compute difference with real Beta
        data.Beta = read.csv(file ="True_Beta.csv", head =FALSE ,sep =",")
        real_Beta = as.matrix(data.Beta)
        MSE = mean ((real_Beta - Beta)^2)
        # print number of iterations
        print(iter)
        return(MSE)
}

print(paste("b= 10 - MSE :", stoch_gradient(10 ,"2 convergence10.png")))
print(paste("b= 25 - MSE :", stoch_gradient(25 ,"2 convergence25.png")))
print(paste("b= 100 - MSE :", stoch_gradient(100 ,"2 convergence100.png")))
print(paste("b= 1 - MSE :", stoch_gradient(1 ,"2 convergence1.png")))