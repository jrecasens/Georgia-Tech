R code for all problems can be found in the Annex.

#### Gradient Descent for Multiple Linear Regression

Please write down the gradient $\nabla f(\boldsymbol{\beta})$.

**Solution**
$$\begin{aligned}
            f(\boldsymbol{\beta})
            &=  \frac{1}{2n} \| \boldsymbol{Y} - \boldsymbol{X}\boldsymbol{\beta}\|_2^2
            \\
            &=  \frac{1}{2n} (\boldsymbol{Y} - \boldsymbol{X}\boldsymbol{\beta})^\top (\boldsymbol{Y} - \boldsymbol{X}\boldsymbol{\beta})
            \\
            &=  \frac{1}{2n} \Big(
                        \boldsymbol{Y}^\top \boldsymbol{Y} 
                        - \boldsymbol{Y}^\top \boldsymbol{X}\boldsymbol{\beta}
                        - \boldsymbol{Y}^\top \boldsymbol{X}\boldsymbol{\beta} 
                        + \boldsymbol{\beta}^\top \boldsymbol{X}^\top \boldsymbol{X} \boldsymbol{\beta} 
                        \Big)
            \\
            &=           \frac{1}{2n} \boldsymbol{Y}^\top \boldsymbol{Y} 
                        - \frac{1}{n} \boldsymbol{Y}^\top \boldsymbol{X}\boldsymbol{\beta}
                        + \frac{1}{2n} \boldsymbol{\beta}^\top \boldsymbol{X}^\top \boldsymbol{X} \boldsymbol{\beta}        
            \end{aligned}$$

            
Derivative with respect to $\boldsymbol{\beta}$:
$$\begin{split}
            \nabla f(\boldsymbol{\beta})
            &= - \frac{1}{n} \boldsymbol{Y}^\top \boldsymbol{X} + \frac{1}{n} \boldsymbol{X}^\top \boldsymbol{X} \boldsymbol{\beta}
            \\
            &=  \frac{1}{n} (\boldsymbol{X}^\top \boldsymbol{X} \boldsymbol{\beta} - \boldsymbol{X}^\top \boldsymbol{Y})
            \\
            &=  \frac{\boldsymbol{X}^\top (\boldsymbol{X} \boldsymbol{\beta} - \boldsymbol{Y})  }{n} 
            \end{split}$$

Please write down the step size you use in every iteration and explain
why you use it.

**Solution**

$f(\boldsymbol{\beta})$ is twice differentiable and its Hessian
$\nabla^2 f(\boldsymbol{\beta})$ is
$\boldsymbol{X}^\top\boldsymbol{X} / n$, which does not depend on
$\boldsymbol{\beta}$. Because this cost function is convex we know that
its Hessian is positive semidefinite, so we have that (using the
definition of greatest curvature):

$$\begin{aligned}
\boldsymbol{0} \leq  \nabla^2 f(\boldsymbol{\beta}) \leq \frac 1n \|\boldsymbol{X}\|_2 \boldsymbol{I}
\end{aligned}$$

Therefore, the smallest Lipschitz constant of $\nabla f$ is the largest
eigenvalue of $\frac 1n \boldsymbol{X}^\top\boldsymbol{X}$. We know that
matrix 2-norm induced by the euclidean vector norm is:

$$\begin{aligned}
\|\boldsymbol{X}\|_2
&= \underset{\|\boldsymbol{y}\|_2 = 1}{\text{minimize}}
& \|\boldsymbol{Ay}\|_2\\
&= \sqrt{{\lambda}_{max}}  \\
\end{aligned}$$

Where ${\lambda}_{max}$ is the largest number $\lambda$ such that
$\boldsymbol{X}^\top\boldsymbol{X} - \lambda\boldsymbol{I}$ is singular.
i.e. ${\lambda}_{max}$ is the largest eigenvalue of
$\boldsymbol{X}^\top\boldsymbol{X}$.

Since we want to take the biggest steps possible, we can compute the
Lipschitz constant as $L = \frac 1n \|\boldsymbol{X}\|_2 ^ 2$ and set a
fixed step size: $$\begin{split}
\alpha
&= \frac 1L =  \frac{n}{\|\boldsymbol{X}\|_2 ^ 2}
\end{split}$$

Please explain the rule of stopping your gradient descent algorithm.

**Solution**

In a descent method, as each new point is generated by the algorithm,
the corresponding value of the objective function decreases in value.
The gradient varies as the search proceeds, tending to zero as we
approach the minimizer.
The stopping rules that we will use:

-   Condition $ \nabla f(\boldsymbol{{\beta}}^{(k+1)}) = 0$. However,
    this condition is not directly suitable as a practical stopping
    criterion because the numerical computation of the gradient will
    rarely be identically equal to zero. A practical criterion is to
    check if the norm $\|\nabla f(\boldsymbol{{\beta}}^{(k+1)})\|_2$ is
    less than a pre-specified threshold.

-   Also we will compute
    $\|f(\boldsymbol{{\beta}}^{(k+1)}) - f(\boldsymbol{{\beta}}^{(k)}\|$,
    and if the difference is less than some threshold, then we stop.

-   Finally, in some cases, to halt numerical algorithms a pre-specified
    number of iterations need to be specified, so if the above are not
    met the algorithm will stop after a finite number of iterations.

Please draw a plot of $f(\boldsymbol{\beta})$ versus number of
iterations to demonstrate the convergence of your algorithm.

**Solution**
After running 50 iterations we plot $f(\beta)$ versus number of
iterations. We can see that the convergence of our algorithm is very
fast:

![image](../R/P1d_Plot)

Please compare the result $\beta^*$ returned by your algorithm with the
true $\beta$ by computing the mean squared error
$\|\beta^* - \beta \|_2^2/30$.

**Solution**
The mean squared error (MSE) is:

$MSE = \|\beta^* - \beta \|_2^2/30 = 0.0000409143$.

This is a very small MSE which was reached very fast. Our step size
choice was very good.

#### Stochastic Gradient Descent for Multiple Linear Regression


Please draw a plot for the value of objective function $g(\beta)$ versus
the number of iterations to demonstrate the convergence.

**Solution**
The one statistical unit gradient is:
$\nabla g_i(\boldsymbol{\beta}) = (\boldsymbol{x}_i^\top \boldsymbol{\beta} - y_i) \boldsymbol{x}_i$
 

The plot of objective function $g(\beta)$ versus the number of
iterations is:

![image](../R/P2_Plot1)

We can see that with almost 1000 iterations the algorithm keeps looking
for an optimal value of the loss function. Convergence is not fast and
is not apparent that it has reached the minimum.

Please draw three plots with $b=10,25,100$ for the value of objective
function $g(\beta)$ versus the number of iterations to demonstrate the
convergence.

**Solution**

We obtain the following three plots for $b=10,25,100$:

![image](../R/P2_Plot10) 
![image](../R/P2_Plot25)
![image](../R/P2_Plot100)

We can see that as $b$ increases the algorithm converges faster. Almost
250 iterations for $b=10$, 100 iterations for $b=25$ and about 20
iterations for $b=100$.

Please compare the result $\beta^*$ returned by your algorithm with the
true $\beta$ by computing the mean squared error
$\|\beta^* - \beta \|_2^2/30$.

**Solution**

The MSE for each $b$ are:

Algorithm & MSE

Stochastic Gradient Descent (b=1) : 0.0021607

Mini-batch Stochastic Gradient Descent (b=10) : 0.0000590

Mini-batch Stochastic Gradient Descent (b=25) : 0.0000656

Mini-batch Stochastic Gradient Descent (b=100) : 0.0000310

The larger the $b$ the smaller is the MSE. Is not only converging faster
but we get better results for larger batches.

#### Online Principal Component Analysis.


Please draw a plot for $dist(w_i,v)$ versus number of iterations, where
$w_i$ is the output of Oja’s algorithm at $i$th iteration and $v$ is the
true top eigenvector of $\Sigma$.

**Solution**

Oja’s algorithm with step size $\eta_i = 1/100$ constant for all
iteration seems to have a fast convergence but it has a jagged looking
plot, meaning that with each constant step it seems to get away and
closer at each iteration.

The final similarity measure $dist(w_n,v)$ is 0.007513874.

![image](../R/P3_Plota)

Setting decreasing step sizes $\eta_i = 1/(100+i)$ at $i$th iteration.
Please draw a plot for $dist(w_i,v)$ versus number of iterations.

**Solution** 
When we use a variable step size $\eta_i = 1/(100+i)$ at each $i$th
iteration the final similarity measure $dist(w_n,v)$ is 0.0006637931.
Almost 10 times better accuracy than the previous result. Even though
the convergence is slower, we have a smoother plot that converges to a
more precise eigenvector $v$.

![image](../R/P3_Plotb)