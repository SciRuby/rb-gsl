=begin
= Least-Squares Fitting
This chapter describes routines for performing least squares fits to 
experimental data using linear combinations of functions. The data may be 
weighted or unweighted, i.e. with known or unknown errors. For weighted data 
the functions compute the best fit parameters and their associated covariance 
matrix. For unweighted data the covariance matrix is estimated from the 
scatter of the points, giving a variance-covariance matrix. 

The functions are divided into separate versions for simple one- or 
two-parameter regression and multiple-parameter fits. 

Contents:
(1) ((<Overview|URL:fit.html#1>))
(2) ((<Linear regression|URL:fit.html#2>))
    (1) ((<Module functions for linear regression|URL:fit.html#2.1>))
(3) ((<Linear fitting without a constant term|URL:fit.html#3>))
(4) ((<Multi-parameter fitting|URL:fit.html#4>))
    (1) ((<GSL::MultiFit::Workspace class|URL:fit.html#4.1>))
    (2) ((<Module functions|URL:fit.html#4.2>))
    (3) ((<Higer level interface|URL:fit.html#4.3>))
    (4) ((<NDLINEAR: multi-linear, multi-parameter least squares fitting|URL:ndlinear.html>)) (GSL extension)
(5) ((<Examples|URL:fit.html#5>))
    (1) ((<Linear regression|URL:fit.html#5.1>))
    (2) ((<Exponential fitting|URL:fit.html#5.2>))
    (3) ((<Multi-parameter fitting|URL:fit.html#5.3>))

== Overview
Least-squares fits are found by minimizing \chi^2 (chi-squared), the weighted 
sum of squared residuals over n experimental datapoints (x_i, y_i) for the 
model Y(c,x), The p parameters of the model are c = {c_0, c_1, Åc}. The weight 
factors w_i are given by w_i = 1/\sigma_i^2, where \sigma_i is the 
experimental error on the data-point y_i. The errors are assumed to be 
gaussian and uncorrelated. For unweighted data the chi-squared sum is computed 
without any weight factors. 

The fitting routines return the best-fit parameters c and their p \times p 
covariance matrix. The covariance matrix measures the statistical errors on 
the best-fit parameters resulting from the errors on the data, \sigma_i, and 
is defined as C_{ab} = <\delta c_a \delta c_b> where < > denotes an average 
over the gaussian error distributions of the underlying datapoints. 

The covariance matrix is calculated by error propagation from the data errors 
\sigma_i. The change in a fitted parameter \delta c_a caused by a small change 
in the data \delta y_i is given by allowing the covariance matrix to be written
in terms of the errors on the data, For uncorrelated data the fluctuations of 
the underlying datapoints satisfy 
<\delta y_i \delta y_j> = \sigma_i^2 \delta_{ij}, giving a corresponding 
parameter covariance matrix of When computing the covariance matrix for 
unweighted data, i.e. data with unknown errors, the weight factors w_i in this 
sum are replaced by the single estimate w = 1/\sigma^2, where \sigma^2 is the 
computed variance of the residuals about the 
best-fit model, \sigma^2 = \sum (y_i - Y(c,x_i))^2 / (n-p). 
This is referred to as the variance-covariance matrix. 

The standard deviations of the best-fit parameters are given by the square 
root of the corresponding diagonal elements of the covariance matrix, 
\sigma_{c_a} = \sqrt{C_{aa}}. The correlation coefficient of the fit 
parameters c_a and c_b is given by \rho_{ab} = C_{ab} / \sqrt{C_{aa} C_{bb}}. 


== Linear regression
The functions described in this section can be used to perform least-squares 
fits to a straight line model, Y = c_0 + c_1 X. For weighted data the best-fit 
is found by minimizing the weighted sum of squared residuals, chi^2,

chi^2 = sum_i w_i (y_i - (c0 + c1 x_i))^2

for the parameters (({c0, c1})). For unweighted data the sum is computed with 
(({w_i = 1})).

=== Module functions for linear regression
--- GSL::Fit::linear(x, y)
    This function computes the best-fit linear regression coefficients (c0,c1) 
    of the model Y = c0 + c1 X for the datasets ((|(x, y)|)), two vectors of 
    equal length with stride 1. This returns an array of 7 elements, 
    (({[c0, c1, cov00, cov01, cov11, chisq, status]})), where (({c0, c1})) are the
    estimated parameters, (({cov00, cov01, cov11})) are the variance-covariance 
    matrix elements, (({chisq})) is the sum of squares of the residuals, and
    (({status})) is the return code from the GSL function (({gsl_fit_linear()})).

--- GSL::Fit::wlinear(x, w, y)
    This function computes the best-fit linear regression coefficients (c0,c1) 
    of the model Y = c_0 + c_1 X for the weighted datasets ((|(x, y)|)). 
    The vector ((|w|)), specifies the weight of each datapoint, which is the 
    reciprocal of the variance for each datapoint in ((|y|)). This returns an
    array of 7 elements, same as the method (({linear})).

--- GSL::Fit::linear_est(x, c0, c1, c00, c01, c11)
--- GSL::Fit::linear_est(x, [c0, c1, c00, c01, c11])
    This function uses the best-fit linear regression coefficients ((|c0,c1|)) and 
    their estimated covariance ((|cov00,cov01,cov11|)) to compute the fitted function 
    and its standard deviation for the model Y = c_0 + c_1 X at the point ((|x|)).
    The returned value is an array of (({[y, yerr]})).

== Linear fitting without a constant term
--- GSL::Fit::mul(x, y)
    This function computes the best-fit linear regression coefficient (({c1})) 
    of the model Y = c1 X for the datasets ((|(x, y)|)), two vectors of 
    equal length with stride 1. This returns an array of 4 elements, 
    (({[c1, cov11, chisq, status]})).

--- GSL::Fit::wmul(x, w, y)
    This function computes the best-fit linear regression coefficient (({c1})) 
    of the model Y = c_1 X for the weighted datasets (({(x, y)})). The vector 
    ((|w|)) specifies the weight of each datapoint. The weight is the reciprocal 
    of the variance for each datapoint in ((|y|)).

--- GSL::Fit::mul_est(x, c1, c11)
--- GSL::Fit::mul_est(x, [c1, c11])
    This function uses the best-fit linear regression coefficient ((|c1|)) 
    and its estimated covariance ((|cov11|)) to compute the fitted function 
    (({y})) and its standard deviation (({y_err})) 
    for the model Y = c_1 X at the point ((|x|)).  
    The returned value is an array of (({[y, yerr]})).

== Multi-parameter fitting
=== GSL::MultiFit::Workspace class
--- GSL::MultiFit::Workspace.alloc(n, p)
    This creates a workspace for fitting a model to ((|n|)) 
    observations using ((|p|)) parameters.

=== Module functions
--- GSL::MultiFit::linear(X, y, work)
--- GSL::MultiFit::linear(X, y)
    This function computes the best-fit parameters (({c})) of the model (({y = X c})) 
    for the observations ((|y|)) and the matrix of predictor variables ((|X|)). 
    The variance-covariance matrix of the model parameters (({cov})) is estimated 
    from the scatter of the observations about the best-fit. The sum of squares 
    of the residuals from the best-fit is also calculated. The returned value is
    an array of 4 elements, (({[c, cov, chisq, status]})), where (({c})) is a
    ((<GSL::Vector|URL:vector.html>)) object which contains the best-fit parameters,
    and (({cov})) is the variance-covariance matrix as a 
    ((<GSL::Matrix|URL:matrix.html>)) object.

    The best-fit is found by singular value decomposition of the matrix ((|X|)) 
    using the workspace provided in ((|work|)) (optional, if not given, it is allocated
    internally). 
    The modified Golub-Reinsch SVD algorithm is used, with column scaling to improve 
    the accuracy of the singular values. Any components which have zero singular 
    value (to machine precision) are discarded from the fit.

--- GSL::MultiFit::wlinear(X, w, y, work)
--- GSL::MultiFit::wlinear(X, w, y)
    This function computes the best-fit parameters (({c})) of the model 
    (({y = X c}))  for the observations ((|y|)) and the matrix of predictor 
    variables ((|X|)).  The covariance matrix of the model parameters 
    (({cov})) is estimated from the  weighted data. The weighted sum of
    squares of the residuals from the best-fit is also calculated. 
    The returned value is an array of 4 elements, 
    (({[c: Vector, cov: Matrix, chisq: Float, status: Fixnum]})).
    The best-fit is found by singular value decomposition of the matrix ((|X|)) 
    using the workspace provided in ((|work|)) (optional). Any components 
    which have 
    zero singular value (to machine precision) are discarded from the fit.

=== Higer level interface

--- GSL::MultiFit::polyfit(x, y, order)
    Finds the coefficient of a polynomial of order ((|order|)) 
    that fits the vector data (((|x, y|))) in a least-square sense.

    Example:
      #!/usr/bin/env ruby
      require("gsl")

      x = Vector[1, 2, 3, 4, 5]
      y = Vector[5.5, 43.1, 128, 290.7, 498.4]
      # The results are stored in a polynomial "coef"
      coef, err, chisq, status = MultiFit.polyfit(x, y, 3) 

      x2 = Vector.linspace(1, 5, 20)
      graph([x, y], [x2, coef.eval(x2)], "-C -g 3 -S 4")

== Examples
=== Linear regression
     #!/usr/bin/env ruby
     require("gsl")
     include GSL::Fit

     n = 4
     x = Vector.alloc(1970, 1980, 1990, 2000)
     y = Vector.alloc(12, 11, 14, 13)
     w = Vector.alloc(0.1, 0.2, 0.3, 0.4)

     #for i in 0...n do
     #   printf("%e %e %e\n", x[i], y[i], 1.0/Math::sqrt(w[i]))
     #end

     c0, c1, cov00, cov01, cov11, chisq = wlinear(x, w, y)

     printf("# best fit: Y = %g + %g X\n", c0, c1);
     printf("# covariance matrix:\n");
     printf("# [ %g, %g\n#   %g, %g]\n", 
             cov00, cov01, cov01, cov11);
     printf("# chisq = %g\n", chisq);

=== Exponential fitting
    #!/usr/bin/env ruby
    require("gsl")

    # Create data
    r = Rng.alloc("knuthran")
    a = 2.0
    b = -1.0
    sigma = 0.01
    N = 10
    x = Vector.linspace(0, 5, N)
    y = a*Sf::exp(b*x) + sigma*r.gaussian

    # Fitting
    a2, b2, = Fit.linear(x, Sf::log(y))
    x2 = Vector.linspace(0, 5, 20)
    A = Sf::exp(a2)
    printf("Expect: a = %f, b = %f\n", a, b)
    printf("Result: a = %f, b = %f\n", A, b2)
    graph([x, y], [x2, A*Sf::exp(b2*x2)], "-C -g 3 -S 4")

=== Multi-parameter fitting
     #!/usr/bin/env ruby
     require("gsl")
     include GSL::MultiFit

     Rng.env_setup()
  
     r = GSL::Rng.alloc(Rng::DEFAULT)
     n = 19
     dim = 3
     X = Matrix.alloc(n, dim)
     y = Vector.alloc(n)
     w = Vector.alloc(n)

     a = 0.1
     for i in 0...n
       y0 = Math::exp(a)
       sigma = 0.1*y0
       val = r.gaussian(sigma)
       X.set(i, 0, 1.0)
       X.set(i, 1, a)
       X.set(i, 2, a*a)
       y[i] = y0 + val
       w[i] = 1.0/(sigma*sigma)
       #printf("%g %g %g\n", a, y[i], sigma)
       a += 0.1
     end

     c, cov, chisq, status = MultiFit.wlinear(X, w, y)

((<prev|URL:multimin.html>))
((<next|URL:nonlinearfit.html>))

((<Reference index|URL:ref.html>))
((<top|URL:index.html>))

=end
