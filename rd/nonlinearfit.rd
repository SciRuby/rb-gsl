=begin
= Nonlinear Least-Squares Fitting
This chapter describes functions for multidimensional nonlinear least-squares 
fitting. The library provides low level components for a variety of iterative 
solvers and convergence tests. These can be combined by the user to achieve 
the desired solution, with full access to the intermediate steps of the 
iteration. Each class of methods uses the same framework, so that you can 
switch between solvers at runtime without needing to recompile your program. 
Each instance of a solver keeps track of its own state, allowing the solvers 
to be used in multi-threaded programs. 

Contents:
(1) ((<Overview|URL:nonlinearfit.html#1>))
(2) ((<Initializing the Solver|URL:nonlinearfit.html#2>))
    (1) ((<GSL::MultiFit::FdfSolver class|URL:nonlinearfit.html#2.1>))
(3) ((<Providing the function to be minimized|URL:nonlinearfit.html#3>))
    (1) ((<GSL::MultiFit::Function_fdf class|URL:nonlinearfit.html#3.1>))
(4) ((<Iteration|URL:nonlinearfit.html#4>))
(5) ((<Search Stopping Parameters|URL:nonlinearfit.html#5>))
(6) ((<Computing the covariance matrix of best fit parameters|URL:nonlinearfit.html#6>))
(7) ((<Higher level interfaces|URL:nonlinearfit.html#7>))
(8) ((<Examples|URL:nonlinearfit.html#8>))
    (1) ((<Fitting to user-defined functions|URL:nonlinearfit.html#8.1>))
    (2) ((<Fitting to built-in functions|URL:nonlinearfit.html#8.2>))

== Overview
The problem of multidimensional nonlinear least-squares fitting requires the 
minimization of the squared residuals of n functions, f_i, in p parameters, 
x_i, All algorithms proceed from an initial guess using the linearization, 
where x is the initial point, p is the proposed step and J is the Jacobian 
matrix J_{ij} = d f_i / d x_j. Additional strategies are used to enlarge the 
region of convergence. These include requiring a decrease in the norm ||F|| 
on each step or using a trust region to avoid steps which fall outside the 
linear regime. 

To perform a weighted least-squares fit of a nonlinear model Y(x,t) to data 
(t_i, y_i) with independent gaussian errors \sigma_i, use function components 
of the following form, Note that the model parameters are denoted by x in this 
chapter since the non-linear least-squares algorithms are described 
geometrically (i.e. finding the minimum of a surface). The independent 
variable of any data to be fitted is denoted by t. 

With the definition above the Jacobian is 
J_{ij} =(1 / \sigma_i) d Y_i / d x_j, where Y_i = Y(x,t_i). 

== Initializing the Solver 

=== GSL::MultiFit::FdfSolver class
--- GSL::MultiFit::FdfSolver.alloc(T, n, p)
    This creates an instance of the (({GSL::MultiFit::FdfSolver})) class of 
    type ((|T|)) for ((|n|)) observations and ((|p|)) parameters. The type ((|T|))
    is given by a (({Fixnum})) constant or a (({String})), 
    * (({GSL::MultiFit::LMSDER})) or (({"lmsder"}))
    * (({GSL::MultiFit::LMDER})) or (({"lmder"}))
    For example, the following code creates an instance of a Levenberg-Marquardt 
    solver for 100 data points and 3 parameters,

        solver = MultiFit::FdfSolver.alloc(MultiFit::LMDER, 100, 3)

--- GSL::MultiFit::FdfSolver#set(f, x)
    This method initializes, or reinitializes, an existing solver ((|self|)) 
    to use the function ((|f|)) and the initial guess ((|x|)). The function ((|f|))
    is an instance of the (({GSL::MultiFit::Function_fdf})) class (see below). The
    initial guess of the parameters ((|x|)) is given by a ((<GSL::Vector|URL:vector.html>)) object.

--- GSL::MultiFit::FdfSolver#name
    This returns the name of the solver ((|self|)) as a String.


--- GSL::MultiFit::FdfSolver#x
--- GSL::MultiFit::FdfSolver#dx
--- GSL::MultiFit::FdfSolver#f
--- GSL::MultiFit::FdfSolver#J
--- GSL::MultiFit::FdfSolver#jacobian
--- GSL::MultiFit::FdfSolver#jac
    Access to the members (see (({gsl_multifit_nlin.h})))

== Providing the function to be minimized
=== GSL::MultiFit::Function_fdf class
--- GSL::MultiFit::Function_fdf.alloc()
--- GSL::MultiFit::Function_fdf.alloc(f, df, p)
--- GSL::MultiFit::Function_fdf.alloc(f, df, fdf, p)
    Constructor for the (({Function_fdf})) class, to a
    function with ((|p|)) parameters, The first two or three arguments are Ruby Proc objects 
    to evaluate the function to minimize and its derivative (Jacobian). 

--- GSL::MultiFit::Function_fdf#set_procs(f, df, p)
--- GSL::MultiFit::Function_fdf#set_procs(f, df, fdf, p)
    This initialize of reinitialize the function ((|self|)) with ((|p|)) parameters
    by two or three Proc objects ((|f, df|)) and ((|fdf|)).

--- GSL::MultiFit::Function_fdf#set_data(t, y)
--- GSL::MultiFit::Function_fdf#set_data(t, y, sigma)
    This sets the data ((|t, y, sigma|)) of length ((|n|)), to the function ((|self|)).
        
== Iteration
--- GSL::MultiFit::FdfSolver#iterate
    THis performs a single iteration of the solver ((|self|)). If the iteration 
    encounters an unexpected problem then an error code will be returned. 
    The solver maintains a current estimate of the best-fit parameters at all 
    times. This information can be accessed with the method (({position})).

--- GSL::MultiFit::FdfSolver#position
    This returns the current position (i.e. best-fit parameters) of the solver 
    ((|self|)), as a (({GSL::Vector})) object.


== Search Stopping Parameters
A minimization procedure should stop when one of the following conditions is true:
  * A minimum has been found to within the user-specified precision.
  * A user-specified maximum number of iterations has been reached.
  * An error has occurred.
The handling of these conditions is under user control. The method below allows
the user to test the current estimate of the best-fit parameters.

--- GSL::MultiFit::FdfSolver#test_delta(epsabs, epsrel)
    This method tests for the convergence of the sequence by comparing the last 
    step with the absolute error ((|epsabs|)) and relative error (((|epsrel|)) 
    to the current position. The test returns (({GSL::SUCCESS})) if the following 
    condition is achieved,
      |dx_i| < epsabs + epsrel |x_i|
    for each component of ((|x|)) and returns (({GSL::CONTINUE})) otherwise.

--- GSL::MultiFit::FdfSolver#test_gradient(g, epsabs)
--- GSL::MultiFit::FdfSolver#test_gradient(epsabs)
    This function tests the residual gradient ((|g|)) against the absolute error 
    bound ((|epsabs|)). If ((|g|)) is not given, it is calculated internally. 
    Mathematically, the gradient should be exactly zero at the minimum. 
    The test returns (({GSL::SUCCESS})) if the following condition is achieved,
      \sum_i |g_i| < epsabs
    and returns (({GSL::CONTINUE})) otherwise. This criterion is suitable for 
    situations where the precise location of the minimum, x, is unimportant provided 
    a value can be found where the gradient is small enough.

--- GSL::MultiFit::FdfSolver#gradient
    This method returns the gradient g of \Phi(x) = (1/2) ||F(x)||^2 from the 
    Jacobian matrix and the function values, using the formula g = J^T f.

--- GSL::MultiFit.test_delta(dx, x, epsabs, epsrel)
--- GSL::MultiFit.test_gradient(g, epsabs)
--- GSL::MultiFit.gradient(jac, f, g)
--- GSL::MultiFit.covar(jac, epsrel)
--- GSL::MultiFit.covar(jac, epsrel, covar)
    Singleton methods of the (({GSL::MultiFit})) module.


== Computing the covariance matrix of best fit parameters
--- GSL::MultiFit.covar(J, epsrel)
--- GSL::MultiFit.covar(J, epsrel, covar)
    This method uses the Jacobian matrix ((|J|)) to compute the covariance 
    matrix of the best-fit parameters. If an existing matrix ((|covar|)) is given,
    it is overwritten, and if not, this method returns a new matrix. 
    The parameter ((|epsrel|)) is used to remove linear-dependent columns when 
    ((|J|)) is rank deficient.

    The covariance matrix is given by,
       covar = (J^T J)^{-1}
    and is computed by QR decomposition of ((|J|)) with column-pivoting. 
    Any columns of R which satisfy
      |R_{kk}| <= epsrel |R_{11}|
    are considered linearly-dependent and are excluded from the covariance matrix 
    (the corresponding rows and columns of the covariance matrix are set to zero).

== Higher level interfaces
--- GSL::MultiFit::FdfSolver.fit(x, y, type[, guess])
--- GSL::MultiFit::FdfSolver.fit(x, w, y, type[, guess])
    This method uses (({FdfSolver})) with the LMSDER algorithm to fit the data
    ((|[x, y]|)) to a function of type ((|type|)). The returned value is 
    an array of 4 elements, (({[coef, err, chisq, dof]})),
    where (({coef})) is an array of the fitting coefficients, (({err})) contains
    errors in estimating (({coef})), (({chisq})) is the 
    chi-squared, and (({dof})) is the degree-of-freedom in the fitting 
    which equals to (data length - number of fitting coefficients). The optional
    argument ((|guess|)) is an array of initial guess of the coefficients.
    The fitting type ((|type|)) is given by a (({String})) as follows.
    * (({"gaussian"})): Gaussian fit, 
      (({y = y0 + A exp(-(x-x0)^2/2/var)})), (({coef = [y0, A, x0, var]}))
    * (({"gaussian_2peaks"})): 2-peak Gaussian fit,
      (({y = y0 + A1 exp(-(x-x1)^2/2/var1) + A2 exp(-(x-x2)^2/2/var2)})), (({coef = [y0, A1, x1, var1, A2, x2, var2]}))
    * (({"exp"})): Exponential fit, 
      (({y = y0 + A exp(-b x)})), (({coef = [y0, A, b]}))
    * (({"dblexp"})): Double exponential fit, 
      (({y = y0 + A1 exp(-b1 x) + A2 exp(-b2 x)})), (({coef = [y0, A1, b1, A2, b2]}))
    * (({"sin"})): Sinusoidal fit, 
      (({y = y0 + A sin(f x + phi)})), (({coef = [y0, A, f, phi]}))
    * (({"lor"})): Lorentzian peak fit, 
      (({y = y0 + A/((x-x0)^2 + B)})), (({coef = [y0, A, x0, B]}))
    * (({"hill"})): Hill's equation fit, 
      (({y = y0 + (m - y0)/(1 + (xhalf/x)^r)})), (({coef = [y0, n, xhalf, r]}))
    * (({"sigmoid"})): Sigmoid (Fermi-Dirac) function fit, 
      (({y = y0 + m/(1 + exp((x0-x)/r))})), (({coef = [y0, m, x0, r]}))
    * (({"power"})): Power-law fit, 
      (({y = y0 + A x^r})), (({coef = [y0, A, r]}))
    * (({"lognormal"})): Lognormal peak fit, 
      (({y = y0 + A exp[ -(log(x/x0)/width)^2 ]})), (({coef = [y0, A, x0, width]}))

    See ((<Linear fitting|URL:fit.html#2.3>)) for linear and polynomical fittings.

== Examples
=== Fitting to user-defined functions

The following example program fits a weighted exponential model with background 
to experimental data, Y = A exp(-lambda t) + b. The first part of the program sets 
up the functions (({procf})) and (({procdf})) to calculate the model and its Jacobian. 
The appropriate fitting function is given by,
  f_i = ((A exp(-lambda t_i) + b) - y_i)/sigma_i
where we have chosen t_i = i. The Jacobian matrix (({jac})) is the derivative of 
these functions with respect to the three parameters (A, lambda, b). It is given by,
  J_{ij} = d f_i / d x_j
where x_0 = A, x_1 = lambda and x_2 = b.


     require("gsl")
     include GSL::MultiFit

     # x: Vector, list of the parameters to determine
     # t, y, sigma: Vectors, observational data
     # f: Vector, function to minimize
     procf = Proc.new { |x, t, y, sigma, f|
       a = x[0]
       lambda = x[1]
       b = x[2]
       n = t.size
       for i in 0...n do
         yi = a*Math::exp(-lambda*t[i]) + b
         f[i] = (yi - y[i])/sigma[i]
       end
     }

     # jac: Matrix, Jacobian
     procdf = Proc.new { |x, t, y, sigma, jac|
       a = x[0]
       lambda = x[1]
       n = t.size
       for i in 0...n do
         ti = t[i]
         si = sigma[i]
         ei = Math::exp(-lambda*ti)
         jac.set(i, 0, ei/si)
         jac.set(i, 1, -ti*a*ei/si)
         jac.set(i, 2, 1.0/si)
       end
     }

     f = GSL::MultiFit::Function_fdf.alloc(procf, procdf, 2)

     # Create data
     r = GSL::Rng.alloc()
     t = GSL::Vector.alloc(n)
     y = GSL::Vector.alloc(n)
     sigma = Vector.alloc(n)
     for i in 0...n do
       t[i] = i
       y[i] = 1.0 + 5*Math::exp(-0.1*t[i]) + r.gaussian(0.1)
       sigma[i] = 0.1
     end

     f.set_data(t, y, sigma)
     x = GSL::Vector.alloc(1.0, 0.0, 0.0)    # initial guess

     solver = GSL::FdfSolver.alloc(FdfSolver::LMSDER, n, np)

     solver.set(f, x)

     iter = 0
     solver.print_state(iter)
     begin
       iter += 1
       status = solver.iterate
       solver.print_state(iter)
       status = solver.test_delta(1e-4, 1e-4)
     end while status == GSL::CONTINUE and iter < 500

     covar = solver.covar(0.0)
     position = solver.position
     chi2 = pow_2(solver.f.dnrm2)
     dof = n - np
     printf("A      = %.5f +/- %.5f\n", position[0], Math::sqrt(chi2/dof*covar[0][0]))
     printf("lambda = %.5f +/- %.5f\n", position[1], Math::sqrt(chi2/dof*covar[1][1]))
     printf("b      = %.5f +/- %.5f\n", position[2], Math::sqrt(chi2/dof*covar[2][2]))


=== Fitting to built-in functions
  #!/usr/bin/env ruby
  require("gsl")
  include MultiFit

  N = 100

  y0 = 1.0
  A = 2.0
  x0 = 3.0
  w = 0.5

  r = Rng.alloc
  x = Vector.linspace(0.01, 10, N)
  sig = 1
  # Lognormal function with noise
  y =  y0 + A*Sf::exp(-pow_2(Sf::log(x/x0)/w)) + 0.1*Ran::gaussian(r, sig, N)

  guess = [0, 3, 2, 1]
  coef, err, chi2, dof = MultiFit::FdfSolver.fit(x, y, "lognormal", guess)
  y0 = coef[0]
  amp = coef[1]
  x0 = coef[2]
  w = coef[3]

  graph(x, y, y0+amp*Sf::exp(-pow_2(Sf::log(x/x0)/w)))


((<prev|URL:fit.html>))
((<next|URL:bspline.html>))

((<Reference index|URL:ref.html>))
((<top|URL:index.html>))


=end
