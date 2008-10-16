=begin
= One dimensional Minimization

This chapter describes routines for finding minima of arbitrary 
one-dimensional functions.


Contents:
(1) ((<Introduction|URL:min.html#1>))
(2) ((<GSL::Min::FMinimizer class|URL:min.html#2>))
(3) ((<Iteration|URL:min.html#3>))
(4) ((<Stopping Parameters|URL:min.html#4>))
(5) ((<Examples|URL:min.html#5>))

== Introduction

The minimization algorithms begin with a bounded region known to contain 
a minimum. The region is described by (({a})) lower bound a and an upper bound 
(({b})), with an estimate of the location of the minimum (({x})).

The value of the function at (({x})) must be less than the value of the 
function at the ends of the interval,
  f(a) > f(x) < f(b)
This condition guarantees that a minimum is contained somewhere within the 
interval. On each iteration a new point (({x'})) is selected using one of the
available algorithms. If the new point is a better estimate of the minimum, 
(({f(x') < f(x)})), then the current estimate of the minimum (({x})) is 
updated. The new point also allows the size of the bounded interval to be 
reduced, by choosing the most compact set of points which satisfies the 
constraint (({f(a) > f(x) < f(b)})). The interval is reduced until it 
encloses the true minimum to a desired tolerance. This provides a best 
estimate of the location of the minimum and a rigorous error estimate.

Several bracketing algorithms are available within a single framework. 
The user provides a high-level driver for the algorithm, and the library 
provides the individual functions necessary for each of the steps. There 
are three main phases of the iteration. The steps are,
  * initialize minimizer (or ((|solver|))) state, (({s})), for algorithm (({T}))
  * update (({s})) using the iteration (({T}))
  * test (({s})) for convergence, and repeat iteration if necessary

The state of the minimizers is held in a (({GSL::Min::FMinimizer})) object. 
The updating procedure use only function evaluations (not derivatives).
The function to minimize is given as an instance of the ((<GSL::Function|URL:function.html>)) class to the minimizer.


== GSL::Min::FMinimizer class
--- GSL::Min::FMinimizer.alloc(t)
    These method create an instance of the (({GSL::Min::FMinimizer})) class of 
    type ((|t|)). The type ((|t|)) is given by a Ruby constant,
      * GSL::Min::FMinimizer::GOLDENSECTION
      * GSL::Min::FMinimizer::BRENT
    ex1)
        include GSL
        s1 = Min::FMinimizer.alloc(Min::FMinimizer::GOLDENSECTION)
    ex2)
        include GSL::Min
        s2 = FMinimizer.alloc(FMinimizer::BRENT)

--- GSL::Min::FMinimizer#set(f, xmin, xlow, xup)
    This method sets, or resets, an existing minimizer ((|self|)) to use 
    the function ((|f|)) (given by a (({GSL::Function}))
    object) and the initial search interval [((|xlow, xup|))], 
    with a guess for the location of the minimum ((|xmin|)).

    If the interval given does not contain a minimum, then the 
    method returns an error code of (({GSL::FAILURE})).

--- GSL::Min::FMinimizer#set_with_values(f, xmin, fmin, xlow, flow, xup, fup)
    This method is equivalent to (({Fminimizer#set})) but uses the values 
    ((|fmin, flowe|)) and ((|fup|)) instead of computing 
    ((|f(xmin), f(xlow)|)) and ((|f(xup)|)).

--- GSL::Min::FMinimizer#name
    This returns the name of the minimizer. 

== Iteration
--- GSL::Min::FMinimizer#iterate
    This method performs a single iteration of the minimizer ((|self|)). 
    If the iteration encounters an unexpected problem then an error code 
    will be returned,
    * (({GSL::EBADFUNC})): the iteration encountered a singular point where the 
      function evaluated to (({Inf})) or (({NaN})).
    * (({GSL::FAILURE})): the algorithm could not improve the current best 
      approximation or bounding interval.
    The minimizer maintains a current best estimate of the position of 
    the minimum at all times, and the current interval bounding the minimum. 
    This information can be accessed with the following auxiliary methods

--- GSL::Min::FMinimizer#x_minimum
    Returns the current estimate of the position of the minimum 
    for the minimizer ((|self|)).

--- GSL::Min::FMinimizer#x_upper
--- GSL::Min::FMinimizer#x_lower
    Return the current upper and lower bound of the interval for the 
    minimizer ((|self|)).

--- GSL::Min::FMinimizer#f_minimum
--- GSL::Min::FMinimizer#f_upper
--- GSL::Min::FMinimizer#f_lower
    Return the value of the function at the current estimate of the 
    minimum and at the upper and lower bounds of interval 
    for the minimizer ((|self|)).

== Stopping Parameters
--- GSL::Min::FMinimizer#test_interval(epsabs, epsrel)
--- GSL::Min.test_interval(xlow, xup, epsabs, epsrel)
    These methoeds test for the convergence of the interval 
    [((|xlow, xup|))] with absolute error ((|epsabs|)) and relative 
    error ((|epsrel|)). The test returns (({GSL::SUCCESS})) 
    if the following condition is achieved,
      |a - b| < epsabs + epsrel min(|a|,|b|) 
    when the interval (({x = [a,b]})) does not include the origin. 
    If the interval includes the origin then (({min(|a|,|b|)})) is 
    replaced by zero (which is the minimum value of |x| over the interval). 
    This ensures that the relative error is accurately estimated for minima 
    close to the origin.

    This condition on the interval also implies that any estimate of the 
    minimum x_m in the interval satisfies the same condition with respect 
    to the true minimum x_m^*,
      |x_m - x_m^*| < epsabs + epsrel x_m^*
    assuming that the true minimum x_m^* is contained within the interval.

== Example
To find the minimum of the function f(x) = cos(x) + 1.0:

     #!/usr/bin/env ruby
     require("gsl")
     include GSL::Min

     fn1 = Function.alloc { |x| Math::cos(x) + 1.0 }

     iter = 0;  max_iter = 500
     m = 2.0             # initial guess
     m_expected = Math::PI
     a = 0.0
     b = 6.0

     gmf = FMinimizer.alloc(FMinimizer::BRENT)
     gmf.set(fn1, m, a, b)

     printf("using %s method\n", gmf.name)
     printf("%5s [%9s, %9s] %9s %10s %9s\n", "iter", "lower", "upper", "min",
            "err", "err(est)")

     printf("%5d [%.7f, %.7f] %.7f %+.7f %.7f\n", iter, a, b, m, m - m_expected, b - a)

     begin
       iter += 1
       status = gmf.iterate
       status = gmf.test_interval(0.001, 0.0)
       puts("Converged:") if status == GSL::SUCCESS
       a = gmf.x_lower
       b = gmf.x_upper
       m = gmf.x_minimum
       printf("%5d [%.7f, %.7f] %.7f %+.7f %.7f\n",
              iter, a, b, m, m - m_expected, b - a);
     end while status == GSL::CONTINUE and iter < max_iter

((<prev|URL:roots.html>))
((<next|URL:multiroot.html>))

((<Reference index|URL:ref.html>))
((<top|URL:index.html>))

=end
