=begin
= One dimensional root-finding and the solver classes
One-dimensional root finding algorithms can be divided into two classes, 
((|root bracketing|)) and ((|root polishing|)). The state for bracketing solvers 
is held in a (({GSL::Root::FSolver})) object. The updating procedure uses only 
function evaluations (not derivatives). The state for root polishing solvers is 
held in a (({GSL::Root::FdfSolver})) object. The updates require both the function 
and its derivative (hence the name fdf) to be supplied by the user.

Contents:
(1) ((<Solver classes|URL:roots.html#1>))
(2) ((<Methods|URL:roots.html#2>))
    (1) ((<Methods to solve equations|URL:roots.html#2.1>))
    (2) ((<GSL::Function_fdf class: Providing the function to solve|URL:roots.html#2.2>))
    (3) ((<Search Stopping Parameters|URL:roots.html#2.3>))
(3) ((<Higher-level interface|URL:roots.html#3>))
(4) ((<Example|URL:roots.html#4>))

== Solver classes

--- GSL::Root::FSolver.alloc(T)
    This creates a equation solver with a root bracketing algorithm of type ((|T|)). 
    The type ((|T|)) is given by a String or a constant,
      * (({"bisection"})) or (({GSL::Root::FSolver::BISECION}))
      * (({"falsepos"})) or (({GSL::Root::FSolver::FALSEPOS}))
      * (({"brent"})) or (({GSL::Root::FSolver::BRENT})).

    * Ex:
         include GSL::Root
         s1 = FSolver.alloc("bisection")
         s2 = FSolver.alloc("brent")
         s3 = FSolver.alloc(FSolver::BISECTION)
         s4 = FSolver.alloc(FSolver::BRENT)

--- GSL::Root::FdfSolver.alloc(T)
    This creates a derivative-based solver of type ((|T|)). 
    The type ((|T|)) is given by a String or a constant,
      * (({"newton"})) or (({GSL::Root::FdfSolver::NEWTON}))
      * (({"secant"})) or (({GSL::Root::FdfSolver::SECANT}))
      * (({"steffenson"})) or (({GSL::Root::FdfSolver::STEFFENSON})).

== Methods

--- GSL::Root::FSolver#set(f, xl, xu)
    This initialize the solver ((|self|)) to use the function ((|f|)), 
    and the initial search 
    interval ((|xl, xu|)). The function to be solved ((|f|)) is given 
    an instanse of the ((<GSL::Function|URL:function.html>)) class.

--- GSL::Root::FdfSolver#set(fdf, r)
    This initializes, or reinitializes, an existing solver ((|self|)) 
    to use the function and derivative ((|fdf|)) and the initial guess ((|r|)).
    Here ((|fdf|)) is a (({GSL::Function_fdf})) object (see below).


=== Methods to solve equations

--- GSL::Root::FSolver#iterate
--- GSL::Root::FdfSolver#iterate
    This performs a single iteration of the solver. If the iteration encounters an unexpected problem then an error code will be returned ( (({GSL::EBADFUNC})) or 
    (({GSL::EZERODIV})) ).

--- GSL::Root::FSolver#root
--- GSL::Root::FdfSolver#root
    Returns the current estimate of the root.


--- GSL::Root::FSolver#name
--- GSL::Root::FdfSolver#name
    This returns the name of the algorithm.

--- GSL::Root::FSolver#x_lower
--- GSL::Root::FSolver#x_upper
    Return the current bracketing interval for the solver.

=== GSL::Function_fdf class: Providing the function to solve
The (({FSolver})) object require an instance of the 
((<GSL::Function|URL:function.html>)) class, which is already introduced elsewhere.
The (({FdfSolver})) which uses the root-polishing algorithm requires not only
the function to solve, but also
procedures to calculate the derivatives. This is
given by the (({GSL::Function_fdf})) class.

--- GSL::Function_fdf.alloc()
--- GSL::Function_fdf.alloc(f, df)
--- GSL::Function_fdf.alloc(f, df, fdf)
    Constructors. Here ((|f, df|)) are Ruby (({Proc})) objects which return a single value.
    The option ((|fdf|)) must return an array which contain the values of the function 
    and its derivative.

--- GSL::Function_fdf#set(f, df)
--- GSL::Function_fdf#set(f, df, fdf)
    This initializes or reinitializes the (({Function_fdf})) object ((|self|)) by
    two or three (({Proc})) objects ((|f, df|)) and ((|fdf|)).
    
    * ex: A quadratic equation a*x*x + b*x + c = 0:

        # Returns a value of the function
        f = Proc.new { |x, params| 
          a = params[0]; b = params[1]; c = params[2]
          (a*x + b)*x + c
        }
        # Calculate the derivative
        df = Proc.new { |x, params| 
          a = params[0]; b = params[1]
          2*a*x + b
        }

        function_fdf = Function_fdf.alloc(f, df)

--- GSL::Function_fdf#set(f, df, params...)
--- GSL::Function_fdf#set(f, df, fdf, params...)
    This sets or resets the procedures and the constant parameters in the function.

--- GSL::Function_fdf#set_params(...)
    This sets or resets the constant parameters in the function.

    * Ex: x*x - 5 == 0
        
        function_fdf.set_params([1, 0, -5])

=== Search Stopping Parameters

--- GSL::Root::test_interval(xl, xu, epsrel, epsabs)
    This function tests for the convergence of the interval 
    ((|[xl, xu]|)) with absolute error ((|epsabs|)) and relative error ((|epsrel|)). 
    The test returns (({GSL::SUCCESS})) if the following condition is achieved,
      |a - b| < epsabs + epsrel min(|a|,|b|) 
    when the interval x = [a,b] does not include the origin. If the interval includes 
    the origin then min(|a|,|b|) is replaced by zero (which is the minimum value of 
    |x| over the interval). This ensures that the relative error is accurately 
    estimated for roots close to the origin.

    This condition on the interval also implies that any estimate of the root r in 
    the interval satisfies the same condition with respect to the true root r0,

      |r - r0| < epsabs + epsrel r0

    assuming that the true root r0 is contained within the interval.

--- GSL::Root::test_delta(x1, x0, epsrel, epsabs)
    This function tests for the convergence of the sequence ..., ((|x0, x1|)) with 
    absolute error ((|epsabs|)) and relative error ((|epsrel|)). 
    The test returns (({GSL::SUCCESS})) if the following condition is achieved,
      |x_1 - x_0| < epsabs + epsrel |x_1|
    and returns (({GSL::CONTINUE})) otherwise.

--- GSL::Root::test_residual(f, epsabs)
    This function tests the residual value ((|f|)) against the absolute error bound 
    ((|epsabs|)). The test returns (({GSL::SUCCESS})) if the following condition is 
    achieved,
      |f| < epsabs
    and returns (({GSL::CONTINUE})) otherwise. This criterion is suitable for 
    situations where the precise location of the root, x, is unimportant provided a 
    value can be found where the residual, |f(x)|, is small enough.

== High-level interface
--- GSL::Root:FSolver.solve(func, [xl, xu], [epsabs = 0, epsrel = 1e-6])
    This method try to find a root of the function ((|func|)) between the interval
    ((|[xl, xu]|)), with the accuracy ((|[epsabs, epsrel]|)) (optional). An array
    of 3 elements is returned, as [((|root, iterations, status|))].

--- GSL::Root:FdfSolver.solve(func, x0, [epsabs = 0, epsrel = 1e-6])
    This method try to find a root of the function ((|func|)) around ((|x0|)),
    with the accuracy ((|[epsabs, epsrel]|)) (optional).
    An array of 3 elements is returned, as [((|root, iterations, status|))].

--- GSL::Function#fsolve([xl, xu])
--- GSL::Function#fsolve(xl, xu)
    These methods try to find a root of (({f(x) = 0})) between the interval
    ((|[xl, xh]|)) using Brent's algorithm.
    An array of 3 elements is returned, as [((|root, iterations, status|))].

    * ex:
        f = Function.alloc { |x| x*x - 5 }
        f.fsolve([0, 5])             <----- 2.23606797749979

== Example
This example is equivalent to the one found in the GSL manual, 
using the Brent's algorithm to solve the equation x^2 - 5 = 0.

  #!/usr/bin/env ruby
  require "gsl"

  #solver = Root::FSolver.alloc("bisection")
  #solver = Root::FSolver.alloc("falsepos")
  solver = Root::FSolver.alloc(Root::FSolver::BRENT)
  puts "using #{solver.name} method"

  func = GSL::Function.alloc { |x, params|      # Define a function to solve
    a = params[0]; b = params[1]; c = params[2]
    (a*x + b)*x + c
  }
  expected = Math::sqrt(5.0)

  func.set_params([1, 0, -5])

  printf("%5s [%9s, %9s] %9s %10s %9s\n",
          "iter", "lower", "upper", "root", 
          "err", "err(est)")

  solver.set(func, 0.0, 5.0)              # initialize the solver
  i = 1
  begin
    status = solver.iterate
    r = solver.root
    xl = solver.x_lower
    xu = solver.x_upper
    status = Root.test_interval(xl, xu, 0, 0.001)   # Check convergence
    if status == GSL::SUCCESS
      printf("Converged:\n")
    end
    printf("%5d [%.7f, %.7f] %.7f %+.7f %.7f\n",
           i, xl, xu, r, r - expected, xu - xl)

    i += 1
  end while status != GSL::SUCCESS

The following is an another version, using the (({FdfSolver})) with the Newton-Raphson 
algorithm.

     #!/usr/bin/env ruby
     require "gsl"

     f = Proc.new { |x, params| 
       a = params[0]; b = params[1]; c = params[2]
       (a*x + b)*x + c
     }

     df = Proc.new { |x, params| 
       a = params[0]; b = params[1]
       2.0*a*x + b
     }

     function_fdf = Function_fdf.alloc(f, df)
     params = [1, 0, -5]
     function_fdf.set_params(params)

     solver = Root::FdfSolver.alloc(Root::FdfSolver::NEWTON)
     puts "using #{solver.name} method"

     expected = Math::sqrt(5.0)
     x = 5.0
     solver.set(function_fdf, x)

     printf("%-5s %10s %10s %10s\n", "iter", "root", "err", "err(est)")
     iter = 0
     begin
       iter += 1
       status = solver.iterate
       x0 = x
       x = solver.root

       status = Root::test_delta(x, x0, 0, 1e-3)
       
       if status == GSL::SUCCESS
         printf("Converged:\n")
       end

       printf("%5d %10.7f %+10.7f %10.7f\n", iter, x, x - expected, x - x0)
     end while status != GSL::SUCCESS

((<prev|URL:dht.html>))
((<next|URL:min.html>))

((<Reference index|URL:ref.html>))
((<top|URL:index.html>))

=end
