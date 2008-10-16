=begin
= Multidimensional Root-Finding
This chapter describes functions for multidimensional root-finding 
(solving nonlinear systems with n equations in n unknowns). 
The library provides low level components for a variety of iterative solvers 
and convergence tests. These can be combined by the user to achieve the 
desired solution, with full access to the intermediate steps of the iteration. 
Each class of methods uses the same framework, so that you can switch between 
solvers at runtime without needing to recompile your program. Each instance of 
a solver keeps track of its own state, allowing the solvers to be used in 
multi-threaded programs.

(1) ((<Overview|URL:multiroot.html#1>))
(2) ((<Initializing the Solver|URL:multiroot.html#2>))
(3) ((<Providing the function to solve|URL:multiroot.html#3>))
(4) ((<Iteration|URL:multiroot.html#4>))
(5) ((<Search Stopping Parameters|URL:multiroot.html#5>))
(6) ((<Higher Level Interface|URL:multiroot.html#6>))
(7) ((<Examples|URL:multiroot.html#7>))
    (1) ((<FSolver|URL:multiroot.html#7.1>))
    (2) ((<FdfSolver|URL:multiroot.html#7.2>))


== Overview
The problem of multidimensional root finding requires the simultaneous 
solution of n equations, f_i, in n variables, x_i, In general there are no 
bracketing methods available for n dimensional systems, and no way of knowing 
whether any solutions exist. All algorithms proceed from an initial guess 
using a variant of the Newton iteration, where x, f are vector quantities and 
J is the Jacobian matrix J_{ij} = d f_i / d x_j. Additional strategies can be 
used to enlarge the region of convergence. These include requiring a decrease 
in the norm |f| on each step proposed by Newton's method, or taking 
steepest-descent steps in the direction of the negative gradient of |f|. 

Several root-finding algorithms are available within a single framework. 
The user provides a high-level driver for the algorithms, and the library 
provides the individual functions necessary for each of the steps. There are 
three main phases of the iteration. The steps are, 

* initialize solver state, ((|s|)), for algorithm ((|T|))
* update ((|s|)) using the iteration ((|T|))
* test ((|s|)) for convergence, and repeat iteration if necessary 

The evaluation of the Jacobian matrix can be problematic, either because 
programming the derivatives is intractable or because computation of the n^2 
terms of the matrix becomes too expensive. For these reasons the algorithms 
provided by the library are divided into two classes according to whether 
the derivatives are available or not. 

The state for solvers with an analytic Jacobian matrix is held in a 
(({GSL::MultiRoot::FdfSolver})) object. The updating procedure requires both 
the function and its derivatives to be supplied by the user. 

The state for solvers which do not use an analytic Jacobian matrix is held in 
a (({GSL::MultiRoot::FSolver})) object. The updating procedure uses only 
function evaluations (not derivatives). The algorithms estimate the matrix J 
or J^{-1} by approximate methods. 


== Initializing the Solver
Two types of solvers are available. The solver itself depends only on the 
dimension of the problem and the algorithm and can be reused for different problems.
The (({FdfSolver})) requires derivatives of the function to solve.

--- GSL::MultiRoot::FSolver.alloc(T, n)
    This creates an instance of the (({FSolver})) class of type ((|T|)) 
    for a system of ((|n|)) dimensions. The type is given by a constant or a string,
      * GSL::MultiRoot:FSolver::HYBRIDS, or "hybrids"
      * GSL::MultiRoot:FSolver::HYBRID, or "hybrid"
      * GSL::MultiRoot:FSolver::DNEWTON, or "dnewton"
      * GSL::MultiRoot:FSolver::BROYDEN, or "broyden"

--- GSL::MultiRoot::FdfSolver.alloc(T, n)
    This creates an instance of the (({FdfSolver})) class of type ((|T|)) 
    for a system of ((|n|)) dimensions. The type is given by a constant,
      * GSL::MultiRoot:FdfSolver::HYBRIDSJ, or "hybridsj"
      * GSL::MultiRoot:FdfSolver::HYBRIDJ, or "hybridj",
      * GSL::MultiRoot:FdfSolver::NEWTON, or "newton",
      * GSL::MultiRoot:FdfSolver::GNEWTON, or "gnewton

--- GSL::MultiRoot::FSolver#set(func, x)
    This method sets, or resets, an existing solver ((|self|)) 
    to use the function ((|func|)) and the initial guess ((|x|)).
    Here ((|x|)) is a (({Vector})), and ((|func|)) 
    is a (({MultiRoot:Function})) object.
--- GSL::MultiRoot::FdfSolver#set(func_fdf, x)
    This method sets, or resets, an existing solver ((|self|)) 
    to use the function ((|func_fdf|)) and the initial guess ((|x|)).
    Here ((|x|)) is a (({Vector})), and ((|func_fdf|)) 
    is a (({MultiRoot:Function_fdf})) object.

--- GSL::MultiRoot::FSolver#name
--- GSL::MultiRoot::FdfSolver#name

== Providing the function to solve
--- GSL::MultiRoot:Function.alloc(proc, dim, params)
    See example below:

      # x: vector, current guess
      # params: a scalar or an array
      # f: vector, function value
      proc = Proc.new { |x, params, f|
        a = params[0]; b = params[1]
        x0 = x[0]; x1 = x[1]
        f[0] = a*(1 - x0)
        f[1] = b*(x1 - x0*x0)
      }

      params = [1.0, 10.0]
      func = MultiRoot::Function.alloc(proc, 2, params)
      fsolver = MultiRoot::FSolver.alloc("broyden", 2)
      x = [-10, -5]    # initial guess
      fsolver.set(func, x)

--- GSL::MultiRoot:Function_fdf.alloc(proc, dim, params)
    See the example below:
    
      procf = Proc.new { |x, params, f|
        a = params[0]; b = params[1]
        x0 = x[0]; x1 = x[1]
        f[0] = a*(1 - x0)
        f[1] = b*(x1 - x0*x0)
      }

      procdf = Proc.new { |x, params, jac|
        a = params[0]; b = params[1]
        jac.set(0, 0, -a)
        jac.set(0, 1, 0)
        jac.set(1, 0, -2*b*x[0])
        jac.set(1, 1, b)
      }

      params = [1.0, 10.0]
      func_fdf = MultiRoot::Function_fdf.alloc(procf, procdf, n, params)

      fdfsolver = MultiRoot::FdfSolver.alloc("gnewton", n)
      x = [-10.0, -5.0]
      fdfsolver.set(func_fdf, x)

== Iteration
--- GSL::MultiRoot::FSolver#interate
--- GSL::MultiRoot::FdfSolver#interate
    These methods perform a single iteration of the solver ((|self|)). 
    If the iteration encounters an unexpected problem then an error code will 
    be returned,
    * GSL_EBADFUNC: the iteration encountered a singular point where the function 
      or its derivative evaluated to Inf or NaN.
    * GSL_ENOPROG: the iteration is not making any progress, preventing the 
      algorithm from continuing.
    The solver maintains a current best estimate of the root at all times. 
    This information can be accessed with the following auxiliary methods.

--- GSL::MultiRoot::FSolver#root
--- GSL::MultiRoot::FdfSolver#root
    These methods return the current estimate of the root (Vector) for the solver ((|self|)).

--- GSL::MultiRoot::FSolver#f
--- GSL::MultiRoot::FdfSolver#f
    These methds return the function value (({f(x)})) (Vector) at the current estimate 
    of the root for the solver ((|self|)).

--- GSL::MultiRoot::FSolver#dx
--- GSL::MultiRoot::FdfSolver#dx
    These method return the last step ((|dx|)) (Vector) taken by the solver ((|self|)).

== Search Stopping Parameters
--- GSL::MultiRoot::FSolver#test_delta(epsabs, epsrel)
--- GSL::MultiRoot::FdfSolver#test_delta(epsabs, epsrel)
    This method tests for the convergence of the sequence by comparing the last step 
    (({dx})) with the absolute error ((|epsabs|)) and relative error ((|epsrel|)) 
    to the current position (({x})). 
    The test returns (({GSL::SUCCESS})) if the following condition is achieved,
      |dx_i| < epsabs + epsrel |x_i|
    for each component of (({x})) and returns (({GSL::CONTINUE})) otherwise.

--- GSL::MultiRoot::FSolver#test_residual(epsabs)
--- GSL::MultiRoot::FdfSolver#test_residual(epsabs)
    This method tests the residual value (({f})) against the absolute error 
    bound ((|epsabs|)). The test returns (({GSL::SUCCESS})) if the following 
    condition is achieved,
       sum_i |f_i| < epsabs
    and returns (({GSL::CONTINUE})) otherwise. This criterion is suitable for 
    situations where the precise location of the root, (({x})), is unimportant 
    provided a value can be found where the residual is small enough.

== Higher Level Interface
--- GSL::MultiRoot::Function#solve(x0, max_iter = 1000, eps = 1e-7, type = "hybrids")
--- GSL::MultiRoot::FSolver#solve(max_iter = 1000, eps = 1e-7)
--- GSL::MultiRoot::FSolver.solve(fsolver, max_iter = 1000, eps = 1e-7)
See sample script (({examples/multiroot/fsolver3.rb})).

== Example

=== FSolver

     proc = Proc.new { |x, params, f|
       a = params[0];  b = params[1]
       x0 = x[0];  x1 = x[1]
       f[0] = a*(1 - x0)
       f[1] = b*(x1 - x0*x0)
     }

     params = [1.0, 10.0]
     func = MultiRoot::Function.alloc(proc, 2, params)

     fsolver = MultiRoot::FSolver.alloc("hybrid", 2)
     x = [-10, -5]
     fsolver.set(func, x)

     iter = 0
     begin
       iter += 1
       status = fsolver.iterate
       root = fsolver.root
       f = fsolver.f
       printf("iter = %3u x = % .3f % .3f f(x) = % .3e % .3e\n",
               iter, root[0], root[1], f[0], f[1])
       status = fsolver.test_residual(1e-7)
     end while status == GSL::CONTINUE and iter < 1000

=== FdfSolver
     n = 2

     procf = Proc.new { |x, params, f|
       a = params[0]; b = params[1]
       x0 = x[0]; x1 = x[1]
       f[0] = a*(1 - x0)
       f[1] = b*(x1 - x0*x0)
     }

     procdf = Proc.new { |x, params, jac|
       a = params[0]; b = params[1]
       jac.set(0, 0, -a)
       jac.set(0, 1, 0)
       jac.set(1, 0, -2*b*x[0])
       jac.set(1, 1, b)
     }

     params = [1.0, 10.0]
     f = MultiRoot::Function_fdf.alloc(procf, procdf, n, params)

     fdfsolver = MultiRoot::FdfSolver.alloc("gnewton", n)

     x = [-10.0, -5.0]

     fdfsolver.set(f, x)

     iter = 0
     begin
       iter += 1
       status = fdfsolver.iterate
       root = fdfsolver.root
       f = fdfsolver.f
       printf("iter = %3u x = % .3f % .3f f(x) = % .3e % .3e\n",
               iter, root[0], root[1], f[0], f[1])
       status = fdfsolver.test_residual(1e-7)
     end while status == GSL::CONTINUE and iter < 1000

((<prev|URL:min.html>))
((<next|URL:multimin.html>))

((<Reference index|URL:ref.html>))
((<top|URL:index.html>))

=end
