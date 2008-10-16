=begin
= Multidimensional Minimization
This chapter describes routines for finding minima of arbitrary 
multidimensional functions. The library provides low level components for a 
variety of iterative minimizers and convergence tests. These can be combined 
by the user to achieve the desired solution, while providing full access to 
the intermediate steps of the algorithms. Each class of methods uses the 
same framework, so that you can switch between minimizers at runtime without 
needing to recompile your program. Each instance of a minimizer keeps track 
of its own state, allowing the minimizers to be used in multi-threaded 
programs. 

Contents:
(1) ((<Overview|URL:multimin.html#1>))
(2) ((<Caveats|URL:multimin.html#2>))
(3) ((<Initializing the Multidimensional Minimizer|URL:multimin.html#3>))
(4) ((<Providing a function to minimize|URL:multimin.html#4>))
(5) ((<Iteration|URL:multimin.html#5>))
(6) ((<Stopping Criteria|URL:multimin.html#6>))
(7) ((<Examples|URL:multimin.html#7>))
    (1) ((<FdfMinimizer|URL:multimin.html#7.1>))
    (2) ((<FMinimizer|URL:multimin.html#7.2>))

== Overview
The problem of multidimensional minimization requires finding a point x such 
that the scalar function, takes a value which is lower than at any neighboring 
point. For smooth functions the gradient g = \nabla f vanishes at the minimum. 
In general there are no bracketing methods available for the minimization of 
n-dimensional functions. The algorithms proceed from an initial guess using a 
search algorithm which attempts to move in a downhill direction. 

Algorithms making use of the gradient of the function perform a 
one-dimensional line minimisation along this direction until the lowest point 
is found to a suitable tolerance. The search direction is then updated with 
local information from the function and its derivatives, and the whole process 
repeated until the true n-dimensional minimum is found. 

The Nelder-Mead Simplex algorithm applies a different strategy. It maintains 
n+1 trial parameter vectors as the vertices of a n-dimensional simplex. 
In each iteration step it tries to improve the worst vertex by a simple 
geometrical transformation until the size of the simplex falls below a given 
tolerance. 

Both types of algorithms use a standard framework. The user provides a 
high-level driver for the algorithms, and the library provides the individual 
functions necessary for each of the steps. There are three main phases of the 
iteration. The steps are, 

* initialize minimizer state, s, for algorithm T 
* update s using the iteration T 
* test s for convergence, and repeat iteration if necessary 

Each iteration step consists either of an improvement to the line-minimisation 
in the current direction or an update to the search direction itself. The 
state for the minimizers is held in a (({GSL::MultiMin::FdfMinimizer})) or 
a (({GSL::MultiMin::FMinimizer})) object.

== Caveats
Note that the minimization algorithms can only search for one local minimum 
at a time. When there are several local minima in the search area, the first 
minimum to be found will be returned; however it is difficult to predict which 
of the minima this will be. In most cases, no error will be reported if you 
try to find a local minimum in an area where there is more than one. 

It is also important to note that the minimization algorithms find local 
minima; there is no way to determine whether a minimum is a global minimum of 
the function in question. 


== Initializing the Multidimensional Minimizer
--- GSL::MultiMin::FdfMinimizer.alloc(type, n)
--- GSL::MultiMin::FMinimizer.alloc(type, n)
    These method create a minimizer of type ((|type|)) for an ((|n|))-dimension function. 
    The type is given by a name of the algorithm, or by a Ruby constant.

    * (({GSL::MultiMin::FdfMinimizer::CONJUGATE_FR})) or (({"conjugate_fr"}))
    * (({GSL::MultiMin::FdfMinimizer::CONJUGATE_PR})) or (({"conjugate_pr"}))
    * (({GSL::MultiMin::FdfMinimizer::VECTOR_BFGS})) or (({"vector_bfgs"}))
    * (({GSL::MultiMin::FdfMinimizer::VECTOR_BFGS2})) or (({"vector_bfgs2"})) (GSL-1.9 or later)
    * (({GSL::MultiMin::FdfMinimizer::STEEPEST_DESCENT})) or (({"steepest_descent"}))
    * (({GSL::MultiMin::FMinimizer::NMSIMPLEX})) or (({"nmsimplex"}))

    * ex:
        include GSL::MultiMin
        m1 = FdfMinimizer.alloc(FdfMinimizer::CONJUGATE_FR, 2)
        m2 = FdfMinimizer.alloc("steepest_descent", 4)
        m3 = FMinimizer.alloc(FMinimizer::NMSIMPLEX, 3)
        m4 = FMinimizer.alloc("nmsimplex", 2)

--- GSL::MultiMin::FdfMinimizer#set(func, x, step_size, tol)
    This method initializes the minimizer ((|self|)) to minimize the function 
    ((|fdf|)) (the (({GSL::MultiMin::Function_fdf})) class, see below) starting from 
    the initial point ((|x|)) ((({GSL::Vector}))). The size of the first trial step is 
    given by ((|step_size|)) ((({Vector}))). The accuracy of the line minimization is 
    specified by ((|tol|)).

--- GSL::MultiMin::FMinimizer#set(func, x, step_size)
    This method initializes the minimizer ((|self|)) to minimize the function ((|func|)), 
    starting from the initial point ((|x|)) (Vector). The size of the initial trial steps 
    is given in vector ((|step_size|)).

--- GSL::MultiMin::FdfMinimizer#name
--- GSL::MultiMin::FMinimizer#name
    These return the name of the minimizer ((|self|)).

== Providing a function to minimize
You must provide a parametric function of (({n})) variables for the minimizers to 
operate on. You may also need to provide a routine which calculates the gradient of the 
function. In order to allow for general parameters the functions are defined by the 
classes, (({GSL::MultiMin::Function_fdf})) and (({GSL::MultiMin::Function})).

--- GSL::MultiMin:Function_fdf.alloc(proc_f, proc_df, n)
--- GSL::MultiMin:Function_fdf.alloc(proc_f, proc_df, proc_fdf, n)
--- GSL::MultiMin:Function_fdf#set_procs(proc_f, proc_df)
--- GSL::MultiMin:Function_fdf#set_procs(proc_f, proc_df, n)
--- GSL::MultiMin:Function_fdf#set_procs(proc_f, proc_df, proc_fdf, n)
--- GSL::MultiMin:Function_fdf#set_params(params)
    See example below.

      include GSL::MultiMin

      my_f = Proc.new { |v, params|
        x = v[0]; y = v[1]
        p0 = params[0]; p1 = params[1]
        10.0*(x - p0)*(x - p0) + 20.0*(y - p1)*(y - p1) + 30.0
      }

      my_df = Proc.new { |v, params, df|
        x = v[0]; y = v[1]
        p0 = params[0]; p1 = params[1]
        df[0] = 20.0*(x-p0)
        df[1] = 40.0*(y-p1)
      }

      my_func = Function_fdf.alloc(my_f, my_df, 2)
      my_func.set_params([1.0, 2.0])      # parameters

--- GSL::MultiMin:Function.alloc(proc_f, n)
--- GSL::MultiMin:Function#set_proc(proc_f)
--- GSL::MultiMin:Function#set_proc(proc_f, n)
--- GSL::MultiMin:Function#set_params(params)
    See example below.

      include GSL::MultiMin

      np = 2
      my_f = Proc.alloc { |v, params|
        x = v[0]; y = v[1]
        p0 = params[0]; p1 = params[1]
        10.0*(x - p0)*(x - p0) + 20.0*(y - p1)*(y - p1) + 30.0
      }

      my_func = Function.alloc(my_f, np)
      my_func.set_params([1.0, 2.0])      # parameters

== Iteration
--- GSL::MultiMin::FdfMinimizer#iterate
--- GSL::MultiMin::FMinimizer#iterate
    These methods perform a single iteration of the minimizer ((|self|)). 
    If the iteration encounters an unexpected problem then an error code will be returned.
    The minimizer maintains a current best estimate of the minimum at all times. 
    This information can be accessed with the following methods,

--- GSL::MultiMin::FdfMinimizer#x
--- GSL::MultiMin::FdfMinimizer#minimum
--- GSL::MultiMin::FdfMinimizer#gradient
--- GSL::MultiMin::FMinimizer#x
--- GSL::MultiMin::FMinimizer#minimum
--- GSL::MultiMin::FMinimizer#size
    These method return the current best estimate of the location of the minimum, 
    the value of the function at that point, its gradient, and minimizer specific 
    characteristic size for the minimizer ((|self|)).

--- GSL::MultiMin::FdfMinimizer#restart
    This method resets the minimizer ((|self|)) to use the current point as a new 
    starting point.

== Stopping Criteria
A minimization procedure should stop when one of the following conditions is true:
  * A minimum has been found to within the user-specified precision.
  * A user-specified maximum number of iterations has been reached.
  * An error has occurred.
The handling of these conditions is under user control. The methods below allow the 
user to test the precision of the current result.

--- GSL::MultiMin::FdfMinimizer#test_gradient(epsabs)
--- GSL::MultiMin::FdfMinimizer.test_gradient(g, epsabs)
    These method test the norm of the gradient ((|g|)) against the absolute tolerance 
    ((|epsabs|)). The gradient of a multidimensional function goes to zero at a minimum. 
    The tests return (({GSL::SUCCESS})) if the following condition is achieved,
      |g| < epsabs
    and returns (({GSL::CONTINUE})) otherwise. A suitable choice of ((|epsabs|)) can 
    be made from the desired accuracy in the function for small variations in (({x})). 
    The relationship between these quantities is given by (({\delta f = g \delta x})).

--- GSL::MultiMin::FdfMinimizer#test_size(epsabs)
--- GSL::MultiMin::FdfMinimizer.test_size(size, epsabs)
    These method test the minimizer specific characteristic ((|size|)) 
    (if applicable to the used minimizer) against absolute tolerance ((|epsabs|)). 
    The tests return ((({GSL::SUCCESS})) if the size is smaller than tolerance, 
    otherwise (({GSL::CONTINUE})) is returned.

== Examples

=== FdfMinimizer
     #!/usr/bin/env ruby
     require("gsl")
     include GSL::MultiMin

     my_f = Proc.new { |v, params|
       x = v[0]; y = v[1]
       p0 = params[0]; p1 = params[1]
       10.0*(x - p0)*(x - p0) + 20.0*(y - p1)*(y - p1) + 30.0
     }

     my_df = Proc.new { |v, params, df|
       x = v[0]; y = v[1]
       p0 = params[0]; p1 = params[1]
       df[0] = 20.0*(x-p0)
       df[1] = 40.0*(y-p1)
     }

     my_func = Function_fdf.alloc(my_f, my_df, 2)
     my_func.set_params([1.0, 2.0])      # parameters

     x = Vector.alloc(5.0, 7.0)          # starting point

     minimizer = FdfMinimizer.alloc("conjugate_fr", 2)
     minimizer.set(my_func, x, 0.01, 1e-4)

     iter = 0
     begin
       iter += 1
       status = minimizer.iterate()
       status = minimizer.test_gradient(1e-3)
       if status == GSL::SUCCESS
         puts("Minimum found at")
       end
       x = minimizer.x
       f = minimizer.f
       printf("%5d %.5f %.5f %10.5f\n", iter, x[0], x[1], f)
     end while status == GSL::CONTINUE and iter < 100

=== FMinimizer
     #!/usr/bin/env ruby
     require("gsl")
     include GSL::MultiMin

     np = 2
 
     my_f = Proc.new { |v, params|
       x = v[0]; y = v[1]
       p0 = params[0]; p1 = params[1]
       10.0*(x - p0)*(x - p0) + 20.0*(y - p1)*(y - p1) + 30.0
     }

     my_func = Function.alloc(my_f, np)
     my_func.set_params([1.0, 2.0])      # parameters

     x = Vector.alloc([5, 7])
     ss = Vector.alloc(np)
     ss.set_all(1.0)

     minimizer = FMinimizer.alloc("nmsimplex", np)
     minimizer.set(my_func, x, ss)

     iter = 0
     begin
       iter += 1
       status = minimizer.iterate()
       status = minimizer.test_size(1e-2)
       if status == GSL::SUCCESS
         puts("converged to minimum at")
       end
       x = minimizer.x
       printf("%5d ", iter);
       for i in 0...np do
         printf("%10.3e ", x[i])
       end
       printf("f() = %7.3f size = %.3f\n", minimizer.fval, minimizer.size);
     end while status == GSL::CONTINUE and iter < 100

((<prev|URL:multiroot.html>))
((<next|URL:fit.html>))

((<Reference index|URL:ref.html>))
((<top|URL:index.html>))

=end
