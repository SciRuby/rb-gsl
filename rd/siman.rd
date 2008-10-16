=begin
= Simulated Annealing
=== Library

== Module and classes
  * GSL::
    * Siman::  (Module)
      * Params (Class)
      * Efunc (Class)
      * Step (Class)
      * Metric (Class)
      * Print (Class)

== (({Siman})) Module
=== Singleton method
--- GSL::Siman.solve(rng, x0_p, efunc, stepper, metric, printer, params)
    This performs a simulated annealing search through a given space. 
    The space is specified by providing the functions ((|efunc|)) and ((|metric|)). 
    The simulated annealing steps are generated using the random number generator 
    ((|rng|)) and the function ((|stepper|)). The starting configuration of the 
    system should be given by a (({Vector})) object ((|x0_p|)).
  
    The parameter ((|params|)) controls the run by providing the temperature 
    schedule and other tunable parameters to the algorithm.

    On exit the best result achieved during the search is placed in ((|x0_p|)). 
    If the annealing process has been successful this should be a good approximation 
    to the optimal point in the space.

    If the function ((|printer|)) is not (({nil})), a debugging log will be printed 
    to stdout with the following columns:
        number_of_iterations   temperature    x    x-(x0_p)   efunc(x)
    and the output of ((|printer|)) itself. If ((|printer|)) is (({nil})) 
    then no information is printed.

== Example

     #!/usr/bin/env ruby
     require("gsl")
     include GSL::Siman

     N_TRIES = 200
     ITERS_FIXED_T = 10
     STEP_SIZE = 10
     K = 1.0
     T_INITIAL = 0.002
     MU_T = 1.005
     T_MIN = 2.0e-6

     params = Siman::Params.alloc(N_TRIES, ITERS_FIXED_T, STEP_SIZE, K, T_INITIAL,
                                MU_T, T_MIN)

     efunc = Efunc.alloc { |vx|
       x = vx[0]
       s = (x - 1.0)
       Math::exp(-s*s)*Math::sin(8*x)
     }

     metric = Metric.alloc { |vx, vy|
       (x[0] - y[0]).abs
     }

     step = Step.alloc { |rng, vx, step_size|
       r = rng.uniform
       old_x = vx[0]
       a =  r * 2 * step_size - step_size + old_x
       vx[0] = a
     }

     simanprint = Print.alloc { |vx|
       printf("%12g", vx[0])
     }

     x = Vector.alloc([15.5])
     Rng.env_setup()
     rng = Rng.alloc()

     #Siman::solve(rng, x, efunc, step, metric, simanprint, params)
     Siman::solve(rng, x, efunc, step, metric, nil, params)
     p x

((<prev|URL:monte.html>))
((<next|URL:odeiv.html>))

((<Reference index|URL:ref.html>))
((<top|URL:index.html>))
=end
