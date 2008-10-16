=begin
= Ordinary Differential Equations
This chapter describes functions for solving ordinary differential equation 
(ODE) initial value problems. The library provides a variety of low-level 
methods, such as Runge-Kutta and Bulirsch-Stoer routines, and higher-level 
components for adaptive step-size control. The components can be combined 
by the user to achieve the desired solution, with full access to any 
intermediate steps. 


Contents:
(1) ((<Classes for ODE solver|URL:odeiv.html#1>))
(2) ((<Class Descriptions|URL:odeiv.html#2>))
    (1) ((<GSL::Odeiv::System : Defining the ODE System|URL:odeiv.html#2.1>))
    (2) ((<GSL::Odeiv::Step : Stepping Algorithms|URL:odeiv.html#2.2>))
    (3) ((<GSL::Odeiv::Control : Adaptive Step-size Control|URL:odeiv.html#2.3>))
    (4) ((<GSL::Odeiv::Evolve : Evolution|URL:odeiv.html#2.4>))
    (5) ((<GSL::Odeiv::Solver : Higher level interface|URL:odeiv.html#2.5>))
(3) ((<Examples|URL:odeiv.html#3>))

== Classes for ODE solver

--- GSL::Odeiv::System
--- GSL::Odeiv::Step
--- GSL::Odeiv::Control
--- GSL::Odeiv::Evolve
    These are GSL structure wrappers.

--- GSL::Odeiv::Solver
    Another higher-level interface to ODE system classes.

== Class Descriptions

=== GSL::Odeiv::System
--- GSL::Odeiv::System.alloc(func, jac, dim)
--- GSL::Odeiv::System.alloc(func, dim)
    Constructor. This defines a general ODE system with the dimension ((|dim|)).
    
        # t: variable (scalar)
        # y: vector
        # dydt: vector
        # params: scalar or an array
 
        func = Proc.new { |t, y, dydt, params|
          mu = params
          dydt[0] = y[1]
          dydt[1] = -y[0] - mu*y[1]*(y[0]*y[0] - 1.0)
        }

        # t: scalar
        # y: vector
        # dfdy: matrix, jacobian
        # dfdt: vector
        # params: scalar of an array
        jac = Proc.new { |t, y, dfdy, dfdt, params|
          mu = params
          dfdy.set(0, 0, 0.0)
          dfdy.set(0, 1, 1.0)
          dfdy.set(1, 0, -2*mu*y[0]*y[1] - 1.0)
          dfdy.set(1, 1, -mu*(y[0]*y[0] - 1.0))
          dfdt[0] = 0.0
          dfdt[1] = 0.0
        }     

       sys = GSL:Odeiv::System.alloc(func, jac, dim)   # for "BSIMP" algorithm

    Note that some of the simpler solver algorithms do not make use of the 
    Jacobian matrix, so it is not always strictly necessary to provide it. 
    Thus the constructor is as
       sys = GSL:Odeiv::System.alloc(func, nil, dim)   # for others, replaced by nil
       sys = GSL:Odeiv::System.alloc(func, dim)        # or omit

--- GSL::Odeiv::System#set(func, jac, parameters...)
    This method sets or resets the procedures to evaluate the function and jacobian,
    and the constant parameters.

--- GSL::Odeiv::System#set_params(...)
    Set the constant parameters of the function.
--- GSL::Odeiv::System#function
--- GSL::Odeiv::System#func
--- GSL::Odeiv::System#jacobian
--- GSL::Odeiv::System#jac
    Return Proc objects

--- GSL::Odeiv::System#dimension
--- GSL::Odeiv::System#dim

=== GSL::Odeiv::Step
The lowest level components are the stepping functions which advance a solution from time ((|t|)) to ((|t+h|)) for a fixed step-size ((|h|)) and estimate the resulting local error.

--- GSL::Odeiv::Step.alloc(T, dim)
    Constructor for a stepping function of an algorithm type ((|T|)) for a system of 
    dimension ((|dim|)). The algorithms are specified by one of the constants under the 
    (({GSL::Odeiv::Step})) class, as

    (1) (({GSL::Odeiv::Step::RK2})), Embedded 2nd order Runge-Kutta with 3rd order error estimate.
    (2) (({GSL::Odeiv::Step::RK4})), 4th order (classical) Runge-Kutta.
    (3) (({GSL::Odeiv::Step::RKF45})), Embedded 4th order Runge-Kutta-Fehlberg method with 5th order error estimate. This method is a good general-purpose integrator.
    (4) (({GSL::Odeiv::Step::RKCK})), Embedded 4th order Runge-Kutta Cash-Karp method with 5th order error estimate.
    (5) (({GSL::Odeiv::Step::RK8PD})), Embedded 8th order Runge-Kutta Prince-Dormand method with 9th order error estimate.
    (6) (({GSL::Odeiv::Step::RK2IMP})), Implicit 2nd order Runge-Kutta at Gaussian points
    (7) (({GSL::Odeiv::Step::RK4IMP})), Implicit 4th order Runge-Kutta at Gaussian points
    (8) (({GSL::Odeiv::Step::BSIMP})), Implicit Bulirsch-Stoer method of Bader and Deuflhard. This algorithm requires the Jacobian.
    (9) (({GSL::Odeiv::Step::GEAR1})), M=1 implicit Gear method
    (10) (({GSL::Odeiv::Step::GEAR2})), M=2 implicit Gear method
    (11) (({GSL::Odeiv::Step::RK2SIMP})) (GSL-1.6)

    * Ex: 
        step = Odeiv::Step.alloc(Odeiv::Step::RKF45, 2)

    The algorithm types can also be given by a String, same as the C struct name,

    (1) "(({rk2}))" or "(({gsl_odeiv_step_rk2}))"
    (2) "(({rk4}))" or "(({gsl_odeiv_step_rk4}))"
    (3) "(({rkf45}))" or "(({gsl_odeiv_step_rkf45}))"
    (4) "(({rkck}))" or "(({gsl_odeiv_step_rkck}))"
    (5) "(({rk8pd}))" or "(({gsl_odeiv_step_rk8pd}))"
    (6) "(({rk2imp}))" or "(({gsl_odeiv_step_rk2imp}))"
    (7) "(({rk4imp}))" or "(({gsl_odeiv_step_rk4imp}))"
    (8) "(({bsimp}))" or "(({gsl_odeiv_step_bsimp}))"
    (9) "(({gear1}))" or "(({gsl_odeiv_step_gear1}))"
    (10) "(({gear2}))" or "(({gsl_odeiv_step_gear2}))"
    (10) "(({rk2simp}))" or "(({gsl_odeiv_step_rk2simp}))" (GSL-1.6)

    * Ex: 
        step = Odeiv::Step.alloc("bsimp", 4)
        step2 = Odeiv::Step.alloc("gsl_odeiv_step_rkck", 3)

--- GSL::Odeiv::Step#reset
    This method resets the stepper. It should be used whenever the next use 
    of s will not be a continuation of a previous step.

--- GSL::Odeiv::Step#name
    Returns the name of the stepper as a String. For example,

      require("gsl")
      include Odeiv
      s = Step.alloc(Step::RK4, 2)
      printf("step method is '%s'\n", s.name)

    would print something like step method is 'rk4'.

--- GSL::Odeiv::Step#order
    Returns the order of the stepper on the previous step. 
    This order can vary if the stepper itself is adaptive.

--- GSL::Odeiv::Step#apply(t, h, y, yerr, dydt_in, dydt_out, sys)
--- GSL::Odeiv::Step#apply(t, h, y, yerr, dydt_in, sys)
--- GSL::Odeiv::Step#apply(t, h, y, yerr, sys)
    This method applies the stepper to the system of equations defined by 
    ((|dydt|)), using the step size ((|h|)) to advance the system from time 
    ((|t|)) and state ((|y|)) to time ((|t+h|)). The new state of the system 
    is stored in ((|y|)) on output, with an estimate of the absolute error in 
    each component stored in ((|yerr|)). If the argument ((|dydt_in|)) is not 
    (({nil})) it should be a ((<GSL::Vector|URL:vector.html>)) object 
    containing the derivatives for the system at time ((|t|)) on input. 
    This is optional as the derivatives will be computed internally if they 
    are not provided, but allows the reuse of existing derivative information. 
    On output the new derivatives of the system at time ((|t+h|)) will be 
    stored in ((|dydt_out|)) if it is not nil.

=== GSL::Odeiv::Control
--- GSL::Odeiv::Control.standard_new(epsabs, epsrel, a_y, a_dydt)
--- GSL::Odeiv::Control.alloc(epsabs, epsrel, a_y, a_dydt)
    The standard control object is a four parameter heuristic based on 
    absolute and relative errors ((|epsabs|)) and ((|epsrel|)), and 
    scaling factors ((|a_y|)) and ((|a_dydt|)) for the system state 
    ((|y(t)|)) and derivatives ((|y'(t)|)) respectively.

--- GSL::Odeiv::Control.y_new(epsabs, epsrel)
    This method creates a new control object which will keep the local error 
    on each step within an absolute error of ((|epsabs|)) and relative error 
    of ((|epsrel|)) with respect to the solution ((|y_i(t)|)). 
    This is equivalent to the standard control object with ((|a_y=1|)) 
    and ((|a_dydt=0|)).

--- GSL::Odeiv::Control.yp_new(epsabs, epsrel)
    This method creates a new control object which will keep the local 
    error on each step within an absolute error of ((|epsabs|)) and 
    relative error of ((|epsrel|)) with respect to the derivatives of the 
    solution ((|y'_i(t)|)). 
    This is equivalent to the standard control object with ((|a_y=0|)) 
    and ((|a_dydt=1|)).

--- GSL::Odeiv::Control.alloc(epsabs, epsrel, a_y, a_dydt, vscale, dim)
    This creates a new control object which uses the same algorithm as 
    (({GSL::Odeiv::Control.standard_new})) but with an absolute error which 
    is scaled for each component by the (({GSL::Vector})) object ((|vscale|)).

--- GSL::Odeiv::Control#init(epsabs, epsrel, a_y, a_dydt)
    This method initializes the controler with the parameters ((|epsabs|)) 
    (absolute error),  ((|epsrel|)) (relative error), ((|a_y|)) 
    (scaling factor for y) and ((|a_dydt|)) (scaling factor for derivatives).

--- GSL::Odeiv::Control#name
--- GSL::Odeiv::Control#hadjust(step, y0, yerr, dydt, h)
    This method adjusts the step-size ((|h|)) using the control function 
    object, and the current values of ((|y|)),  ((|yerr|)) and ((|dydt|)). 
    The stepping function ((|step|)) is also needed to determine the order 
    of the method. On output, an array of two elements [((|hadj, status|))]
    is returned: If the error in the y-values  ((|yerr|)) is found to be 
    too large then the step-size ((|h|)) is reduced and the method returns 
    [((|hadj, status|))=(({GSL::ODEIV::HADJ_DEC}))]. 
    If the error is sufficiently small then 
    ((|h|)) may be increased and [((|hadj, status|))=(({GSL::ODEIV::HADJ_INC}))] 
    is returned. 
    The method returns [((|hadj, status|))=(({GSL::ODEIV::HADJ_NIL}))] if the step-size is 
    unchanged. The goal of the method is to estimate the largest step-size 
    which satisfies the user-specified accuracy requirements for the current 
    point.

=== GSL::Odeiv::Evolve
The higher level of the system is the (({GSL::Evolve})) class which combines the 
results of a stepper and controler to reliably advance the solution forward 
over an interval (({(t_0, t_1)})). If the controler signals that the step-size 
should be decreased the (({GSL::Evolve})) object backs out of the current step and 
tries the proposed smaller step-size. This process is continued until an 
acceptable step-size is found.

--- GSL::Odeiv::Evolve.alloc(dim)
    These create a (({GSL::Evolve})) object for a system of ((|dim|)) dimensions.

--- GSL::Odeiv::Evolve#reset
    This method resets the GSL::Evolve object. It should be used whenever 
    the next use of e will not be a continuation of a previous step.

--- GSL::Odeiv::Evolve#apply(evolve, control, step, sys, t, t1, h, y)
    This method advances the system ((|sys|)) from time ((|t|)) and position 
    ((|y|)) using the stepping function ((|step|)). The initial step-size is 
    taken as ((|h|)). The maximum time ((|t1|)) is guaranteed not to be exceeded by 
    the time-step. On output, an array of three elements is returned,
    [((|tnext, hnext, status|))], where ((|tnext|)) is the time advanced, 
    ((|hnext|)) is the step-size
    for the next step, and ((|status|)) is an error code retunred by (({gsl_odeiv_evolve_apply()})) function.
    On the final time-step the value of ((|tnext|)) will be set to 
    ((|t1|)) exactly.

--- GSL::Odeiv::Evolve#count

=== GSL::Odeiv::Solver
This is the highest level interface to solve ODE system, 
which contains System, Step, Control, and Evolve classes.

--- GSL::Odeiv::Solver.alloc(T, cary, fac, dim)
    This creates a ODE solver with the algorithm type ((|T|)) for the system of dimention ((|dim|)). Here ((|cary|)) is an array as an argument for the (({GSL::Odeive:Control})) constructor.

    * Ex1
        solver = Odeiv::Solver.alloc(Odeiv::Step::RKF45, [1e-6, 0.0], func, dim)

      * Type: RKF45,
      * Control: epsabs = 1e-6, epsrel = 0.0, a_y = 1, a_dydt = 0
      * System: function = ((|func|)), jacobian = ((|nil|))
      * Dimension: dim

    * Ex2:
        solver = Odeiv::Solver.alloc(Odeiv::Step::BSIMP, [1e-6, 0.0, 1, 0], func, jac, dim)

      * Type: BSIMP,
      * Control: epsabs = 1e-6, epsrel = 0.0, a_y = 1, a_dydt = 0
      * System: function = ((|func|)), jacobian = ((|jac|))
      * Dimension: dim

--- GSL::Odeiv:::Solver#reset
    Reset the solver elements (step, evolve)

--- GSL::Odeiv:::Solver#step
--- GSL::Odeiv:::Solver#control
--- GSL::Odeiv:::Solver#evolve
--- GSL::Odeiv:::Solver#system
    Access to the solver elements.

--- GSL::Odeiv::System#set_params(...)
    Set the constant parameters of the function to solve.

--- GSL::Odeiv:::Solver#apply(t, t1, h, y)
    This method advances the system from time ((|t|)) and position ((|y|)) ((({GSL::Vector})) object) using the stepping function. On output, the new time and position are returned as an array [((|tnext, hnext, status|))], i.e. ((|t, y|)) themselves are not modified by this method. The maximum time ((|t1|)) is guaranteed not to be exceeded by the time-step. On the final time-step the value of ((|tnext|)) will be set to ((|t1|)) exactly.

== Example

The following program solves the second-order nonlinear Van der Pol oscillator equation,
as found in the GSL manual, x"(t) + \mu x'(t) (x(t)^2 - 1) + x(t) = 0, 

This can be converted into a first order system suitable for use with the routines described in this chapter by introducing a separate variable for the velocity, y = x'(t),

  * x' = y
  * y' = -x + \mu y (1-x^2)

      require("gsl")
      include Odeiv

      dim = 2    # dimension of the system

      # Proc object to calculate the derivatives
      func = Proc.new { |t, y, dydt, mu|
        dydt[0] = y[1]
        dydt[1] = -y[0] - mu*y[1]*(y[0]*y[0] - 1.0)
      }

      # Create the solver
      solver = Solver.alloc(Step::RKF45, [1e-6, 0.0], func, dim)
      mu = 10.0
      solver.set_params(mu)

      t = 0.0       # initial time
      t1 = 100.0    # end time
      h = 1e-6      # initial step
      y = Vector.alloc([1.0, 0.0])    # initial value

      while t < t1
        t, h, status = solver.apply(t, t1, h, y)

        break if status != GSL::SUCCESS

        printf("%.5e %.5e %.5e %.5e\n", t, y[0], y[1], h)
      end


((<prev|URL:siman.html>))
((<next|URL:interp.html>))

((<Reference index|URL:ref.html>))
((<top|URL:index.html>))

=end
