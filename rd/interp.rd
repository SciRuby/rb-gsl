=begin
= Interpolation
This chapter describes functions for performing interpolation. 
The library provides a variety of interpolation methods, including 
Cubic splines and Akima splines. The interpolation types are interchangeable, 
allowing different methods to be used without recompiling. Interpolations can 
be defined for both normal and periodic boundary conditions. Additional 
functions are available for computing derivatives and integrals of 
interpolating functions. 

(1) ((<Interpolation classes|URL:interp.html#1>))
(2) ((<Initializing interpolation objects|URL:interp.html#2>))
(3) ((<Index Look-up and Acceleration|URL:interp.html#3>))
(4) ((<Evaluation of Interpolating Functions|URL:interp.html#4>))
(5) ((<Higher level interface: GSL::Spline class|URL:interp.html#5>))
    (1) ((<Class initialization|URL:interp.html#5.1>))
    (2) ((<Evaluation|URL:interp.html#5.2>))
    (3) ((<Finding and acceleration|URL:interp.html#5.3>))

== Interpolation Classes
* GSL
  * Interp (class)
    * Accel (class)
  * Spline (class)
  
== Initializing interpolation objects

--- GSL::Interp.alloc(T, n)
--- GSL::Interp.alloc(T, x, y)
--- GSL::Interp.alloc(x, y)
    These methods create an  interpolation object of type ((|T|)) for ((|n|)) 
    data-points.

    The library provides six types, which are specifiled by an identifier of a 
    constant or a string:

    * Interp::LINEAR or "linear"

      Linear interpolation. This interpolation method does not require any additional memory.   
    * Interp::POLYNOMIAL or "polynomial"

      Polynomial interpolation. This method should only be used for interpolating small numbers of points because polynomial interpolation introduces large oscillations, even for well-behaved datasets. The number of terms in the interpolating polynomial is equal to the number of points.

    * Interp::CSPLINE or "cspline"

      Cubic spline with natural boundary conditions.
    * Interp::CSPLINE_PERIODIC or "gsl_cspline_periodic" or "cspline_periodic"

      Cubic spline with periodic boundary conditions
    * Interp::AKIMA or "akima"    

      Non-rounded Akima spline with natural boundary conditions. This method uses the non-rounded corner algorithm of Wodicka.
    * Interp::AKIMA_PERIODIC or "akima_periodic"    

      Non-rounded Akima spline with periodic boundary conditions. This method uses the non-rounded corner algorithm of Wodicka.
    
      * ex: For cubic spline for 10 points,
          sp = Interp.alloc("cspline", 10)

--- GSL::Interp#init(xa, ya)
    This method initializes the interpolation object interp for the data 
    ((|(xa,ya)|)) where ((|xa|)) and ((|ya|)) are vectors. 
    The interpolation object ((({GSL::Interp}))) does not save the data 
    vectors ((|xa, ya|)) and only stores the static state computed from the data. 
    The ((|xa|)) vector is always assumed to be strictly ordered; the behavior 
    for other arrangements is not defined.


--- GSL::Interp#name
    This returns the name of the interpolation type used by ((|self|)).



--- GSL::Interp#min_size
    This returns the minimum number of points required by the interpolation 
    type of ((|self|)). For example, Akima spline interpolation requires 
    a minimum of 5 points.

== Index Look-up and Acceleration 
--- GSL::Interp.bsearch(xa, x, index_lo, index_hi)
    This returns the index i of the vector ((|xa|)) such that 
    (({xa[i] <= x < x[i+1]})). The index is searched for in the range 
    ((|[index_lo,index_hi]|)).


--- GSL::Interp#accel
    In C level, the library requires a (({gsl_interp_accel})) object, 
    but it is hidden in Ruby/GSL. It is automatically allocated 
    when a (({GSL::Interp})) object is created, stored in it, 
    and destroyed when the (({Interp})) object 
    is cleaned by the Ruby GC. 
    This method is used to access to the (({Interp::Accel})) object
    stored in ((|self|)).

--- GSL::Interp#find(xa, x)
--- GSL::Interp#accel_find(xa, x)
--- GSL::Interp::Accel#find(xa, x)
    This method performs a lookup action on the data array ((|xa|)). 
    This is how lookups are performed during evaluation 
    of an interpolation. The function returns an index (({i})) such that 
    (({xa[i] <= x < xa[i+1]})).


== Evaluation of Interpolating Functions 

--- GSL::Interp#eval(xa, ya, x)
--- GSL::Interp#eval_e(xa, ya, x)
    These methods return the interpolated value for a given point ((|x|)), 
    using the interpolation object ((|self|)), data vectors ((|xa|)) and ((|ya|)).
    The data ((|x|)) can be a (({Numeric, Vector, Matrix})) or an (({NArray})).
--- GSL::Interp#eval_deriv(xa, ya, x)
--- GSL::Interp#eval_deriv_e(xa, ya, x)
    These methods return the derivative of an interpolated function for a 
    given point ((|x|)), using the interpolation object ((|self|)), 
    data vectors ((|xa|)) and ((|ya|)).

--- GSL::Interp#eval_deriv2(xa, ya, x)
--- GSL::Interp#eval_deriv2_e(xa, ya, x)
    These methods return the second derivative of an interpolated function 
    for a given point ((|x|)), using the interpolation object ((|self|)), 
    data vectors ((|xa|)) and ((|ya|)).

--- GSL::Interp#eval_integ(xa, ya, a, b)
--- GSL::Interp#eval_integ_e(xa, ya, a, b)
    These methods return the numerical integral result of an interpolated 
    function over the range ((|[a, b]|)), using the interpolation object ((|self|)), 
    data vectors ((|xa|)) and ((|ya|)).

== Higher level interface: GSL::Spline class
=== Class initialization

--- GSL::Spline.alloc(T, n)
--- GSL::Spline.alloc(T, x, y)
--- GSL::Spline.alloc(x, y, T)
    This creates a (({GSL::Spline})) object of type ((|T|)) for ((|n|)) 
    data-points. The type ((|T|)) is the same as (({GSL::Interp})) class.

    These two are equivalent.
    * (({GSL::Spline.alloc})) and (({GSL::Spline#init}))
        sp = GSL::Spline.alloc(T, n)
        sp.init(x, y)                 # x and y are vectors of length n
    * (({GSL::Spline.alloc})) with two vectors
        sp = GSL::Spline.alloc(T, x, y)
    If ((|T|)) is not given, "cspline" is used.

--- GSL::Spline#init(xa, ya)
    This initializes a (({GSL::Spline})) object ((|self|)) for the data 
    (((|xa, ya|))) where ((|xa|)) and ((|ya|)) are Ruby arrays of equal sizes 
    or (({GSL::Vector})).

--- GSL::Spline#name
    This returns the name of the spline type used by ((|self|)).

=== Evaluation
--- GSL::Spline#eval(x)
    This returns the interpolated value for a given point ((|x|)).
    The data ((|x|)) can be a (({Numeric, Vector, Matrix})) or an (({NArray})).

    NOTE: In a GSL-C program, a (({gsl_interp_accel})) object is required to use
    the function (({gsl_spline_eval})).
    In Ruby/GSL, the (({gsl_interp_accel})) is hidden, it is automatically 
    allocated when a (({GSL::Spline})) object is created, 
    and also destroyed when the (({Spline})) object 
    is cleaned by the Ruby GC. The accel object can be accessed via the method 
    (({GSL::Spline#accel})).

--- GSL::Spline#eval_deriv(x)
    This returns the derivative of an interpolated function for a given point ((|x|)), usingthe data arrays ((|xa|)) and ((|ya|)) set by ((|init|)).

--- GSL::Spline#eval_deriv2(x)
    This returns the second derivative at ((|x|)).

--- GSL::Spline#eval_integ(a, b)
    Returns the numerical integral over the range [((|a, b|))].

=== Finding and acceleration
--- GSL::Spline#find(xa, x)
--- GSL::Spline#accel_find(xa, x)
    This method performs a lookup action on the data array ((|xa|)). 
    This is how lookups are performed during evaluation 
    of an interpolation. The function returns an index (({i})) such that 
    (({xa[i] <= x < xa[i+1]})).

See also the GSL manual and the examples in (({examples/}))

((<prev|URL:odeiv.html>))
((<next|URL:diff.html>))

((<Reference index|URL:ref.html>))
((<top|URL:index.html>))

=end
  
