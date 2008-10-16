=begin
= Numerical Integration
Contents:
(1) ((<Introduction|URL:integration.html#1>))
(2) ((<QNG non-adaptive Gauss-Kronrod integration|URL:integration.html#2>))
(3) ((<QAG adaptive integration|URL:integration.html#3>))
    (1) ((<GSL::Integration::Workspace class|URL:integration.html#3.1>))
    (2) ((<Methods|URL:integration.html#3.2>))
(4) ((<QAGS adaptive integration with singularities|URL:integration.html#4>))
(5) ((<QAGP adaptive integration with known singular points|URL:integration.html#5>))
(6) ((<QAGI adaptive integration on infinite intervals|URL:integration.html#6>))
(7) ((<QAWC adaptive integration for Cauchy principal values|URL:integration.html#7>))
(8) ((<QAWS adaptive integration for singular functions|URL:integration.html#8>))
(9) ((<QAWO adaptive integration for oscillatory functions|URL:integration.html#9>))
(10) ((<QAWF adaptive integration for Fourier integrals|URL:integration.html#10>))

== Introduction
This section describes how to compute numerical integration of a function 
in one dimension. In Ruby/GSL, all the GSL routines for numerical integration 
is provided as methods of ((<GSL::Function|URL:function.html>)) objects.
For example, a (({GSL::Function})) object which represents the sine function 
((|sin(x)|)) can be expressed as
  f = GSL::Function.alloc { |x| sin(x) }
To compute numerical integration of ((|sin(x)|)) over the range ((|(a, b)|)), 
one can use the methods ((|integrate_xxx|)) or simply ((|xxx|)), as
  * (({f.integrate_xxx([a, b])})), or (({f.xxx([a, b])}))
  * (({f.integrate_xxx(a, b)})), or (({f.xxx(a, b)}))

== QNG non-adaptive Gauss-Kronrod integration
--- GSL::Function#integration_qng([a, b], [epsabs = 0.0, epsrel = 1e-10])
--- GSL::Function#qng(...)
--- GSL::Integration::qng(...)
    These methods apply the Gauss-Kronrod integration rules in succession until 
    an estimate of the integral of the reciever function (a (({GSL::Function})) 
    object) over ((|(a,b)|)) is achieved within the desired absolute and relative 
    error limits, ((|epsabs|)) and ((|epsrel|)) (these are optional, the default
    values are 0,0 and 1e-10 respectively). These methods return an array of
    four elements (({[result, err, neval, status]})), those are the final 
    approximation 
    of the integration, an estimate of the absolute error, the number of function 
    evaluation, and the status which is returned by the GSL 
    (({integration_qng()})) function.

    * Ex: Integrate sin(x) over ((|x = 0 -- 2|)) with accuracies ((|epsabs = 0, epsrel = 1.0e-7|)).

        require 'gsl'

        f = GSL::Function.alloc { |x| sin(x) }
        ans = f.integration_qng([0, 2], [0, 1.0e-7])   # or shortly f.qng(...)
        p ans[0]   <- result

    For all the methods described in this section, the arguments ((|[epsabs, epsrel]|)) are optional, and the default values are ((|[epsabs = 0.0, epsrel = 1e-10]|)).

== QAG adaptive integration
The QAG algorithm is a simple adaptive integration procedure. 
The integration region is divided into subintervals, and on each iteration 
the subinterval with the largest estimated error is bisected. 
This reduces the overall error rapidly, as the subintervals become concentrated
around local difficulties in the integrand. These subintervals are managed by 
a GSL::Integration::Workspace object, which handles the memory for the 
subinterval ranges, results and error estimates.

=== GSL::Integration::Workspace class
--- GSL::Integration::Workspace.alloc(n = 1000)
    This creates a workspace sufficient to hold n double precision intervals, 
    their integration results and error estimates.

--- GSL::Integration::Workspace#limit
--- GSL::Integration::Workspace#size

==== Algorithms which require the workspace
The algorithms described below require (({gsl_integration_workspace})) struct
in C. In Ruby/GSL, the corresponding methods require 
a (({GSL::Integration::Workspace})) object in their arguments. But it is also
possible to use these methods without workspace arguments: if it 
is not given, a workspace is created/destroyed internally. Thus
method calls are as

          f = GSL::Function.alloc { |x| Math::sin(x)/x }
          p f.qag([a, b])

or
          w = GSL::Integration::Workspace.alloc(limit)
          p f.qag([a, b], w)

Explicit uses of a (({Workspace})) object reduce C function calls for memory 
allocations of workspace objects. 

=== Methods
--- GSL::Function#integration_qag([a, b], key = GSL::Integration::GAUSS61)
--- GSL::Function#integration_qag([a, b], key, w)
--- GSL::Function#integration_qag([a, b], w)
--- GSL::Function#integration_qag([a, b], [epsabs, epsrel], key)
--- GSL::Function#integration_qag([a, b], [epsabs, epsrel], key, w)
--- GSL::Function#qag(...)
--- GSL::Integration::qag(...)
    These methods apply an integration rule adaptively until an estimate of the 
    integral of the reciever function over ((|(a,b)|)) is achieved within the 
    desired absolute and relative error limits, ((|epsabs|)) and ((|epsrel|)).
    One can give a (({GSL::Integration::Workspace})) object ((|w|)) with the
    last argument (option: if not given, the workspace is internally allocated and
    freed). The method returns an array with four elements 
    (({[result, err, neval, status]})). 
    The integration rule is determined by the value of key, which should be 
    chosen from the following symbolic names, 

       GSL::Integration::GAUSS15  (key = 1)
       GSL::Integration::GAUSS21  (key = 2)
       GSL::Integration::GAUSS31  (key = 3)
       GSL::Integration::GAUSS41  (key = 4)
       GSL::Integration::GAUSS51  (key = 5)
       GSL::Integration::GAUSS61  (key = 6)

    corresponding to the 15, 21, 31, 41, 51 and 61 point Gauss-Kronrod rules. The 
    higher-order rules give better accuracy for smooth functions, 
    while lower-order rules save time when the function contains local 
    difficulties, such as discontinuities. 

== QAGS adaptive integration with singularities
The presence of an integrable singularity in the integration region causes 
an adaptive routine to concentrate new subintervals around the singularity. 
As the subintervals decrease in size the successive approximations to the 
integral converge in a limiting fashion. This approach to the limit can be 
accelerated using an extrapolation procedure. 
The QAGS algorithm combines adaptive bisection with the Wynn epsilon-algorithm 
to speed up the integration of many types of integrable singularities.

--- GSL::Function#integration_qags([a, b], [epsabs = 0.0, epsrel = 1e-10], limit)
--- GSL::Function#integration_qags([a, b], [epsabs, epsrel], limit, w)
--- GSL::Function#integration_qags([a, b], [epsabs, epsrel], w)
--- GSL::Function#qags(...)
--- GSL::Integration::qags(...)
    These methods apply the Gauss-Kronrod 21-point integration rule 
    adaptively until an estimate of the integral over ((|(a,b)|)) is 
    achieved within the desired absolute and relative error limits, 
    ((|epsabs|)) and ((|epsrel|)). The results are extrapolated using the 
    epsilon-algorithm, which accelerates the convergence of the integral 
    in the presence of discontinuities and integrable singularities. 
    The maximum number of subintervals is given by ((|limit|)).

    * ex:

        proc = Proc.new{ |x, alpha|     # integrant
          log(alpha*x)/sqrt(x)
        }

        # create the function, with the parameter alpha = 1.0
        f = GSL::Function.alloc(proc, 1.0) 

        p f.integration_qags(0, 1)

== QAGP adaptive integration with known singular points
--- GSL::Function#integration_qagp(pts, [epsabs = 0.0, epsrel = 1e-10], limit = 1000, w)
--- GSL::Function#qagp(...)
--- GSL::Integration::qagp(...)
    These methods apply the adaptive integration algorithm QAGS taking 
    account of the user-supplied locations of singular points. The array 
    ((|pts|)) (a Ruby array or a GSL::Vector object) should contain the 
    endpoints of the integration ranges defined by the integration region a
    nd locations of the singularities. For example, to integrate over the 
    region ((|(a,b)|)) with break-points at x_1, x_2, x_3 
    (where a < x_1 < x_2 < x_3 < b) the following ((|pts|)) array should be used

          pts[0] = a
          pts[1] = x_1
          pts[2] = x_2
          pts[3] = x_3
          pts[4] = b

    If you know the locations of the singular points in the integration region then this routine will be faster than QAGS.

    * ex:
       f454 = Function.alloc{ |x|
         x2 = x*x
         x3 = x2*x
         x3*log(((x2-1)*(x2-2)).abs)
       }
       pts = [0, 1, sqrt(2), 3]     # range: [0, 3], singular points: [1, sqrt(2)]
       p f454.qagp(pts, 0.0, 1e-3)  # <---- [52.7408061167272, 0.000175570384826074, 20, 0]
                                    # Expect: 61 log(2) + (77/4) log(7) - 27 = 52.7408061167272

== QAGI adaptive integration on infinite intervals
--- GSL::Function#integration_qagi([epsabs = 0.0, epsrel = 1e-10], limit = 1000, w)
--- GSL::Function#qagi(...)
--- GSL::Integration::qagi(...)
    These methods compute the integral of the function over the infinite 
    interval (-infty,+infty). The integral is mapped onto the interval
    (0,1] using the transformation x = (1-t)/t. It is then integrated using 
    the QAGS algorithm. The normal 21-point Gauss-Kronrod rule of QAGS is 
    replaced by a 15-point rule, because the transformation can generate an 
    integrable singularity at the origin. In this case a lower-order rule is 
    more efficient.

    * ex
        f = Function.alloc{ |x| Math::exp(-x*x) }
        exact = Math::sqrt(Math::PI)

        result, = f.qagi
        puts("exp(-x*x), x = -infty --- +infty")
        printf("exact  = %.18f\n", exact)
        printf("result = %.18f\n\n", result)

--- GSL::Function#integration_qagiu(a, epsabs = 0.0, epsrel = 1e-10, limit = 1000)
--- GSL::Function#integration_qagiu(a, epsabs = 0.0, epsrel = 1e-10, w)
--- GSL::Function#qagiu(...)
--- GSL::Integration::qagiu(...)
    These methods compute the integral of the function over the 
    semi-infinite interval (a,+infty).
    
--- GSL::Function#integration_qagil(b, epsabs = 0.0, epsrel = 1e-10, limit = 1000)
--- GSL::Function#integration_qagil(b, epsabs = 0.0, epsrel = 1e-10, w)
--- GSL::Function#integration_qagil(b, [epsabs, epsrel], limit, w)
--- GSL::Function#qagil(...)
--- GSL::Integration::qagil(...)
    These methods compute the integral of the function over the 
    semi-infinite interval (-infty,b).

== QAWC adaptive integration for Cauchy principal values
--- GSL::Function#integration_qawc([a, b], c, [epsabs = 0.0, epsrel = 1e-10], limit. 1000)
--- GSL::Function#qawc(...)
--- GSL::Function#qawc(...)
    These methods compute the Cauchy principal value of the integral over 
    ((|(a,b)|)), with a singularity at ((|c|)). The adaptive bisection algorithm 
    of QAG is used, with modifications to ensure that subdivisions do not occur 
    at the singular point ((|x = c|)). When a subinterval contains the point 
    ((|x = c|)) or is close to it then a special 25-point modified 
    Clenshaw-Curtis rule is used to control the singularity. Further away from 
    the singularity the algorithm uses an ordinary 15-point Gauss-Kronrod 
    integration rule.

    * ex:
         require 'gsl'
         f459 = Function.alloc { |x| 1.0/(5.0*x*x*x + 6.0) }

         p f459.qawc([-1.0, 5.0], 0, [0.0, 1e-3]) # Expect: log(125/631)/18

== QAWS adaptive integration for singular functions
The QAWS algorithm is designed for integrands with algebraic-logarithmic 
singularities at the end-points of an integration region. 
In order to work efficiently the algorithm requires a precomputed 
table of Chebyshev moments.

--- GSL::Function#integration_qaws([a, b], table, [epsabs = 0.0, epsrel = 1e-10], limit = 1000)
--- GSL::Function#integration_qaws(a, b, table, epsabs, epsrel, limit, w)
--- GSL::Function#qaws(...)
--- GSL::Integration::qaws(...)
    These methods compute the integral of the function over the interval 
    ((|(a,b)|)) with the singular weight function 
       (x-a)^alpha (b-x)^beta log^mu (x-a) log^nu (b-x)
    The parameters (({[alpha, beta, mu, nu]})) is given by a Ruby array
    ((|table|)), or by a (({GSL::Integration::QAWS_Table})) object.

    The adaptive bisection algorithm of QAG is used. When a subinterval contains one of the endpoints then a special 25-point modified Clenshaw-Curtis rule is used to control the singularities. For subintervals which do not include the endpoints an ordinary 15-point Gauss-Kronrod integration rule is used.

    * ex1:
           f458 = Function.alloc { |x|
             if x.zero?
               val = 0.0
             else
               u = log(x)
               v = 1.0 + u*u
               val = 1.0/(v*v)
             end
             val
           }
           table = [0.0, 0.0, 1, 0]
           p f458.qaws([0.0, 1.0], table, [0.0, 1e-10])  # Expect: -0.1892752

    * ex2:
         table = Integration::QAWS_Table.alloc(0.0, 0.0, 1, 0)
         p f458.qaws([0.0, 1.0], table, [0.0, 1e-10])

== QAWO adaptive integration for oscillatory functions
The QAWO algorithm is designed for integrands with an oscillatory factor, 
sin(omega x) or cos(omega x). In order to work efficiently the algorithm 
requires a table of Chebyshev moments.


--- GSL::Function#integration_qawo(a, [epsabs = 0.0, epsrel = 1e-10], limit = 1000, w, table, )
--- GSL::Function#qawo(...)
--- GSL::Integration::qawo(...)
    This method uses an adaptive algorithm to compute the integral over 
    ((|[a,b]|)) with the weight function sin(omega x) or cos(omega x) defined by
    the table ((|table|)).

    * ex1:
          require("gsl")
          f456 = Function.alloc { |x|
            if x.zero?
              val = 0.0
            else
              val = Math::log(x)
            end
            val
          }
          table = [10.0*Math::PI, 1.0, Integration::SINE, 1000]
          p f456.qawo(0.0, [0.0, 1e-10], table)

    * ex2:
          table = Integration::QAWO_Table.alloc(10.0*Math::PI, 1.0, Integration::SINE, 1000)
          p f456.qawo(0.0, [0.0, 1e-10], table)

== QAWF adaptive integration for Fourier integrals
--- GSL::Function#integration_qawf(a, epsabs = 1e-10, limit = 1000, w, wc, table)
--- GSL::Function#integration_qawf(a, epsabs = 1e-10, limit = 1000, table)
--- GSL::Function#integration_qawf(a, epsabs = 1e-10, table)
--- GSL::Function#integration_qawf(a, table = 1000, table)
--- GSL::Function#integration_qawf(a, table)
--- GSL::Function#qawf(...)
--- GSL::Integration::qawf(...)
    This method attempts to compute a Fourier integral of the function over 
    the semi-infinite interval [a,+infty).

        I = \int_a^{+infty} dx f(x) sin(omega x)
        I = \int_a^{+infty} dx f(x) cos(omega x)

    The parameter omega is taken from the table ((|table|)) (the length ((|L||)) 
    can take any value, since it is overridden by this function to a value 
    appropriate for the fourier integration).

    * ex:
         f457 = Function.alloc { |x|
           if x.zero?
             val = 0.0
           else
             val = 1.0/Math::sqrt(x)
           end
           val
         }
         table = [PI/2.0, 1.0, GSL::Integration::COSINE, 1000]
         p f457.qawf(0.0, 1e-10, table)     #  0.999999999927979, Expect 1

    In other style:

         limit = 1000
         table = Integration::QAWO_Table.alloc(PI/2.0, 1.0, GSL::Integration::COSINE, 1000)
         w = Integration::Workspace.alloc          # default n is 1000
         wc = Integration::Workspace.alloc(limit)

         p f457.qawf(0.0, table)
         p f457.qawf(0.0, 1e-10, table)
         p f457.qawf(0.0, 1e-10, limit, table)
         p f457.qawf(0.0, limit, table)
         p f457.qawf(0.0, 1e-10, limit, w, wc, table)
         p f457.qawf(0.0, w, wc, table)
         p f457.qawf(0.0, limit, w, wc, table)
         p f457.qawf(0.0, limit, w, table)       # Error
         p f457.qawf(0.0, limit, wc, table)      # Error

((<prev|URL:wavelet.html>))
((<next|URL:rng.html>))

((<Reference index|URL:ref.html>))
((<top|URL:index.html>))

=end
