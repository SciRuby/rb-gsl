=begin
= Monte Carlo Integration

== The GSL::Monte::Function class
The function to be integrated has its own datatype, the (({GSL::Monte::Function})) class.

--- GSL::Munte::Function.alloc(proc, dim, params)
--- GSL::Munte::Function.alloc(proc, dim)
    Constructor. The following example shows how to use this:

    * ex:
        proc_f = Proc.new { |x, dim, params|
          a = params[0]; b = params[1]; c = params[2]
          if dim != 2; raise("dim != 2"); end
          a*x[0]*x[0] + b*x[0]*x[1] + c*x[1]*x[1]
        }
        dim = 2
        mf = Monte::Function.alloc(proc_f, dim)
        mf.set_params([3, 2, 1])    

--- GSL::Munte::Function#set(proc, dim, params)
--- GSL::Munte::Function#set(proc, dim)
--- GSL::Munte::Function#set(proc)
--- GSL::Munte::Function#set_proc(proc)
--- GSL::Munte::Function#set_proc(proc, dim)
--- GSL::Munte::Function#set_params(params)
--- GSL::Munte::Function#params
--- GSL::Munte::Function#eval
--- GSL::Munte::Function#call

== Monte Carlo plans, alrgorithms 
=== PLAIN Monte Carlo
--- GSL::Monte::Plain.alloc(dim)
--- GSL::Monte::Plain#init
=== Miser
--- GSL::Monte::Miser.alloc(dim)
--- GSL::Monte::Miser#init
=== Vegas
--- GSL::Monte::Vegas.alloc(dim)
--- GSL::Monte::Vegas#init

== Integration
--- GSL:Monte::Function#integrate(xl, xu, dim, calls, rng, s)
--- GSL:Monte::Function#integrate(xl, xu, dim, calls, s)
--- GSL:Monte::Function#integrate(xl, xu, calls, rng, s)
--- GSL:Monte::Function#integrate(xl, xu, calls, s)
    This method performs Monte-Carlo integration of the function ((|self|)) 
    using the algorithm ((|s|)), over the ((|dim|))-dimensional hypercubic 
    region defined by the lower and upper 
    limits in the arrays ((|xl|)) and ((|xu|)), each of size ((|dim|)). 
    The integration uses a fixed number of function calls ((|calls|)).
    The argument ((|rng|)) is a random number generator (optional). If it is not
    given, a new generator is created internally and freed when the calculation
    finishes.

    See sample scripts (({sample/monte*.rb})) for more details.

== Accessing internal state of the Monte Carlo classes
--- GSL::Monte::Miser#estimate_frac
--- GSL::Monte::Miser#estimate_frac=
--- GSL::Monte::Miser#min_calls
--- GSL::Monte::Miser#min_calls=
--- GSL::Monte::Miser#min_call_per_bisection
--- GSL::Monte::Miser#min_calls_per_bisection=
--- GSL::Monte::Miser#alpha
--- GSL::Monte::Miser#alpha=
--- GSL::Monte::Miser#dither
--- GSL::Monte::Miser#dither=
--- GSL::Monte::Vegas#alpha
--- GSL::Monte::Vegas#result
--- GSL::Monte::Vegas#sigma
--- GSL::Monte::Vegas#chisq
--- GSL::Monte::Vegas#iterations
--- GSL::Monte::Vegas#iterations=
--- GSL::Monte::Vegas#alpha
--- GSL::Monte::Vegas#alpha=
--- GSL::Monte::Vegas#stage
--- GSL::Monte::Vegas#stage=
--- GSL::Monte::Vegas#mode
--- GSL::Monte::Vegas#mode=
--- GSL::Monte::Vegas#verbose
--- GSL::Monte::Vegas#verbose=

== Example

     #!/usr/bin/env ruby
     require("gsl")
     include GSL::Monte
     include Math

     proc_f = Proc.new { |k, dim, params|
       pi = Math::PI
       a = 1.0/(pi*pi*pi)
       a/(1.0 - cos(k[0])*cos(k[1])*cos(k[2]))
     }

     def display_results(title, result, error)
       exact = 1.3932039296856768591842462603255

       diff = result - exact
       printf("%s ==================\n", title);
       printf("result = % .6f\n", result);
       printf("sigma  = % .6f\n", error);
       printf("exact  = % .6f\n", exact);
       printf("error  = % .6f = %.1g sigma\n", diff, diff.abs/error)
     end

     dim = 3
     xl = Vector.alloc(0, 0, 0)
     xu = Vector.alloc(PI, PI, PI)
     G = Monte::Function.alloc(proc_f, dim)
     calls = 500000
     r = GSL::Rng.alloc(Rng::DEFAULT)

     plain = Monte::Plain.alloc(dim)
     result, error = G.integrate(xl, xu, dim, calls, r, plain)
     display_results("plain", result, error)

     miser = Monte::Miser.alloc(dim)
     result, error = G.integrate(xl, xu, dim, calls, r, miser)
     display_results("miser", result, error)

     vegas = Monte::Vegas.alloc(dim)
     result, error = G.integrate(xl, xu, dim, 10000, r, vegas)
     display_results("vegas warm-up", result, error)
     puts("converging...");
     begin
       result, error = G.integrate(xl, xu, dim, calls/5, r, vegas)
       printf("result = % .6f sigma = % .6f chisq/dof = %.1f\n", 
               result, error, vegas.chisq)
     end while (vegas.chisq-1.0).abs > 0.5
     display_results("vegas final", result, error)

((<prev|URL:ntuple.html>))
((<next|URL:siman.html>))

((<Reference index|URL:ref.html>))
((<top|URL:index.html>))

=end
