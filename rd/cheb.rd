=begin
= Chebyshev Approximations
This chapter describes routines for computing Chebyshev approximations to 
univariate functions. A Chebyshev approximation is a truncation of the series 
f(x) = \sum c_n T_n(x), 
where the Chebyshev polynomials T_n(x) = \cos(n \arccos x) 
provide an orthogonal basis of polynomials on the interval [-1,1] 
with the weight function 1 / \sqrt{1-x^2}. 
The first few Chebyshev polynomials are, 
T_0(x) = 1, T_1(x) = x, T_2(x) = 2 x^2 - 1. 
For further information see Abramowitz & Stegun, Chapter 22. 

(1) ((<GSL::Cheb class|URL:cheb.html#1>))
(2) ((<Chebyshev Series Evaluation|URL:cheb.html#2>))
(3) ((<Derivatives and Integrals|URL:cheb.html#3>))
(4) ((<Examples|URL:cheb.html#4>))

== (({GSL::Cheb})) class

--- GSL::Cheb.alloc(n) 
    This create an instance of the GSL::Cheb class for a Chebyshev series of order n.


--- GSL::Cheb#init(f, a, b)
    This computes the Chebyshev approximation the function ((|f|)) over the range (((|a,b|))) to the previously specified order. Where ((|f|)) is a ((<GSL::Function|URL:function.html>)) object. The computation of the Chebyshev approximation is an O(n^2) process, and requires ((|n|)) function evaluations.

    * ex: Approximate a step function defined in (0, 1) by a Chebyshev series of order 40.
        f = GSL::Function.alloc { |x|
          if x < 0.5
            0.25
          else
            0.75
          end
        }

        cs = GSL::Cheb.alloc(40)
        cs.init(f, 0, 1)

== Chebyshev Series Evaluation
--- GSL::Cheb#eval(x)
    This evaluates the Chebyshev series at a given point ((|x|)).

--- GSL::Cheb#eval_n(n, x)
    This evaluates the Chebyshev series at a given point ((|x|)), to (at most) the given order ((|n|)).

== Derivatives and Integrals 

--- GSL::Cheb#calc_deriv()
--- GSL::Cheb#deriv()
    This computes the derivative of the series, and returns a new GSL::Cheb object which contains the computed derivative. The reciever is not changed.

--- GSL::Cheb#calc_integ()
--- GSL::Cheb#integ()
    This computes the integral of the series, and returns a new GSL::Cheb object which contains the computed integral coefficients. The reciever is not changed.

== Example
   #!/usr/bin/env ruby
   require("gsl")

   f = GSL::Function.alloc { |x|
     if x < 0.5
       0.25
     else
       0.75
     end
   }

   n = 1000
   order = 40
   cs = GSL::Cheb.alloc(order)
   cs.init(f, 0, 1)

   x = Vector.linspace(0, 1, n)
   ff = f.eval(x)
   r10 = cs.eval_n(10, x)
   r40 = cs.eval(x)
   GSL::graph(x, ff, r10, r40)


See also the example scripts in (({examples/cheb/})).

((<prev|URL:diff.html>))
((<next|URL:sum.html>))

((<Reference index|URL:ref.html>))
((<top|URL:index.html>))
=end
