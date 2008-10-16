=begin

= Polynomials
Contents:
(1) ((<Polynomial Evaluation|URL:poly.html#1>))
(2) ((<Solving polynomial equations|URL:poly.html#2>))
    (1) ((<Quadratic Equations|URL:poly.html#2.1>))
    (2) ((<Cubic Equations|URL:poly.html#2.2>))
    (3) ((<General Polynomial Equations|URL:poly.html#2.3>))
(3) ((<GSL::Poly Class|URL:poly.html#3>))
    (1) ((<Constructors|URL:poly.html#3.1>))
    (2) ((<Methods|URL:poly.html#3.2>))
(4) ((<Polynomial Fitting|URL:poly.html#4>))
(5) ((<Divided-difference representations|URL:poly.html#5>))
(6) ((<Extensions|URL:poly.html#6>))
    (1) ((<Special Polynomials|URL:poly.html#6.1>))
    (2) ((<Polynomial Operations|URL:poly.html#6.2>))

== Polynomial Evaluation
--- GSL::Poly.eval(c, x)
    Evaluates the polynomial ((|c[0] + c[1]x + c[2]x^2 + ...|)). 
    The polynomial coefficients ((|c|)) can be an (({Array})), 
    a (({GSL::Vector})), or an (({NArray})). The evaluation point ((|x|))
    is a (({Numeric})), (({Array})), (({GSL::Vector})) or (({NArray})).
  
    Ex)
       irb(main):001:0> require("rbgsl")
       => true
       irb(main):002:0> GSL::Poly.eval([1, 2, 3], 2)
       => 17.0
       irb(main):003:0> GSL::Poly.eval(GSL::Vector[1, 2, 3], 2)
       => 17.0
       irb(main):004:0> GSL::Poly.eval(NArray[1.0, 2, 3], 2)
       => 17.0
       irb(main):005:0> GSL::Poly.eval([1, 2, 3], [1, 2, 3])
       => [6.0, 17.0, 34.0]
       irb(main):006:0> GSL::Poly.eval([1, 2, 3], GSL::Vector[1, 2, 3])
       => GSL::Vector
       [ 6.000e+00 1.700e+01 3.400e+01 ]
       irb(main):007:0> GSL::Poly.eval([1, 2, 3], NArray[1.0, 2, 3])
       => NArray.float(3):
       [ 6.0, 17.0, 34.0 ]

== Solving polynomial equations
=== Quadratic Equations
--- GSL::Poly::solve_quadratic(a, b, c)
--- GSL::Poly::solve_quadratic([a, b, c])
    Find the real roots of the quadratic equation,
      a x^2 + b x + c = 0
    The coefficients are given by 3 numbers, or a Ruby array, 
    or a (({GSL::Vector})) object. The roots are returned as a (({GSL::Vector})).

    * Ex: z^2 - 3z + 2 = 0
        irb(main):006:0> GSL::Poly::solve_quadratic(1, -3, 2)
        => GSL::Vector: 
        [ 1.000e+00 2.000e+00 ]


--- GSL::Poly::complex_solve_quadratic(a, b, c)
--- GSL::Poly::complex_solve_quadratic([a, b, c])
    Find the complex roots of the quadratic equation,
      a z^2 + b z + z = 0
    The coefficients are given by 3 numbers or a Ruby array, or a 
    (({GSL::Vector})).
    The roots are returned as a (({GSL::Vector::Complex})) of two elements.
    
    * Ex: z^2 - 3z + 2 = 0
       irb(main):001:0> require("gsl")
       => true
       irb(main):002:0> GSL::Poly::complex_solve_quadratic(1, -3, 2)
       [ [1.000e+00 0.000e+00] [2.000e+00 0.000e+00] ]   
       => #<GSL::Vector::Complex:0x764014>
       irb(main):003:0> GSL::Poly::complex_solve_quadratic(1, -3, 2).real  <--- Real part
       => GSL::Vector::View: 
       [ 1.000e+00 2.000e+00 ]

=== Cubic Equations
--- GSL::Poly::solve_cubic(same as solve_quadratic)
    This method finds the real roots of the cubic equation,
      x^3 + a x^2 + b x + c = 0

--- GSL::Poly::complex_solve_cubic(same as solve_cubic)
    This method finds the complex roots of the cubic equation,
      z^3 + a z^2 + b z + c = 0

=== General Polynomial Equations
--- GSL::Poly::complex_solve(c0, c1, c2,,, )
--- GSL::Poly::solve(c0, c1, c2,,, )
    Find the complex roots of the polynomial equation. Note that 
    the coefficients are given by "ascending" order.

    * Ex: x^2 - 3 x + 2 == 0 
         irb(main):004:0> GSL::Poly::complex_solve(2, -3, 1)    <--- different from Poly::quadratic_solve
         [ [1.000e+00 0.000e+00] [2.000e+00 0.000e+00] ]
         => #<GSL::Vector::Complex:0x75e614>

== GSL::Poly Class
This class expresses polynomials of arbitrary orders.

=== Constructors
--- GSL::Poly.alloc(c0, c1, c2, ....)
--- GSL::Poly[c0, c1, c2, ....]
    This creates an instance of the (({GSL::Poly})) class, 
    which represents a polynomial
        c0 + c1 x + c2 x^2 + ....
    This class is derived from (({GSL::Vector})).

    * Ex: x^2 - 3 x + 2
        poly = GSL::Poly.alloc([2, -3, 1])

=== Instance Methods
--- GSL::Poly#eval(x)
--- GSL::Poly#at(x)
    Evaluates the polynomial 
    c[0] + c[1] x + c[2] x^2 + ... + c[len-1] x^{len-1} 
    using Horner's method for stability. The argument ((|x|)) is a
    (({Numeric})), (({GSL::Vector, Matrix})) or an (({Array})).

--- GSL::Poly#solve_quadratic
    Solve the quadratic equation.

    * Ex:  z^2 - 3 z + 2 = 0:
        irb(main):004:0> a = GSL::Poly[2, -3, 1]
        => GSL::Poly: 
        [ 2.000e+00 -3.000e+00 1.000e+00 ]
        irb(main):005:0> a.solve_quadratic
        => GSL::Vector: 
        [ 1.000e+00 2.000e+00 ]

--- GSL::Poly#solve_cubic
    Solve the cubic equation.

--- GSL::Poly#complex_solve
--- GSL::Poly#solve
--- GSL::Poly#roots
    These methods find the complex roots of the quadratic equation,
      c0 + c1 z + c2 z^2 + .... = 0

    * Ex: z^2 - 3 z + 2 = 0:
        irb(main):009:0> a = GSL::Poly[2, -3, 1]
        => GSL::Poly: 
        [ 2.000e+00 -3.000e+00 1.000e+00 ]
        irb(main):010:0> a.solve
        [ [1.000e+00 0.000e+00] [2.000e+00 0.000e+00] ]
        => #<GSL::Vector::Complex:0x35db28>

== Polynomial fitting
--- GSL::Poly.fit(x, y, order)
--- GSL::Poly.wfit(x, w, y, order)
    Finds the coefficient of a polynomial of order ((|order|)) 
    that fits the vector data (((|x, y|))) in a least-square sense.
    This provides a higher-level interface to the method
    ((<GSL::Multifit#linear|URL:fit.html>)) in a case of polynomial fitting.

    Example:
      #!/usr/bin/env ruby
      require("rbgsl")

      x = GSL::Vector[1, 2, 3, 4, 5]
      y = GSL::Vector[5.5, 43.1, 128, 290.7, 498.4]
      # The results are stored in a polynomial "coef"
      coef, cov, chisq, status = Poly.fit(x, y, 3) 

      x2 = GSL::Vector.linspace(1, 5, 20)
      graph([x, y], [x2, coef.eval(x2)], "-C -g 3 -S 4")

== Divided-difference representations

--- GSL::Poly::dd_init(xa, ya)
    This method computes a divided-difference representation of the 
    interpolating polynomial for the points ((|(xa, ya)|)).

--- GSL::Poly::DividedDifference#eval(x)
    This method evaluates the polynomial stored in divided-difference form 
    ((|self|)) at the point ((|x|)).

--- GSL::Poly::DividedDifference#taylor(xp)
    This method converts the divided-difference representation of a polynomial 
    to a Taylor expansion. On output the Taylor coefficients of the polynomial 
    expanded about the point ((|xp|)) are returned.

== Extensions
=== Special Polynomials
--- GSL::Poly.hermite(n)
    This returns coefficients of the ((|n|))-th order Hermite polynomial, ((|H(x; n)|)).
    For order of ((|n|)) >= 3, this method uses the recurrence relation
         H(x; n+1) = 2 x H(x; n) - 2 n H(x; n-1)
    * Ex:
        irb(main):013:0> GSL::Poly.hermite(2)
        => GSL::Poly::Int: 
        [ -2 0 4 ]                            <----- 4x^2 - 2
        irb(main):014:0> GSL::Poly.hermite(5)
        => GSL::Poly::Int: 
        [ 0 120 0 -160 0 32 ]                 <----- 32x^5 - 160x^3 + 120x
        irb(main):015:0> GSL::Poly.hermite(7)
        => GSL::Poly::Int: 
        [ 0 -1680 0 3360 0 -1344 0 128 ]

--- GSL::Poly.cheb(n)
--- GSL::Poly.chebyshev(n)
    Return the coefficients of the ((|n|))-th order Chebyshev polynomial, ((|T(x; n|)).
    For order of ((|n|)) >= 3, this method uses the recurrence relation
         T(x; n+1) = 2 x T(x; n) - T(x; n-1)

--- GSL::Poly.cheb_II(n)
--- GSL::Poly.chebyshev_II(n)
    Return the coefficients of the ((|n|))-th order Chebyshev polynomial of type II,
    ((|U(x; n|)).
         U(x; n+1) = 2 x U(x; n) - U(x; n-1)

--- GSL::Poly.bell(n)
    Bell polynomial

--- GSL::Poly.bessel(n)
    Bessel polynomial

--- GSL::Poly.laguerre(n)
    Retunrs the coefficients of the ((|n|))-th order Laguerre polynomial
    multiplied by n!.

    Ex:
         rb(main):001:0> require("rbgsl")
         => true
         irb(main):002:0> GSL::Poly.laguerre(0)
         => GSL::Poly::Int: 
         [ 1 ]                              <--- 1
         irb(main):003:0> GSL::Poly.laguerre(1)
         => GSL::Poly::Int: 
         [ 1 -1 ]                           <--- -x + 1
         irb(main):004:0> GSL::Poly.laguerre(2)
         => GSL::Poly::Int: 
         [ 2 -4 1 ]                         <--- (x^2 - 4x + 2)/2!
         irb(main):005:0> GSL::Poly.laguerre(3)
         => GSL::Poly::Int:                   
         [ 6 -18 9 -1 ]                     <--- (-x^3 + 9x^2 - 18x + 6)/3!
         irb(main):006:0> GSL::Poly.laguerre(4)
         => GSL::Poly::Int:  
         [ 24 -96 72 -16 1 ]                <--- (x^4 - 16x^3 + 72x^2 - 96x + 24)/4!
 
=== Polynomial Operations
--- GSL::Poly#conv
--- GSL::Poly#deconv
--- GSL::Poly#reduce
--- GSL::Poly#deriv
--- GSL::Poly#integ
--- GSL::Poly#compan

((<prev|URL:complex.html>))
((<next|URL:sf.html>))

((<Reference index|URL:ref.html>))
((<top|URL:index.html>))

=end

