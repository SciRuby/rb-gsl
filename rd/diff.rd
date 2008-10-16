=begin
= Numerical Differentiation
The functions described in this chapter compute numerical derivatives by 
finite differencing. An adaptive algorithm is used to find the best choice 
of finite difference and to estimate the error in the derivative. 

Contentes:
(1) ((<Deriv methods|URL:diff.html#1>))
(2) ((<Diff methods|URL:diff.html#2>))

== Deriv methods (for GSL 1.4.90 or later)
Numerical derivatives should now be calculated using the
(({GSL::Deriv.forward, GSL::Deriv.central})) and (({GSL::Deriv.backward})) methods,
which accept a step-size argument in addition to the position x.  The
original (({GSL::Diff})) methods (without the step-size) are deprecated.

--- GSL::Deriv.central(f, x, h = 1e-8)
--- GSL::Function#deriv_central(x, h = 1e-8)
    These methods compute the numerical derivative of the function ((|f|)) 
    at the point ((|x|)) using an adaptive central difference algorithm with a 
    step-size of ((|h|)). If a scalar ((|x|)) is given, the derivative and an 
    estimate of its absolute error are returned as an array, [((|result, abserr, status|))].
    If a vector/matrix/array ((|x|)) is given, an array of two elements are returned,
    [((|result, abserr|))], here each them is also a vector/matrix/array of the same
    dimension of ((|x|)).

    The initial value of ((|h|)) is used to estimate an optimal step-size, 
    based on the scaling of the truncation error and round-off error in the 
    derivative calculation. The derivative is computed using a 5-point rule for 
    equally spaced abscissae at x-h, x-h/2, x, x+h/2, x, with an error estimate 
    taken from the difference between the 5-point rule and the corresponding 3-point 
    rule x-h, x, x+h. Note that the value of the function at x does not contribute 
    to the derivative calculation, so only 4-points are actually used.

--- GSL::Deriv.forward(f, x, h = 1e-8)
--- GSL::Function#deriv_forward(x, h = 1e-8)
    These methods compute the numerical derivative of the function ((|f|)) at 
    the point ((|x|)) using an adaptive forward difference algorithm with a step-size 
    of ((|h|)). The function is evaluated only at points greater than ((|x|)), 
    and never at ((|x|)) itself. The derivative and an estimate of its absolute error 
    are returned as an array, [((|result, abserr|))]. 
    These methods should be used if f(x) has a 
    discontinuity at ((|x|)), or is undefined for values less than ((|x|)).

    The initial value of ((|h|)) is used to estimate an optimal step-size, based on the 
    scaling of the truncation error and round-off error in the derivative calculation. 
    The derivative at x is computed using an "open" 4-point rule for equally spaced 
    abscissae at x+h/4, x+h/2, x+3h/4, x+h, with an error estimate taken from the 
    difference between the 4-point rule and the corresponding 2-point rule x+h/2, x+h.

--- GSL::Deriv.backward(f, x, h)
--- GSL::Function#deriv_backward(x, h)
    These methods compute the numerical derivative of the function ((|f|)) at the 
    point ((|x|)) using an adaptive backward difference algorithm with a step-size 
    of ((|h|)). The function is evaluated only at points less than ((|x|)), 
    and never at ((|x|)) itself. The derivative and an estimate of its absolute error 
    are returned as an array, [((|result, abserr|))]. 
    This function should be used if f(x) has a discontinuity at ((|x|)), 
    or is undefined for values greater than ((|x|)).

    These methods are equivalent to calling the method (({forward})) 
    with a negative step-size.

== Diff Methods (obsolete)

--- GSL::Diff.central(f, x)
--- GSL::Function#diff_central(x)
    These compute the numerical derivative of the function ((|f|)) ( ((<GSL::Function|URL:function.html>)) object) at the point ((|x|)) 
    using an adaptive central difference algorithm. The result is returned as an array 
    which contains the derivative and an estimate of its absolute error.

--- GSL::Diff.forward(f, x)
--- GSL::Function#diff_forward(x)
    These compute the numerical derivative of the function at the point x using an adaptive forward difference algorithm. 

--- GSL::Diff.backward(f, x)
--- GSL::Function#diff_backward(x)
    These compute the numerical derivative of the function at the point x using an adaptive backward difference algorithm. 

== Example

     #!/usr/bin/env ruby
     require "rbgsl"

     f = GSL::Function.alloc { |x|
       pow(x, 1.5)
     }

     printf ("f(x) = x^(3/2)\n");

     x = 2.0
     h = 1e-8
     result, abserr = f.deriv_central(x, h)
     printf("x = 2.0\n");
     printf("f'(x) = %.10f +/- %.10f\n", result, abserr);
     printf("exact = %.10f\n\n", 1.5 * Math::sqrt(2.0));

     x = 0.0
     result, abserr = Deriv.forward(f, x, h) # equivalent to f.deriv_forward(x, h)
     printf("x = 0.0\n");
     printf("f'(x) = %.10f +/- %.10f\n", result, abserr);
     printf("exact = %.10f\n", 0.0);

The results are

     f(x) = x^(3/2)
     x = 2.0
     f'(x) = 2.1213203120 +/- 0.0000004064
     exact = 2.1213203436

     x = 0.0
     f'(x) = 0.0000000160 +/- 0.0000000339
     exact = 0.0000000000

((<prev|URL:interp.html>))
((<next|URL:cheb.html>))

((<Reference index|URL:ref.html>))
((<top|URL:index.html>))

=end
