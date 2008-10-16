=begin
= GSL::Function class

== Class Methods

--- GSL::Function.alloc
    Constructor.
   
    * ex:
       require("gsl")
       f = GSL::Function.alloc { |x| sin(x) }

    The value of the function is calculated by the method (({Function#eval})), as

       p f.eval(x)

    The function can have parameters of arbitrary numbers. Here is an
    example in case of exponential function (({f(x; a, b) = a*exp(-b*x)})).
      
       f = GSL::Function.alloc { |x, params|    # x: a scalar, params: an array
         a = params[0]; b = params[1]
         a*exp(-b*x)
       }
    To evaluate the function (({f(x) = 2*exp(-3*x)})), 
       f.set_params([2, 3])
       f.eval(x)

== Methods

--- GSL::Function#eval(x)
--- GSL::Function#call(x)
--- GSL::Function#at(x)
--- GSL::Function#[x]
    These methods return a value of the function at ((|x|)).
      p f.eval(2.5)
      p f.call(2.5)
      p f[2.5]
    The argument ((|x|)) can be a scalar, a Vector, Matrix, Array or Range.

--- GSL::Function#set { |x| ... }
--- GSL::Function#set(proc, params)
    This method sets or resets the procedure of ((|self|)), as

      f = GSL::Function.alloc { |x| sin(x) }
      p f.eval(1.0)               <- sin(1.0)
      f.set { |x| cos(x) }
      p f.eval(1.0)               <- cos(1.0)

--- GSL::Function#set_params(params)
    This set the constant parameters of the function.

== Graph
--- GSL::Function#graph(x[, options])
    This method uses (({GNU graph})) to plot the function ((|self|)).
    The argument ((|x|)) is given by a (({GSL::Vector})) or an (({Array})).

    Ex: Plot sin(x)
         f = Function.alloc { |x| Math::sin(x) }
         x = Vector.linspace(0, 2*M_PI, 50)
         f.graph(x, "-T X -g 3 -C -L 'sin(x)'")


== Example
A quadratic function, f(x) = x^2 + 2x + 3.

    irb(main):001:0> require("gsl")
    => true
    irb(main):002:0> f = Function.alloc { |x, param| x*x + param[0]*x + param[1] } 
    => #<GSL::Function:0x6e8eb0>
    irb(main):003:0> f.set_params(2, 3)
    => #<GSL::Function:0x6e8eb0>
    irb(main):004:0> f.eval(2)                             <--- Scalar
    => 11
    irb(main):005:0> f.eval(1..4)                          <--- Range
    => [6.0, 11.0, 18.0, 27.0]
    irb(main):006:0> f.eval([1, 2, 3])                     <--- Array
    => [6.0, 11.0, 18.0]
    irb(main):007:0> f.eval(Matrix.alloc([1, 2], [3, 4]))    <--- GSL::Matrix
    [ 6.000e+00 1.100e+01 
      1.800e+01 2.700e+01 ]
    => #<GSL::Matrix:0x6dd1b4>

((<back|URL:index.html>))
=end
