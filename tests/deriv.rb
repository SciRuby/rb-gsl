#!/usr/bin/env ruby
# Ruby/GSL implementation of GSL "deriv/test.c"
require("gsl")
require("./gsl_test.rb")
include GSL::Test

include Math

f1 = GSL::Function.alloc { |x| exp(x) }
df1 = GSL::Function.alloc { |x| exp(x) }

f2 = GSL::Function.alloc { |x| 
  if x >= 0.0; x*sqrt(x);
  else; 0.0; end
}
df2 = GSL::Function.alloc { |x| 
  if x >= 0.0; 1.5*sqrt(x);
  else; 0.0; end
}

f3 = GSL::Function.alloc { |x| 
  if x != 0.0; sin(1.0/x);
  else; 0.0; end
}
df3 = GSL::Function.alloc { |x| 
  if x != 0.0; -cos(1.0/x)/(x*x);
  else; 0.0; end
}

f4 = GSL::Function.alloc { |x| exp(-x*x) }
df4 = GSL::Function.alloc { |x| -2.0*x*exp(-x*x) }

f5 = GSL::Function.alloc { |x| x*x }
df5 = GSL::Function.alloc { |x| 2.0*x }

f6 = GSL::Function.alloc { |x| 1.0/x }
df6 = GSL::Function.alloc { |x| -1.0/(x*x) }

def test_deriv(deriv, f, df, x, desc)
  expected = df.eval(x)
  h = 1e-4
  case deriv
  when "central"
    result, abserr = f.deriv_central(x, h)
  when "forward"
    result, abserr = f.deriv_forward(x, h)
  when "backward"
    result, abserr = f.deriv_backward(x, h)
  else
    raise("undefined operation")
  end
  GSL::Test::test_abs(result, expected, GSL::MIN(h, expected.abs) + GSL::DBL_EPSILON, desc)
  if abserr < (result - expected).abs
    GSL::Test::test_factor(abserr, (result - expected).abs, 2, desc + " error estimate")
  elsif result == expected or expected == 0.0
    GSL::Test::test_abs(abserr, 0.0, 1e-6, desc + " abserr")
  else
    d = (result - expected).abs
    GSL::Test::test_abs(abserr, (result - expected).abs, 1e6*d, desc + " abserr")
  end
end

test_deriv("central", f1, df1, 1.0, "exp(x), x=1, central deriv")
test_deriv("forward", f1, df1, 1.0, "exp(x), x=1, forward deriv")
test_deriv("backward", f1, df1, 1.0, "exp(x), x=1, backward deriv")

test_deriv("central", f2, df2, 0.1, "x^(3/2), x=0.1, central deriv")
test_deriv("forward", f2, df2, 0.1, "x^(3/2), x=0.1, forward deriv")
test_deriv("backward", f2, df2, 0.1, "x^(3/2), x=0.1, backward deriv")

test_deriv("central", f3, df3, 0.45, "sin(1/x), x=0.45, central deriv")
test_deriv("forward", f3, df3, 0.45, "sin(1/x), x=0.45, forward deriv")
test_deriv("backward", f3, df3, 0.45, "sin(1/x), x=0.45, backward deriv")

test_deriv("central", f4, df4, 0.5, "exp(-x^2), x=0.5, central deriv")
test_deriv("forward", f4, df4, 0.5, "exp(-x^2), x=0.5, forward deriv")
test_deriv("backward", f4, df4, 0.5, "exp(-x^2), x=0.5, backward deriv")

test_deriv("central", f5, df5, 0.0, "x^2, x=0, central deriv")
test_deriv("forward", f5, df5, 0.0, "x^2, x=0, forward deriv")
test_deriv("backward", f5, df5, 0.0, "x^2, x=0, backward deriv")

test_deriv("central", f6, df6, 10.0, "1/x, x=10, central deriv")
test_deriv("forward", f6, df6, 10.0, "1/x, x=10, forward deriv")
test_deriv("backward", f6, df6, 10.0, "1/x, x=10, backward deriv")
