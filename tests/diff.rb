#!/usr/bin/env ruby
# Ruby/GSL implementation of GSL "diff/test.c"
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

def test_diff(diff, f, df, x, desc)
  expected = df.eval(x)
  case diff
  when "central"
    result, abserr = f.diff_central(x)
  when "forward"
    result, abserr = f.diff_forward(x)
  when "backward"
    result, abserr = f.diff_backward(x)
  else
    raise("undefined operation")
  end
  GSL::Test::test_abs(result, expected, abserr, desc)
  desc2 = sprintf("%s, valid error estimate", desc)
  GSL::Test::test((result - expected).abs > abserr, desc2)
end

test_diff("central", f1, df1, 1.0, "exp(x), x=1, central diff")
test_diff("forward", f1, df1, 1.0, "exp(x), x=1, forward diff")
test_diff("backward", f1, df1, 1.0, "exp(x), x=1, backward diff")

test_diff("central", f2, df2, 0.1, "x^(3/2), x=0.1, central diff")
test_diff("forward", f2, df2, 0.1, "x^(3/2), x=0.1, forward diff")
test_diff("backward", f2, df2, 0.1, "x^(3/2), x=0.1, backward diff")

test_diff("central", f3, df3, 0.45, "sin(1/x), x=0.45, central diff")
test_diff("forward", f3, df3, 0.45, "sin(1/x), x=0.45, forward diff")
test_diff("backward", f3, df3, 0.45, "sin(1/x), x=0.45, backward diff")

test_diff("central", f4, df4, 0.5, "exp(-x^2), x=0.5, central diff")
test_diff("forward", f4, df4, 0.5, "exp(-x^2), x=0.5, forward diff")
test_diff("backward", f4, df4, 0.5, "exp(-x^2), x=0.5, backward diff")

test_diff("central", f5, df5, 0.0, "x^2, x=0, central diff")
test_diff("forward", f5, df5, 0.0, "x^2, x=0, forward diff")
test_diff("backward", f5, df5, 0.0, "x^2, x=0, backward diff")

test_diff("central", f6, df6, 10.0, "1/x, x=10, central diff")
test_diff("forward", f6, df6, 10.0, "1/x, x=10, forward diff")
test_diff("backward", f6, df6, 10.0, "1/x, x=10, backward diff")
