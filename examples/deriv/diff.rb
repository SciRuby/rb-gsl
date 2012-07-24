#!/usr/bin/env ruby
require("gsl")

f = GSL::Function.alloc { |x|
  GSL::pow(x, 1.5)
}

printf("f(x) = x^(3/2)\n");

x = 2.0
result, abserr,status = f.diff_central(x)
printf("x = 2.0\n");
printf("f'(x) = %.10f +/- %.5f\n", result, abserr);
printf("exact = %.10f\n\n", 1.5 * Math::sqrt(2.0));

x = 0.0
result, abserr, status = f.diff_forward(x)
printf("x = 0.0\n");
printf("f'(x) = %.10f +/- %.5f\n", result, abserr);
printf("exact = %.10f\n", 0.0);

f2 = GSL::Function.alloc { |x, a| Math::exp(a*x) }

f2.set_params(2)
p f2.diff_central(0)

f2.set_params(3)
p f2.diff_central(0)

f2.set_params(123)
p f2.diff_central(0)

p GSL::Diff.central(f2, 0)
p GSL::Diff.forward(f2, 0)
p GSL::Diff.backward(f2, 0)
