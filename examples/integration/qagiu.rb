#!/usr/bin/env ruby
require("gsl")
include GSL
include Math

printf("QAGIU (integrate [0:+infty])\n\n")

w = Integration::Workspace.alloc(1000)
xmin = 0.0

f1 = Function.alloc{ |x, a|
  1.0/(pow_4(x) + pow_4(a))
}

printf("Case 1: f(x; a) = 1/(x^4 + a^4), I(a) = pi/2sqrt2/a^3\n")
a = 1.0
printf("        Expected: I(1) = %10.9f\n", M_PI/2/M_SQRT2)
f1.set_params(a)
printf("        QAGIU Result:    %10.9f\n\n", f1.qagiu(xmin, w)[0])

f2 = Function.alloc{ |x, a|
  x*x/(pow_4(x) + pow_4(a))
}

printf("Case 2: f(x; a) = 1/(x^4 + a^4), I(a) = pi/2sqrt2/a\n")
a = 2.0
printf("        Expected: I(2) = %10.9f\n", M_PI/2/M_SQRT2/a)
f2.set_params(a)
printf("        QAGIU Result:    %10.9f\n\n", f2.qagiu(xmin, w)[0])
