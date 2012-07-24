#!/usr/bin/env ruby
require("gsl")
include GSL
include Math

printf("QAGI (integrate [-infty:+infty])\n\n")

w = Integration::Workspace.alloc(1000)

f = Function.alloc{ |x, a|
  cos(x)/(x*x + a*a)
}

a = 1.0
f.set_params(a)
printf("Case 1: f(x; a) = cos(x)/(x^2 + a^2), I(a) = pi/a e^(-a)\n")
printf("        Expected: I(1) = %10.9f\n", M_PI*exp(-a))
printf("        QAGI Result:     %10.9f\n", f.qagi([0, 1e-4], w)[0])

f = Function.alloc{ |x, a|
  exp(a*x)/(1 + exp(x))
}

a = 0.3
f.set_params(a)
printf("Case 2: f(x; a) = e^{ax}/(1 + e^x), I(a) = pi/sin(a pi)\n")
printf("        Expected: I(0.3) = %10.9f\n", M_PI/sin(a*M_PI))
printf("        QAGI Result:       %10.9f\n", f.qagi([0, 1e-2], w)[0])
