#!/usr/bin/env ruby
require("gsl")
require("./gsl_test2.rb")
include GSL::Test
exit unless GSL::Poly.method_defined?("complex_solve_quartic")

z = GSL::Poly.complex_solve_quartic(0.0, 0.0, 0.0, -81.0)
puts("Four roots, x^4 - 81")
GSL::Test::test_rel(z[0].re, -3.0, 1e-9, "z0.real")
GSL::Test::test_rel(z[0].im,  0.0, 1e-9, "z0.imag")
GSL::Test::test_rel(z[1].re,  0.0, 1e-9, "z1.real")
GSL::Test::test_rel(z[1].im, -3.0, 1e-9, "z1.imag")
GSL::Test::test_rel(z[2].re,  0.0, 1e-9, "z2.real")
GSL::Test::test_rel(z[2].im,  3.0, 1e-9, "z2.imag")
GSL::Test::test_rel(z[3].re,  3.0, 1e-9, "z3.real")
GSL::Test::test_rel(z[3].im,  0.0, 1e-9, "z3.imag")

sol = 3.0/Math.sqrt(2.0)
z = GSL::Poly.complex_solve_quartic(0.0, 0.0, 0.0, 81.0)
puts("Four roots, x^4 + 81")
GSL::Test::test_rel(z[0].re, -sol, 1e-9, "z0.real")
GSL::Test::test_rel(z[0].im, -sol, 1e-9, "z0.imag")
GSL::Test::test_rel(z[1].re, -sol, 1e-9, "z1.real")
GSL::Test::test_rel(z[1].im,  sol, 1e-9, "z1.imag")
GSL::Test::test_rel(z[2].re,  sol, 1e-9, "z2.real")
GSL::Test::test_rel(z[2].im, -sol, 1e-9, "z2.imag")
GSL::Test::test_rel(z[3].re,  sol, 1e-9, "z3.real")
GSL::Test::test_rel(z[3].im,  sol, 1e-9, "z3.imag")

