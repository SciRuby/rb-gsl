#!/usr/bin/env ruby
require("gsl")
require("../gsl_test2.rb")
include GSL::Test
include Math

eps = 100.0 * GSL::DBL_EPSILON
GSL::IEEE::env_setup()

c = GSL::Poly.alloc(1.0, 0.5, 0.3)
x = 0.5
y = c.eval(x)
GSL::Test::test_rel(y, 1 + 0.5 * x + 0.3 * x * x, eps,
                    "gsl_poly_eval({1, 0.5, 0.3}, 0.5)")

d = GSL::Poly.alloc( 1, -1, 1, -1, 1, -1, 1, -1, 1, -1, 1)
x = 1.0
#y = d.eval(x)
y = GSL::Poly.eval(d, x)
#y = GSL::Poly.eval(d, 11, x)
GSL::Test::test_rel(y, 1.0, eps,
                    "gsl_poly_eval({1,-1, 1, -1, 1, -1, 1, -1, 1, -1, 1}, 1.0)")

x0, x1 = GSL::Poly.solve_quadratic(4.0, -20.0, 25.0).to_a
GSL::Test::test_rel(x0, 2.5, 1e-9, "x0, (2x - 5)^2 = 0")
GSL::Test::test_rel(x1, 2.5, 1e-9, "x1, (2x - 5)^2 = 0")

x0, x1 = GSL::Poly.solve_quadratic(4.0, 7.0, 0.0).to_a
GSL::Test::test_rel(x0, -1.75, 1e-9, "x0, x(4x + 7) = 0")
GSL::Test::test_rel(x1, 0.0, 1e-9, "x1, x(4x + 7) = 0")

x0, x1 = GSL::Poly.solve_quadratic(5.0, 0.0, -20.0).to_a
test_rel(x0, -2.0, 1e-9, "x0, 5 x^2 = 20")
test_rel(x1, 2.0, 1e-9, "x1, 5 x^2 = 20")

x0, x1, x2 = GSL::Poly.solve_cubic(0.0, 0.0, -27.0).to_a
test_rel(x0, 3.0, 1e-9, "x0, x^3 = 27")

x0, x1, x2 = GSL::Poly.solve_cubic(-51.0, 867.0, -4913.0).to_a
test_rel(x0, 17.0, 1e-9, "x0, (x-17)^3=0")
test_rel(x1, 17.0, 1e-9, "x1, (x-17)^3=0")
test_rel(x2, 17.0, 1e-9, "x2, (x-17)^3=0")

x0, x1, x2 = GSL::Poly.solve_cubic(-57.0, 1071.0, -6647.0).to_a
test_rel(x0, 17.0, 1e-9, "x0, (x-17)(x-17)(x-23)=0")
test_rel(x1, 17.0, 1e-9, "x1, (x-17)(x-17)(x-23)=0")
test_rel(x2, 23.0, 1e-9, "x2, (x-17)(x-17)(x-23)=0")

x0, x1, x2 = GSL::Poly.solve_cubic(-11.0, -493.0, +6647.0).to_a
test_rel(x0, -23.0, 1e-9, "x0, (x+23)(x-17)(x-17)=0")
test_rel(x1, 17.0, 1e-9, "x1, (x+23)(x-17)(x-17)=0")
test_rel(x2, 17.0, 1e-9, "x2, (x+23)(x-17)(x-17)=0")

x0, x1, x2 = GSL::Poly.solve_cubic(-143.0, 5087.0, -50065.0).to_a
test_rel(x0, 17.0, 1e-9, "x0, (x-17)(x-31)(x-95)=0")
test_rel(x1, 31.0, 1e-9, "x1, (x-17)(x-31)(x-95)=0")
test_rel(x2, 95.0, 1e-9, "x2, (x-17)(x-31)(x-95)=0")

x0, x1, x2 = GSL::Poly.solve_cubic(-109.0, 803.0, 50065.0).to_a
test_rel(x0, -17.0, 1e-9, "x0, (x+17)(x-31)(x-95)=0")
test_rel(x1, 31.0, 1e-9, "x1, (x+17)(x-31)(x-95)=0")
test_rel(x2, 95.0, 1e-9, "x2, (x+17)(x-31)(x-95)=0")

#z0, z1 = GSL::Poly.complex_solve_quadratic(4.0, -20.0, 26.0).to_a
r = GSL::Poly::Complex.solve_quadratic(4.0, -20.0, 26.0)
z0 = r[0]
z1 = r[1]
test_rel(z0.re, 2.5, 1e-9, "z0.real, (2x - 5)^2 = -1")
test_rel(z0.im, -0.5, 1e-9, "z0.imag, (2x - 5)^2 = -1")
test_rel(z1.re, 2.5, 1e-9, "z1.real, (2x - 5)^2 = -1")
test_rel(z1.im, 0.5, 1e-9, "z1.imag, (2x - 5)^2 = -1")

z = GSL::Poly.complex_solve_quadratic(4.0, -20.0, 25.0)
test_rel(z[0].re, 2.5, 1e-9, "z0.real, (2x - 5)^2 = 0")
test_rel(z[0].im, 0.0, 1e-9, "z0.imag (2x - 5)^2 = 0")
test_rel(z[1].re, 2.5, 1e-9, "z1.real, (2x - 5)^2 = 0")
test_rel(z[1].im, 0.0, 1e-9, "z1.imag (2x - 5)^2 = 0")
test(z[0].re != z[1].re ? 1 : 0,
                        "z0.real == z1.real, (2x - 5)^2 = 0")
test(z[1].im != z[1].im ? 1 : 0,
                        "z0.imag == z1.imag, (2x - 5)^2 = 0")

z = GSL::Poly.complex_solve_quadratic(4.0, -20.0, 21.0)
test_rel(z[0].re, 1.5, 1e-9, "z0.real, (2x - 5)^2 = 4")
test_rel(z[0].im, 0.0, 1e-9, "z0.imag, (2x - 5)^2 = 4")
test_rel(z[1].re, 3.5, 1e-9, "z1.real, (2x - 5)^2 = 4")
test_rel(z[1].im, 0.0, 1e-9, "z1.imag, (2x - 5)^2 = 4")

z = GSL::Poly.complex_solve_quadratic(4.0, 7.0, 0.0)
test_rel(z[0].re, -1.75, 1e-9, "z[0].real, x(4x + 7) = 0")
test_rel(z[0].im, 0.0, 1e-9, "z[0].imag, x(4x + 7) = 0")
test_rel(z[1].re, 0.0, 1e-9, "z[1].real, x(4x + 7) = 0")
test_rel(z[1].im, 0.0, 1e-9, "z[1].imag, x(4x + 7) = 0")

z =GSL::Poly.complex_solve_quadratic(5.0, 0.0, -20.0)
test_rel(z[0].re, -2.0, 1e-9, "z[0].real, 5 x^2 = 20")
test_rel(z[0].im, 0.0, 1e-9, "z[0].imag, 5 x^2 = 20")
test_rel(z[1].re, 2.0, 1e-9, "z[1].real, 5 x^2 = 20")
test_rel(z[1].im, 0.0, 1e-9, "z[1].imag, 5 x^2 = 20")

z = GSL::Poly.complex_solve_quadratic(5.0, 0.0, 20.0)
test_rel(z[0].re, 0.0, 1e-9, "z[0].real, 5 x^2 = -20")
test_rel(z[0].im, -2.0, 1e-9, "z[0].imag, 5 x^2 = -20")
test_rel(z[1].re, 0.0, 1e-9, "z[1].real, 5 x^2 = -20")
test_rel(z[1].im, 2.0, 1e-9, "z[1].imag, 5 x^2 = -20")

z = GSL::Poly.complex_solve_cubic(0.0, 0.0, -27.0)
test_rel(z[0].re, -1.5, 1e-9, "z[0].real, x^3 = 27");
test_rel(z[0].im, -1.5 * sqrt(3.0), 1e-9,
              "z[0].imag, x^3 = 27");
test_rel(z[1].re, -1.5, 1e-9, "z[1].real, x^3 = 27");
test_rel(z[1].im, 1.5 * sqrt(3.0), 1e-9, "z[1].imag, x^3 = 27");
test_rel(z[2].re, 3.0, 1e-9, "z[2].real, x^3 = 27");
test_rel(z[2].im, 0.0, 1e-9, "z[2].imag, x^3 = 27")

z = GSL::Poly.complex_solve_cubic(-1.0, 1.0, 39.0)
test_rel(z[0].re, -3.0, 1e-9, "z[0].real, (x+3)(x^2+1) = 0");
test_rel(z[0].im, 0.0, 1e-9, "z[0].imag, (x+3)(x^2+1) = 0");
test_rel(z[1].re, 2.0, 1e-9, "z[1].real, (x+3)(x^2+1) = 0");
test_rel(z[1].im, -3.0, 1e-9, "z[1].imag, (x+3)(x^2+1) = 0");
test_rel(z[2].re, 2.0, 1e-9, "z[2].real, (x+3)(x^2+1) = 0");
test_rel(z[2].im, 3.0, 1e-9, "z[2].imag, (x+3)(x^2+1) = 0")

z = GSL::Poly.complex_solve_cubic(-51.0, 867.0, -4913.0)
test_rel(z[0].re, 17.0, 1e-9, "z[0].real, (x-17)^3=0");
test_rel(z[0].im, 0.0, 1e-9, "z[0].imag, (x-17)^3=0");
test_rel(z[1].re, 17.0, 1e-9, "z[1].real, (x-17)^3=0");
test_rel(z[1].im, 0.0, 1e-9, "z[1].imag, (x-17)^3=0");
test_rel(z[2].re, 17.0, 1e-9, "z[2].real, (x-17)^3=0");
test_rel(z[2].im, 0.0, 1e-9, "z[2].imag, (x-17)^3=0")

z = GSL::Poly.complex_solve_cubic(-57.0, 1071.0, -6647.0)
test_rel(z[0].re, 17.0, 1e-9, "z[0].real, (x-17)(x-17)(x-23)=0");
test_rel(z[0].im, 0.0, 1e-9, "z[0].imag, (x-17)(x-17)(x-23)=0");
test_rel(z[1].re, 17.0, 1e-9, "z[1].real, (x-17)(x-17)(x-23)=0");
test_rel(z[1].im, 0.0, 1e-9, "z[1].imag, (x-17)(x-17)(x-23)=0");
test_rel(z[2].re, 23.0, 1e-9, "z[2].real, (x-17)(x-17)(x-23)=0");
test_rel(z[2].im, 0.0, 1e-9, "z[2].imag, (x-17)(x-17)(x-23)=0")

z = GSL::Poly.complex_solve_cubic(-11.0, -493.0, +6647.0)
test_rel(z[0].re, -23.0, 1e-9,
          "z[0].real, (x+23)(x-17)(x-17)=0");
test_rel(z[0].im, 0.0, 1e-9, "z[0].imag, (x+23)(x-17)(x-17)=0");
test_rel(z[1].re, 17.0, 1e-9, "z[1].real, (x+23)(x-17)(x-17)=0");
test_rel(z[1].im, 0.0, 1e-9, "z[1].imag, (x+23)(x-17)(x-17)=0");
test_rel(z[2].re, 17.0, 1e-9, "z[2].real, (x+23)(x-17)(x-17)=0");
test_rel(z[2].im, 0.0, 1e-9, "z[2].imag, (x+23)(x-17)(x-17)=0");

z = GSL::Poly.complex_solve_cubic(-143.0, 5087.0, -50065.0)
test_rel(z[0].re, 17.0, 1e-9, "z[0].real, (x-17)(x-31)(x-95)=0");
test_rel(z[0].im, 0.0, 1e-9, "z[0].imag, (x-17)(x-31)(x-95)=0");
test_rel(z[1].re, 31.0, 1e-9, "z[1].real, (x-17)(x-31)(x-95)=0");
test_rel(z[1].im, 0.0, 1e-9, "z[1].imag, (x-17)(x-31)(x-95)=0");
test_rel(z[2].re, 95.0, 1e-9, "z[2].real, (x-17)(x-31)(x-95)=0");
test_rel(z[2].im, 0.0, 1e-9, "z[2].imag, (x-17)(x-31)(x-95)=0")

a = GSL::Poly.alloc(-120, 274, -225, 85, -15, 1)
w = GSL::Poly::Complex::Workspace.alloc(a.size)
z = GSL::Poly.complex_solve(a, 6, w)
#z = GSL::Poly.complex_solve(a, w)
#z = GSL::Poly.complex_solve(a)
test_rel(z[0].re, 1.0, 1e-9, "z[0].real, 5th-order polynomial");
test_rel(z[0].im, 0.0, 1e-9, "z[0].imag, 5th-order polynomial");
test_rel(z[1].re, 2.0, 1e-9, "z[1].real, 5th-order polynomial");
test_rel(z[1].im, 0.0, 1e-9, "z[1].imag, 5th-order polynomial");
test_rel(z[2].re, 3.0, 1e-9, "z[2].real, 5th-order polynomial");
test_rel(z[2].im, 0.0, 1e-9, "z[2].imag, 5th-order polynomial");
test_rel(z[3].re, 4.0, 1e-9, "z3.real, 5th-order polynomial");
test_rel(z[3].im, 0.0, 1e-9, "z3.imag, 5th-order polynomial");
test_rel(z[4].re, 5.0, 1e-9, "z4.real, 5th-order polynomial");
test_rel(z[4].im, 0.0, 1e-9, "z4.imag, 5th-order polynomial")

a = GSL::Poly.alloc(1, 0, 0, 0, 1, 0, 0, 0, 1)
w = GSL::Poly::Complex::Workspace.alloc(a.size)
c = 0.5
s = sqrt(3)/2
z = GSL::Poly.complex_solve(a, w)
test_rel(z[0].re, -s, 1e-9, "z[0].real, 8th-order polynomial");
test_rel(z[0].im, c, 1e-9, "z[0].imag, 8th-order polynomial");
test_rel(z[1].re, -s, 1e-9, "z[1].real, 8th-order polynomial");
test_rel(z[1].im, -c, 1e-9, "z[1].imag, 8th-order polynomial");
test_rel(z[2].re, -c, 1e-9, "z[2].real, 8th-order polynomial");
test_rel(z[2].im, s, 1e-9, "z[2].imag, 8th-order polynomial");
test_rel(z[3].re, -c, 1e-9, "z3.real, 8th-order polynomial");
test_rel(z[3].im, -s, 1e-9, "z3.imag, 8th-order polynomial");
test_rel(z[4].re, c, 1e-9, "z4.real, 8th-order polynomial");
test_rel(z[4].im, s, 1e-9, "z4.imag, 8th-order polynomial");
test_rel(z[5].re, c, 1e-9, "z5.real, 8th-order polynomial");
test_rel(z[5].im, -s, 1e-9, "z5.imag, 8th-order polynomial");
test_rel(z[6].re, s, 1e-9, "z6.real, 8th-order polynomial");
test_rel(z[6].im, c, 1e-9, "z6.imag, 8th-order polynomial");
test_rel(z[7].re, s, 1e-9, "z7.real, 8th-order polynomial");
test_rel(z[7].im, -c, 1e-9, "z7.imag, 8th-order polynomial");

xa = GSL::Poly.alloc(0.16, 0.97, 1.94, 2.74, 3.58, 3.73, 4.70)
ya = GSL::Poly.alloc(0.73, 1.11, 1.49, 1.84, 2.30, 2.41, 3.07)
dd_expected = GSL::Vector.alloc(7.30000000000000e-01,
                               4.69135802469136e-01,
                              -4.34737219941284e-02,
                               2.68681098870099e-02,
                              -3.22937056934996e-03,
                               6.12763259971375e-03,
                              -6.45402453527083e-03)
dd = GSL::Poly.dd_init(xa, ya)

for i in 0...7
  GSL::Test::test_rel(dd[i], dd_expected[i], 1e-10, "divided difference dd[#{i}]")
end
#p dd.class

for i in 0...7
  y = dd.eval(xa, xa[i]);
  test_rel(y, ya[i], 1e-10, "divided difference y[#{i}]");
end

coeff = dd.taylor(1.5, xa)
#coeff = dd.taylor(1.5, 7, GSL::Vector.alloc(7))
#coeff = dd.taylor(1.5, 7)
#coeff = dd.taylor(1.5, GSL::Vector.alloc(7))

#p coeff.class

for i in 0...7
  y = coeff.eval(xa[i] - 1.5)
  test_rel(y, ya[i], 1e-10, "taylor expansion about 1.5 y[#{i}]");
end
