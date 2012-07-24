#!/usr/bin/env ruby
# Ruby/GSL implementation of GSL "cheb/test.c"
require("gsl")
require("./gsl_test.rb")
include GSL::Test
include Math

f_T0 = GSL::Function.alloc { |x| 1.0 }
f_T1 = GSL::Function.alloc { |x| x }
f_T2 = GSL::Function.alloc { |x| 2.0*x*x - 1.0 }
f_sin = GSL::Function.alloc { |x| sin(x) }

tol = 100.0*GSL::DBL_EPSILON

order = 40

cs = GSL::Cheb.alloc(order)
csd = GSL::Cheb.alloc(order)
csi = GSL::Cheb.alloc(order)

GSL::ieee_env_setup()

cs.init(f_T0, -1.0, 1.0)
for i in 0...order do
  c_exp = i == 0 ? 2.0 : 0.0
  desc = sprintf("c[%d] for T_0(x)", i)
  GSL::Test.test_abs(cs.c[i], c_exp, tol, desc)
end

cs.init(f_T1, -1.0, 1.0)
for i in 0...order do
  c_exp = i == 1 ? 1.0 : 0.0
  desc = sprintf("c[%d] for T_1(x)", i)
  GSL::Test.test_abs(cs.c[i], c_exp, tol, desc)
end

cs.init(f_T2, -1.0, 1.0)
for i in 0...order do
  c_exp = i == 2 ? 1.0 : 0.0
  desc = sprintf("c[%d] for T_2(x)", i)
  GSL::Test.test_abs(cs.c[i], c_exp, tol, desc)
end

cs.init(f_sin, -M_PI, M_PI)
GSL::Test.test_abs(cs.c[0], 0.0, tol, "c[0] for F_sin(x)")
GSL::Test.test_abs(cs.c[1], 5.69230686359506e-01, tol, "c[1] for F_sin(x)")
GSL::Test.test_abs(cs.c[2], 0.0, tol, "c[2] for F_sin(x)")
GSL::Test.test_abs(cs.c[3], -6.66916672405979e-01, tol, "c[3] for F_sin(x)")
GSL::Test.test_abs(cs.c[4], 0.0, tol, "c[4] for F_sin(x)")
GSL::Test.test_abs(cs.c[5], 1.04282368734237e-01, tol, "c[5] for F_sin(x)")

x = -M_PI
while x < M_PI
  r = cs.eval(x)
  desc = sprintf("GSL::Cheb#eval, sin\(%.3g\)", x)
  GSL::Test.test_abs(r, sin(x), tol, desc)
  x += M_PI/100.0
end

x = -M_PI
while x < M_PI
  r, e = cs.eval_err(x)
  desc = sprintf("GSL::Cheb#eval_err, sin\(%.3g\)", x)
  GSL::Test.test_abs(r, sin(x), tol, desc)
  desc = sprintf("GSL::Cheb#eval_err, error sin\(%.3g\)", x)
  GSL::Test.test_factor((r-sin(x)).abs + GSL::DBL_EPSILON, e, 10.0, desc)
  x += M_PI/100.0
end

x = -M_PI
while x < M_PI
  r = cs.eval_n(25, x)
  desc = sprintf("GSL::Cheb#eval_n, sin\(%.3g\)", x)
  GSL::Test.test_abs(r, sin(x), tol, desc)
  x += M_PI/100.0
end

x = -M_PI
while x < M_PI
  r, e = cs.eval_n_err(25, x)
  desc = sprintf("GSL::Cheb#eval_n_err, sin\(%.3g\)", x)
  GSL::Test.test_abs(r, sin(x), tol, desc)
  desc = sprintf("GSL::Cheb#eval_n_err, error sin\(%.3g\)", x)
  GSL::Test.test_factor((r-sin(x)).abs + GSL::DBL_EPSILON, e, 10.0, desc)
  x += M_PI/100.0
end

csd = cs.calc_deriv
#cs.calc_deriv(csd)
#csd = GSL::Cheb.calc_deriv(cs)
#GSL::Cheb.calc_deriv(csd, cs)
x = -M_PI
while x < M_PI
  r = csd.eval(x)
  desc = sprintf("GSL::Cheb#eval, deriv sin\(%.3g\)", x)
  GSL::Test.test_abs(r, cos(x), 1600*tol, desc)
  x += M_PI/100.0
end

csi = cs.calc_integ
#cs.calc_integ(csi)
#csi = GSL::Cheb.calc_integ(cs)
#GSL::Cheb.calc_integ(csi, cs)
x = -M_PI
while x < M_PI
  r = csi.eval(x)
  desc = sprintf("GSL::Cheb#eval, integ sin\(%.3g\)", x)
  GSL::Test.test_abs(r, -(1+cos(x)), tol, desc)
  x += M_PI/100.0
end


