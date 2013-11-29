#!/usr/bin/env ruby
# Ruby/GSL implementation of GSL "sum/test.c"
require("gsl")
require("./gsl_test2.rb")
include GSL::Test
include Math

N = 50

def check_trunc(t, expected, desc)
  w = GSL::Sum::Levin_utrunc.alloc(N)
  sum_accel, _err = w.accel(t)
  desc2 = sprintf("trunc result, %s", desc)
  GSL::Test::test_rel(sum_accel, expected, 1e-8, desc2)
end

def check_full(t, expected, desc)
  w = GSL::Sum::Levin_u.alloc(N)
  sum_accel, err_est, = w.accel(t)
  desc2 = sprintf("full result, %s", desc)
  GSL::Test::test_rel(sum_accel, expected, 1e-8, desc2)

  sd_est = -log10(err_est/sum_accel.abs);
  sd_actual = -log10(GSL::DBL_EPSILON + ((sum_accel - expected)/expected).abs);
  desc2 = sprintf("full significant digits, %s (%g vs %g)", 
                  desc, sd_est, sd_actual)
  GSL::Test::test((sd_est > sd_actual + 1.0) ? 1 : 0, desc2)
end

GSL::IEEE::env_setup()

Zeta_2 = M_PI*M_PI/6.0
t = GSL::Vector.alloc(N)
for n in 0...N
  np1 = n + 1.0
  t[n] = 1.0/(np1*np1)
end
check_trunc(t, Zeta_2, "zeta(2)")
check_full(t, Zeta_2, "zeta(2)")

x = 10.0
y = exp(x)
t[0] = 1.0
for n in 1...N
  t[n] = t[n-1]*(x/n)
end
check_trunc(t, y, "exp(10)");
check_full(t, y, "exp(10)");

x = -10.0
y = exp(x)
t[0] = 1.0
for n in 1...N
  t[n] = t[n-1]*(x/n)
end
check_trunc(t, y, "exp(-10)");
check_full(t, y, "exp(-10)");

x = 0.5
y = -log(1-x)
t[0] = x
for n in 1...N
  t[n] = t[n-1]*(x*n)/(n + 1.0)
end
check_trunc(t, y, "-log(1/2)")
check_full(t, y, "-log(1/2)")

x = -1.0
y = -log(1-x)
t[0] = x
for n in 1...N
  t[n] = t[n-1]*(x*n)/(n + 1.0)
end
check_trunc(t, y, "-log(2)")
check_full(t, y, "-log(2)")

result = 0.192594048773
t[0] = 3.0 / (M_PI * M_PI)
for n in 1...N
  t[n] = -t[n - 1] * (4.0 * (n + 1.0) - 1.0) / (M_PI * M_PI)
end
check_trunc(t, result, "asymptotic series")
check_full(t, result, "asymptotic series")

result = 0.5772156649015328606065120900824;
t[0] = 1.0;
for n in 1...N
  t[n] = 1/(n+1.0) + log(n/(n+1.0))
end
check_trunc(t, result, "Euler's constant")
check_full(t, result, "Euler's constant")

result = 0.6048986434216305 
for n in 0...N
  t[n] = (n%2 == 1 ? -1 : 1) * 1.0 / sqrt(n + 1.0)
end
check_trunc(t, result, "eta(1/2)")
check_full(t, result, "eta(1/2)")
