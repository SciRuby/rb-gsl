#!/usr/bin/env ruby
# Ruby/GSL implementation of GSL "odeiv-initval/test.c"
require("gsl")
require("./gsl_test2.rb")
include GSL::Test
include Math

rhs_linear = Proc.new { |t, y, f|
  f[0] = 0.0
  f[1] = y[0]
  GSL::SUCCESS
}

jac_linear = Proc.new { |t, y, dfdy, dfdt|
  dfdy.set(0, 0, 0.0)
  dfdy.set(0, 1, 0.0)
  dfdy.set(1, 0, 1.0)
  dfdy.set(1, 1, 0.0)
  dfdt[0] = 0.0
  dfdt[1] = 0.0
  GSL::SUCCESS
}

Rhs_func_lin = GSL::Odeiv::System.alloc(rhs_linear, jac_linear, 2)

rhs_sin = Proc.new { |t, y, f|
  f[0] = -y[1]
  f[1] = y[0]
  GSL::SUCCESS
}

jac_sin = Proc.new { |t, y, dfdy, dfdt|
  dfdy.set(0, 0, 0.0)
  dfdy.set(0, 1, -1.0)
  dfdy.set(1, 0, 1.0)
  dfdy.set(1, 1, 0.0)
  dfdt[0] = 0.0
  dfdt[1] = 0.0
  GSL::SUCCESS
}

Rhs_func_sin = GSL::Odeiv::System.alloc(rhs_sin, jac_sin, 2)

rhs_exp = Proc.new { |t, y, f|
  f[0] = y[1]
  f[1] = y[0]
  GSL::SUCCESS
}

jac_exp = Proc.new { |t, y, dfdy, dfdt|
  dfdy.set(0, 0, 0.0)
  dfdy.set(0, 1, 1.0)
  dfdy.set(1, 0, 1.0)
  dfdy.set(1, 1, 0.0)
  dfdt[0] = 0.0
  dfdt[1] = 0.0
  GSL::SUCCESS
}

Rhs_func_exp = GSL::Odeiv::System.alloc(rhs_exp, jac_exp, 2)

rhs_stiff = Proc.new { |t, y, f|
  f[0] = 998.0 * y[0] + 1998.0 * y[1]
  f[1] = -999.0 * y[0] - 1999.0 * y[1]
  GSL::SUCCESS
}

jac_stiff = Proc.new { |t, y, dfdy, dfdt|
  dfdy.set(0, 0, 998.0);
  dfdy.set(0, 1, 1998.0);
  dfdy.set(1, 0, -999.0);
  dfdy.set(1, 1, -1999.0);
  dfdt[0] = 0.0;
  dfdt[1] = 0.0;
  GSL::SUCCESS;
}

Rhs_func_stiff = GSL::Odeiv::System.alloc(rhs_stiff, jac_stiff, 2)

def test_stepper_linear(type, hstart, base_prec)
  h = hstart*1.0
  delmax = 0.0
  count = 0
  stepper = GSL::Odeiv::Step.alloc(type, 2)
  y = GSL::Vector.alloc(1.0, 0.0)
  yerr = GSL::Vector.alloc(2)
  t = 0.0
  s = 0
  while t < 4.0
    _status = stepper.apply(t, h, y, yerr, Rhs_func_lin)
    del = ((y[1] - (t + h))/y[1]).abs
    delmax = MAX_DBL(del, delmax)
    if del > (count + 1.0)*base_prec
      printf("  LINEAR(%20.17g)  %20.17g  %20.17g  %8.4g\n", t + h, y[1], t + h, del);
      s += 1
    end
    count += 1
    t += h
  end
  GSL::Test::test(s, "#{stepper.name}, linear [0,4], max relative error = #{delmax}")
end

def test_stepper_sin(type, hstart, base_prec)
  h = hstart*1.0
  delmax = 0.0
  count = 0
  stepper = GSL::Odeiv::Step.alloc(type, 2)
  y = GSL::Vector.alloc(1.0, 0.0)
  yerr = GSL::Vector.alloc(2)
  t = 0.0
  s = 0
  while t < M_PI
    sin_th = sin(t + h)
    _status = stepper.apply(t, h, y, yerr, Rhs_func_sin)
    del = ((y[1] - sin_th)/sin_th).abs
    delmax = GSL::MAX_DBL(del, delmax)
    if t < 0.5*M_PI
      stat = (del > (count+1.0)*base_prec) ? 1 : 0
    elsif t < 0.7*M_PI
      stat = (del > 1.0e+04*base_prec) ? 1 : 0
    elsif t < 0.9*M_PI
      stat = (del > 1.0e+06*base_prec) ? 1 : 0
    else
      stat = (del > 1.0e+09*base_prec) ? 1 : 0
    end
    if stat != 0
      printf("  SIN(%22.18g)  %22.18g  %22.18g  %10.6g\n", t + h, y[1], sin_th, del);
    end
    s += stat
    count += 1
    t += h
  end
  if delmax > 1.0e+09*base_prec
    s += 1
    printf("  SIN(0 .. M_PI)  delmax = %g\n", delmax)  
  end
  GSL::Test::test(s, "#{stepper.name}, sine [0,pi], max relative error = #{delmax}")

  delmax = 0.0
  while t < 3*M_PI
    _status = stepper.apply(t, h, y, yerr, Rhs_func_sin)
    del = (y[1] - sin(t)).abs
    delmax = GSL::MAX_DBL(del, delmax)
    count += 1
    t += h
  end
  if del > count * 2.0 * base_prec
    s += 1
    printf("  SIN(%22.18g)  %22.18g  %22.18g  %10.6g\n", t + h, y[1], sin(t), del)
  end
  test(s, "#{stepper.name}, sin [pi,3*pi], max absolute error = #{delmax}")
end

def test_stepper_exp(type, hstart, base_prec)
  s = 0
  y = GSL::Vector.alloc(1.0, 1.0)
  yerr = GSL::Vector.alloc(2)
  delmax = 0.0
  count = 0
  h = hstart*1.0
  stepper = GSL::Odeiv::Step.alloc(type, 2)
  t = 0.0
  while t < 5.0
    ex = exp(t + h)
    _status = stepper.apply(t, h, y, yerr, Rhs_func_exp)
    del = ((y[1] - ex)/y[1]).abs
    delmax = GSL::MAX_DBL(del, delmax)
    if del > (count+1.0)*2.0*base_prec
      printf("  EXP(%20.17g)  %20.17g  %20.17g  %8.4g\n", t + h, y[1], ex, del);
      s += 1
    end
    count += 1
    t += h
  end
  GSL::Test::test(s, "#{stepper.name}, exponential [0,5], max relative error = #{delmax}")
end

def test_stepper_stiff(type, hstart, base_prec)
  h = hstart*1.0
  s = 0
  delmax = 0.0
  count = 0
  stepper = GSL::Odeiv::Step.alloc(type, 2)
  y = GSL::Vector.alloc(1.0, 0.0)
  yerr = GSL::Vector.alloc(2)
  t = 0.0
  while t < 5.0
    _status = stepper.apply(t, h, y, yerr, Rhs_func_stiff)
    if t > 0.04
      arg = t + h
      e1 = exp(-arg)
      e2 = exp(-1000.0*arg)
      u = 2.0*e1 - e2
      del = ((y[0] - u)/y[0]).abs
      delmax = GSL::MAX_DBL(del, delmax)
      if del > (count + 1.0)*100.0*base_prec
        printf("  STIFF(%20.17g)  %20.17g  %20.17g  %8.4g\n", arg, y[0], u, del)
        s += 1
      end
    end
    count += 1
    t += h
  end
  GSL::Test::test(s, "#{stepper.name}, stiff [0,5], max relative error = #{delmax}")
end

def test_stepper_err(type, hstart, base_prec)
  h = hstart*1.0
  y = GSL::Vector.alloc(1.0, 0.0)
  yerr = GSL::Vector.alloc(2)
  delmax = 0.0
  errmax = 0.0
  count = 0
  s = 0

  stepper = GSL::Odeiv::Step.alloc(type, 2)

  t = 0.0
  while t < M_PI
    y1_t = y[1]
    dy_exp = cos(t)*sin(h)-2*sin(t)*pow(sin(h/2),2.0);
    _status = stepper.apply(t, h, y, yerr, Rhs_func_sin)
    dy_t = y[1] - y1_t;
    del = (dy_t - dy_exp).abs
    if t > 0.1 and t < 0.2
      stat = (del > 10.0*yerr[1].abs + GSL::DBL_EPSILON*y[1].abs) ? 1 : 0
      if stat != 0
        delmax = del
        errmax = yerr[1].abs
        printf("SIN(%.18e) %.5e %.5e %e %e %e\n", t + h, y[1], dy_t, dy_exp, del, yerr[1]);

        s += stat;
        break;
      end
    end
    count += 1
    t += h
  end
  test(s, "#{stepper.name}, sine [0,pi], accurary of estimate error = #{delmax} vs #{errmax}")
end

def test_evolve_system_flat(step, sys, t0, t1, hstart, y, yfin, err_target, desc)
  s = 0
  h = hstart*1.0
  t = t0*1.0
  e = GSL::Odeiv::Evolve.alloc(sys.dimension)
  while t < t1
    t, h, _status = e.apply(nil, step, sys, t, t1, h, y)
  end
  frac = ((y[1] - yfin[1])/yfin[1]).abs + ((y[0] - yfin[0])/yfin[0]).abs
  if frac > 2.0*e.count*err_target
    printf("FLAT t = %.5e  y0 = %g y1= %g\n", t, y[0], y[1]);
    s += 1
  end
  GSL::Test::test(s, "#{step.name}, #{desc}, evolve, no control, max relative error = #{frac}")
end

def test_evolve_system(type, sys, t0, t1, hstart, y, yfin, err_target, desc)
  s = 0
  t = t0*1.0
  h = hstart*1.0
  step = GSL::Odeiv::Step.alloc(type, sys.dimension)
  c = GSL::Odeiv::Control.standard_alloc(0.0, err_target, 1.0, 1.0)
  e = GSL::Odeiv::Evolve.alloc(sys.dimension)
  while t < t1
    t, h, _status = e.apply(c, step, sys, t, t1, h, y)
  end
  frac = ((y[1] - yfin[1])/yfin[1]).abs + ((y[0] - yfin[0])/yfin[0]).abs
  if frac > 2.0*e.count*err_target
    printf("SYS t = %.5e h = %g y0 = %g y1= %g\n", t, h, y[0], y[1]);
    s += 1
  end
  GSL::Test::test(s, "#{step.name}, #{desc}, evolve, standard control, relative error = #{frac}")
end

def test_evolve_exp(type, hstart, err)
  h = hstart*1.0
  y = GSL::Vector.alloc(1.0, 1.0)
  eee = exp(5.0)
  yfin = GSL::Vector.alloc(eee, eee)
  test_evolve_system(type, Rhs_func_exp, 0.0, 5.0, h, y, yfin, err, "exp [0,5]")
end

def test_evolve_sin(type, hstart, err)
  h = hstart*1.0
  y = GSL::Vector.alloc(1.0, 0.0)
  yfin = GSL::Vector.alloc(cos(2.0), sin(2.0))
  test_evolve_system(type, Rhs_func_sin, 0.0, 2.0, h, y, yfin, err, "sin [0,2]")
end

def test_evolve_stiff1(type, hstart, err)
  h = hstart*1.0
  y = GSL::Vector.alloc(1.0, 0.0)
  arg = 1.0
  e1 = exp(-arg)
  e2 = exp(-1000.0*arg)
  yfin = GSL::Vector.alloc(2.0*e1 - e2, -e1 + e2)
  test_evolve_system(type, Rhs_func_stiff, 0.0, 1.0, h, y, yfin, err, "stiff [0,1]")
end

def test_evolve_stiff5(type, hstart, err)
  h = hstart*1.0
  y = GSL::Vector.alloc(1.0, 0.0)
  arg = 5.0
  e1 = exp(-arg)
  e2 = exp(-1000.0*arg)
  yfin = GSL::Vector.alloc(2.0*e1 - e2, -e1 + e2)
less  test_evolve_system(type, Rhs_func_stiff, 0.0, 5.0, h, y, yfin, err, "stiff [0,5]")
end

ptypes = [
  {"type" => "rk2",
   "h" => 1e-3},
  {"type" => "rk2imp",
   "h" => 1e-3},
  {"type" => "rk4",
   "h" => 1e-3},
  {"type" => "rk4imp",
   "h" => 1e-3},
  {"type" => "rkf45",
   "h" => 1e-3},
  {"type" => "rk8pd",
   "h" => 1e-3},
  {"type" => "rkck",
   "h" => 1e-3},
  {"type" => "bsimp",
   "h" => 1e-3},
  {"type" => "gear1",
   "h" => 1e-3},
#  {"type" => "gear2",
#   "h" => 1e-3}
]

GSL::IEEE::env_setup()

ptypes.each do |hash|
  test_stepper_err(hash["type"], hash["h"], GSL::SQRT_DBL_EPSILON)
end

ptypes.each do |hash|
#  test_stepper_linear(hash["type"], hash["h"], GSL::SQRT_DBL_EPSILON)
#  test_stepper_sin(hash["type"], hash["h"], GSL::SQRT_DBL_EPSILON)
  test_stepper_sin(hash["type"], hash["h"]/10, 1e-8)  
#  test_stepper_exp(hash["type"], hash["h"], GSL::SQRT_DBL_EPSILON)
#  test_stepper_stiff(hash["type"], hash["h"], GSL::SQRT_DBL_EPSILON)
end

ptypes.each do |hash|
#  test_evolve_exp(hash["type"], hash["h"], GSL::SQRT_DBL_EPSILON)
  test_evolve_sin(hash["type"], hash["h"], GSL::SQRT_DBL_EPSILON)
#  test_evolve_stiff1(hash["type"], hash["h"], GSL::SQRT_DBL_EPSILON)
#  test_evolve_stiff5(hash["type"], hash["h"], GSL::SQRT_DBL_EPSILON)
end
