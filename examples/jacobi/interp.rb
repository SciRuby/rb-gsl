#!/usr/bin/env ruby
require("gsl")

if ARGV.size != 1
  puts("Usage: integrate n\nn is the number of quadrature points.")
  exit
end

Q = ARGV[0].to_i
if Q < 1
  puts("Usage: integrate n\nn is the number of quadrature points.")
  exit
end

def fun(x)
  Math::cos(3.0*x)
end

def dfun(x)
  -3.0*Math::sin(3.0*x)
end

f = GSL::Vector.alloc(4*Q)
xp = GSL::Vector.linspace(0.1, 0.9, 9)

quad = Jac::Quadrature.alloc(Q)
quad.zwd(Jac::GJ, 0.0, 0.0)
quad.interpmat_alloc(xp)

for i in 0...Q do
  f[i] = fun(quad.x[i])
end

fout = quad.interpolate(f)

printf("X \t f(x) \t Error\n");

for i in 0...xp.size do
  printf("%f \t %f \t %e\n", xp[i], fout[i], fout[i] - fun(xp[i]));
end

quad.interpmat_free

