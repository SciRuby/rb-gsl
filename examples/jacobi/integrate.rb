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

f = GSL::Vector.alloc(4*Q)

quad = Jac::Quadrature.alloc(Q)
quad.zwd(Jac::GJ, 0.0, 0.0)

for i in 0...Q do
  f[i] = fun(quad.x[i])
end

integr = quad.integrate(f)
exact = 2.0/3.0*Math::sin(3.0)

printf("Integral of cos(3x) from -1 to 1: %f\n", exact);
printf("jac_integrate result: %f\n", integr)
printf("Error: %e\n", integr - exact);

