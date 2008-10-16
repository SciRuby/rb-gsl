#!/usr/bin/env ruby
require('gsl')

puts("Solve 2 - 3*x + x*x = 0")

p GSL::Poly.solve_quadratic([1, -3, 2])
p GSL::Poly.solve_quadratic(1, -3, 2)
z = GSL::Poly.complex_solve_quadratic(1, -3, 2)
printf("%f %f\n", z[0].re, z[0].im)
printf("%f %f\n", z[1].re, z[1].im)
#p GSL::Poly.complex_solve_quadratic([1, -3, 2])
#z = GSL::Poly.solve([2, -3, 1])

poly = GSL::Poly.alloc(2, -3, 1)
z = poly.solve
printf("%f %f\n", z[0].re, z[0].im)
printf("%f %f\n", z[1].re, z[1].im)

z = GSL::Poly.solve([2, -3, 1])
printf("%f %f\n", z[0].re, z[0].im)
printf("%f %f\n", z[1].re, z[1].im)



__END__
