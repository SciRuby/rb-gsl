#!/usr/bin/env ruby
require('gsl')

puts("Solve x^3 - 1 == 0")
puts("x = 1, (-1 +/- i sqrt(3))/2")
p GSL::Poly.complex_solve_cubic(0, 0, -1)


__END__
