#!/usr/bin/env ruby
require('gsl')

poly = GSL::Poly.alloc(-1, 0, 0, 0, 0, 1)

w = GSL::Poly::Complex::Workspace.alloc(6)
z = poly.solve(w)

for i in 0...5 do
  printf("z%d = %+.18f %+.18f\n", i, z[i].re, z[i].im)
end

__END__
