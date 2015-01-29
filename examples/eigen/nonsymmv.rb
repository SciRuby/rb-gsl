#!/usr/bin/env ruby
require("gsl")

m = GSL::Matrix[[1, 2, 3], [4, 5, 0], [6, 0, 0]]
evec, eval = m.eigen_nonsymmv
p evec
p eval

m = GSL::Matrix.vandermonde([-1, -2, 3, 4])
w = GSL::Eigen::Nonsymmv.alloc(4)
eval, evec = GSL::Eigen::nonsymmv(m, w)
p eval
p evec

=begin
This can be compared with the corresponding output from GNU OCTAVE,

  octave> [v,d] = eig(vander([-1 -2 3 4]));
octave> diag(d)
ans =

  -6.4139 + 0.0000i
   5.5456 + 3.0854i
   5.5456 - 3.0854i
   2.3228 + 0.0000i

octave> v
v =

 Columns 1 through 3:

  -0.09988 + 0.00000i  -0.04350 - 0.00755i  -0.04350 + 0.00755i
  -0.11125 + 0.00000i   0.06399 - 0.14224i   0.06399 + 0.14224i
   0.29250 + 0.00000i  -0.51518 + 0.04142i  -0.51518 - 0.04142i
   0.94451 + 0.00000i  -0.84059 + 0.00000i  -0.84059 - 0.00000i

 Column 4:

  -0.14493 + 0.00000i
   0.35660 + 0.00000i
   0.91937 + 0.00000i
   0.08118 + 0.00000i
=end
