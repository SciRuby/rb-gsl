#!/usr/bin/env ruby
require("gsl")

m = GSL::Matrix.vandermonde([1, 2, 3, 4])
p m
eval, evec = m.eigen_nonsymmv

p eval.real
p evec.real

GSL::Eigen::nonsymmv_sort(eval, evec, GSL::Eigen::SORT_ABS_ASC)
p eval.real
p evec.real


=begin
Octave result:

octave:1> m = vander([1 2 3 4])
m =

   1   1   1   1
   8   4   2   1
  27   9   3   1
  64  16   4   1

octave:2> [v, d] = eig(m)
v =

   0.106363   0.137820  -0.129196   0.052254
   0.228569   0.036136   0.407586  -0.375568
   0.462688  -0.289537   0.376326   0.794157
   0.849919  -0.946503  -0.821925  -0.474902

d =

   15.48974    0.00000    0.00000    0.00000
    0.00000   -7.70629    0.00000    0.00000
    0.00000    0.00000    1.29422    0.00000
    0.00000    0.00000    0.00000   -0.07768
=end
