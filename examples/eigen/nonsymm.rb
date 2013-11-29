#!/usr/bin/env ruby
require("gsl")

m = GSL::Matrix[[1, 2], [3, 4]]

w = GSL::Eigen::Nonsymm::Workspace.alloc(2)
v = GSL::Vector::Complex.alloc(2)

#w.params(0, 1)

#p GSL::Eigen.nonsymm(m)
#p GSL::Eigen.nonsymm(m, v, w)

p m.eigen_nonsymm
#p m.eigen_nonsymm(v, w)

#p m.eigen_nonsymm_Z
#p m


=begin
Octave result:

octave:1> m = [1 2; 3 4]
m =

  1  2
  3  4

octave:2> eig(m)
ans =

  -0.37228
   5.37228


=end
