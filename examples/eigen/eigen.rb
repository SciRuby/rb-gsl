#!/usr/bin/env ruby
require("gsl")

m = GSL::Matrix[[1.0, 1/2.0, 1/3.0, 1/4.0], [1/2.0, 1/3.0, 1/4.0, 1/5.0],
           [1/3.0, 1/4.0, 1/5.0, 1/6.0], [1/4.0, 1/5.0, 1/6.0, 1/7.0]]

eigval = m.eigen_symm
p eigval

eigval, eigvec = m.eigen_symmv

p eigval == GSL::Eigen.symm(m)
val, vec = GSL::Eigen.symmv(m)
p vec

i = 0
vec.each_col do |v|
  a = (m*v)/v
  if a != val[i]
    puts("error")
  end
  i += 1
end

# Diagonalization
b = eigvec.inv*m*eigvec
b.clean!(1e-10)
p b
d = b.diagonal
p d == eigval

prod = eigval.prod
p GSL.equal?(prod, m.det, 1e-10)

