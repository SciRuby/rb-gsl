#!/usr/bin/env ruby
require("gsl")

m = NMatrix[[1.0, 1/2.0, 1/3.0, 1/4.0], [1/2.0, 1/3.0, 1/4.0, 1/5.0],
            [1/3.0, 1/4.0, 1/5.0, 1/6.0], [1/4.0, 1/5.0, 1/6.0, 1/7.0]]

p GSL::Eigen.symm(m)
p GSL::Eigen.symmv(m)

