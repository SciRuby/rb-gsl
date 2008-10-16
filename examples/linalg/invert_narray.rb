#!/usr/bin/env ruby
require("gsl")
include GSL

na = NArray[[1.0, 4], [2, 3]]
lu, perm, signum = Linalg::LU.decomp(na)

inv = Linalg::LU.invert(lu, perm)
p inv

