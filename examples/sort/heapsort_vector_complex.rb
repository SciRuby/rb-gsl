#!/usr/bin/env ruby
require("gsl")

v = GSL::Vector::Complex.alloc(5)
v[0] = [1, 2]
v[1] = [3, 4]
v[2] = [0, 1]
v[3] = [15, 3]
v[4] = [5, 7]

p v.heapsort { |a, b|
  a.abs <=> b.abs
}

p v.heapsort { |a, b|
  b.abs <=> a.abs
}

p GSL.heapsort(v) { |a, b|
  b.abs <=> a.abs
}
