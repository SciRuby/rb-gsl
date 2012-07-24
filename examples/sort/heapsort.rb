#!/usr/bin/env ruby
require("gsl")

v = GSL::Vector.alloc(1, 2, 6, 7, 8, 9, 3, 4, 5)

p v.heapsort { |a, b|
  b <=> a
}
p v

p v.heapsort! { |a, b|
  b <=> a
}
p v

v = GSL::Vector.alloc(1, 2, 6, 7, 8, 9, 3, 4, 5)
p v.heapsort_index { |a, b|
  a <=> b
}

p GSL.heapsort(v) { |a, b|
  b <=> a
}
