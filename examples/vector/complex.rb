#!/usr/bin/env ruby
require("gsl")

v = GSL::Vector::Complex.alloc(5)
p v

v[2] = [3, 4]
p v[2]

i = 0
v.each do |elm|
  elm.re += i
  i += 1
end

v.each do |elm|
  p elm
end

p v[3]

v.set_all([2, 4.7])

v2 = v.subvector(1, 3)
p v2
p v2.size
v2.each do |vv|
  p vv
end

p v.real
p v.real.class

v.each do |elm|
  p elm
end

a = v.to_a
a.each do |c|
  p c
end
