#!/usr/bin/env ruby
require("gsl")

p = GSL::Permutation::alloc(10)
p.init
p.swap(4, 7)
#p p.to_a
p.print

#p p.inverse.to_a
#p p.next.to_a
#p p.next.to_a
#p p.next.to_a

p.printf("%u \n")

