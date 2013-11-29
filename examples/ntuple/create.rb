#!/usr/bin/env ruby
require("gsl")

v = GSL::Vector.alloc(3)
n = GSL::Ntuple.alloc("test.dat", v)

GSL::Rng.env_setup()
r = GSL::Rng.alloc
for i in 0...10000 do
  for j in 0...3 do
    v[j] = r.gaussian()
  end
  n.write
end

puts("Data file test.dat is created.")
puts("Try project.rb.")
