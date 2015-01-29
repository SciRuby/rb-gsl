#!/usr/bin/env ruby
require("gsl")

r = GSL::Rng.alloc
10.times {
  p r.uniform
}

__END__
#include GSL
#require 'gsl/gsl_rng'

r = GSL::Random::Rng.alloc
r.set(2)
p r.get
p r.get
p r.max
p r.get
p r.get

r2 = GSL::Random::Rng.alloc
p r2.uniform
p r2.uniform

__END__
rf = GSL::Random::Rng.alloc
p rf.uniform
p rf.uniform

p rf.name


p gsl_rng_types_setup
