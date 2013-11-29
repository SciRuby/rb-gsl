#!/usr/bin/env ruby
require("gsl")

r1 = GSL::Rng.alloc(GSL::Rng::RNGEXTRA_RNG1)
p r1.name
p r1.uniform
p r1.gaussian

r2 = GSL::Rng.alloc("rngextra_rng2")
p r2.name
p r2.get
p r2.get
p r2.uniform
p r2.poisson(3)
