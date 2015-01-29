#!/usr/bin/env ruby
require("gsl")

N = 10
pp = GSL::Permutation.alloc(N)
GSL::Rng::env_setup()
r = GSL::Rng.alloc("gsl_rng_default")

puts("initial permutation:")
pp.init
pp.fprintf(STDOUT, " %u")
printf("\n")

puts(" random permutation:");
r.shuffle(pp)
pp.fprintf(STDOUT, " %u")
printf("\n");

puts("inverse permutation:");
q = pp.inverse
q.fprintf(STDOUT, " %u")
printf ("\n");
