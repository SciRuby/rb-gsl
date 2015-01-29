#!/usr/bin/env ruby
require("gsl")
include GSL

names = ["default", "mt19937", "mt19937_1999", "mt19937_1998", "ranlxs0", "ranlxs1",
"ranlxs2", "ranlxd1", "ranlxd2", "ranlux", "ranlux389", "cmrg", "mrg", "taus", "taus2",
"gfsr4", "rand", "random_bsd", "random8_bsd", "random32_bsd", "random64_bsd",
"random128_bsd", "random256_bsd", "random_libc5", "random_glibc2", "rand48", "ran0",
"ran1", "ran2", "ran3", "ranf", "ranmar", "r250", "tt800", "vax", "transputer",
"randu", "minstd", "uni", "uni32", "slatec", "zuf", "borosh13", "coveyou",
"fishman18", "fishman20", "fishman2x", "knuthran2", "knuthran", "lecuyer21",
"waterman14"]

names.each do |name|
  r = Rng.alloc(name)
  printf("%s %s\n", name, r.name)
end

names.each do |name|
  name2 = "gsl_rng_" + name
  r = Rng.alloc(name)
  printf("%s %s\n", name, r.name)
end


r = Rng.alloc(Rng::KNUTHRAN)
p r.name
