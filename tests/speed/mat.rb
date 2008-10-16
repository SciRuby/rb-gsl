#!/usr/bin/env ruby
require("gsl")
require("gslbench")
include GSL::Bench

REPEAT = 1

a = Matrix[0...40000, 200, 200]
b = Matrix[0...40000, 200, 200]

puts("calculating ...")
Bench::bench_time(REPEAT) { c = a*b }

