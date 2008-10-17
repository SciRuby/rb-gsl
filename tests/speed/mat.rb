#!/usr/bin/env ruby
require("gsl")
require("gslbench")
include GSL::Bench

REPEAT = 1

a = GSL::Matrix[0...40000, 200, 200]
b = GSL::Matrix[0...40000, 200, 200]

puts("calculating ...")
bench_time(REPEAT) { c = a*b }

