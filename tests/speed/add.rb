#!/usr/bin/env ruby
require("gsl")
require("gslbench")
include GSL::Bench

REPEAT = 100
SIZE = 100000

a = Vector[0...SIZE]
b = Vector[0...SIZE]

#a = NArray.float(SIZE).indgen!
#b = NArray.float(SIZE).indgen!

puts("calculating ...")
Bench::bench_time(REPEAT) { c = a + b }
