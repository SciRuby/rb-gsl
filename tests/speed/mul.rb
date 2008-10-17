#!/usr/bin/env ruby
require("gsl")
require("gslbench")
include GSL::Bench

REPEAT = 100
SIZE = 100000

a = GSL::Vector[0...SIZE]
b = GSL::Vector[0...SIZE]

#a = NArray.float(SIZE).indgen!
#b = NArray.float(SIZE).indgen!

puts("calculating ...")
bench_time(REPEAT) { c = a*b }
