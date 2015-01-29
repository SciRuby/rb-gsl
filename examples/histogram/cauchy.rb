#!/usr/bin/env ruby
# Usage from command line:
#  % gsl-randist 0 10000 cauchy 30 | ./hist1d.rb -100 100 200

require("gsl")

if ARGV.size != 3
  puts("Usage: gsl-histogram xmin xmax n")
  puts("  Computes a histogram of the data")
  puts("  on stdin using n bins from xmin to xmax")
end

a = ARGV.shift.to_f
b = ARGV.shift.to_f
n = ARGV.shift.to_i

h = GSL::Histogram.alloc(n, a, b)

while line = STDIN.gets
  x = line.chomp.to_f
  h.increment(x)
end

h.graph("-C -g 3 -L 'gsl-randist 0 10000 cauchy 30'")



