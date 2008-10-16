#!/usr/bin/env ruby
# This is equivalent to the gsl-histogram program

require("gsl")
require("getopts")

getopts(nil, "DISPLAY_STATS")

case ARGV.size
when 2
  a = ARGV[0].to_f
  b = ARGV[1].to_f
  n = (b - a).to_i
when 3
  a = ARGV[0].to_f
  b = ARGV[1].to_f
  n = ARGV[2].to_i
else
  puts("Usage: : gsl-histogram.rb [--DISPLAY_STATS] xmin xmax [n]")
  puts("Computes a histogram of the data on stdin using n bins from xmin to xmax.")
  puts("If n is unspecified then bins of integer width are used.")
  exit
end

h = GSL::Histogram.alloc(n)
h.set_ranges_uniform(a, b)

while line = STDIN.gets
  x = line.chomp.split[0].to_f
  h.increment(x)
end

if $OPT_DISPLAY_STATS
  printf("# mean = %g\n", h.mean)
  printf("# sigma = %g\n", h.sigma)
end

h.fprintf(STDOUT, "%g", "%g")

exit
