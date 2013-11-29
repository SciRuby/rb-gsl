#!/usr/bin/env ruby
require("gsl")
pp = 0.5
nn = 3.5
IO.popen("graph -T X -C -g 3 -x 0 10 -y 0 0.3 -L 'Negative Binomial Distribution, p = 0.5, n = 3.5'", "w") do |io|
  for i in 0..10 do
    y = GSL::Ran::negative_binomial_pdf(i, pp, nn)
    io.printf("%d %e\n%d %e\n", i, y, i+1, y)
  end
end
