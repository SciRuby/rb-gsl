#!/usr/bin/env ruby
require("gsl")
n1 = 5
n2 = 20
t = 3
IO.popen("graph -T X -C -g 3 -x 0 10 -y 0 0.7 -L 'Hypergeometric Distribution, n1=5, n2=20, t=3'", "w") do |io|
  for i in 0..10 do
    y = GSL::Ran::hypergeometric_pdf(i, n1, n2, t)
    io.printf("%d %e\n%d %e\n", i, y, i+1, y)
  end
end
