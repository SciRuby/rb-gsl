#!/usr/bin/env ruby
require("gsl")
include GSL

dim = 2
#q = QRng.alloc(QRng::SOBOL, dim)
q = QRng.alloc("sobol", dim)
#q = QRng.alloc("niederreiter_2", dim)
#q = QRng.alloc(QRng::NIEDERREITER_2, dim)

v = Vector.alloc(dim)
IO.popen("graph -T X -C --title-font-size 0.04 -L 'Distribution of first 1024 points from the quasi-random Sobol sequence' -m -1 -S 2", "w") do |io|
  for i in 0..1024 do
    #       v = q.get()    # by creating a alloc vector
    q.get(v)       # by using an existing vector (efficient)
    io.printf("%e %e\n", v[0], v[1])
  end
end

