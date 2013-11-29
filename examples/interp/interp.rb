#!/usr/bin/env ruby
require("gsl")
include GSL

n = 10
x = Vector[n]
y = Vector[n]
#x = NArray.float(n)
#y = NArray.float(n)
File.open("data0.dat", "w") do |f|
  for i in 0...n do
    a = i.to_f
    x[i] = i + 0.5*Math::sin(a)
    y[i] = i + Math::cos(a*a)
    f.printf("%e %e\n", x[i], y[i])
  end
end

interp = Interp.alloc("akima", n)
interp.init(x, y)
#interp = Interp.alloc(x, y)
p interp.class
p interp.name

File.open("data1.dat", "w") do |f|
  xi = x[0]
  while xi < x[9]
    yi = interp.eval(x, y, xi)
    f.printf("%e %e\n", xi, yi)
    xi += 0.01
  end
end

system("graph -T X -g 3 -C -m -1 -S 4 data0.dat -S 1 data1.dat")
File.delete("data0.dat")
File.delete("data1.dat")

