#!/usr/bin/env ruby
require("gsl")
include Math
include GSL

n = 10
xa = Vector[n]
ya = Vector[n]
for i in 0...n do
  a = i.to_f
  xa[i] = i + 0.5*sin(a)
  ya[i] = i + cos(a*a)
end

def spline_compare(type, xa, ya, filename)
  n = xa.size
  spline = Spline.alloc(type, xa, ya)
  p spline.name

  xi = xa[0]
  File.open(filename, "w") do |file|
    while xi < xa[n-1]
      yi = spline.eval(xi)
      file.printf("%e %e\n", xi, yi)
      xi += 0.01
    end
  end
end

types = ["linear", "polynomial", "cspline", "cspline_periodic",
         "akima", "akima_periodic"]

types.each do |t|
  filename = t + ".dat"
  spline_compare(t, xa, ya, filename)
end

#puts("6 interpolation types are examined.")
#puts("Try \"gnuplot -persist spline-compare.gp\"")
system("gnuplot -persist demo.gp")

types.each do |t|
  filename = t + ".dat"
  File.delete(filename)
end
