#!/usr/bin/env ruby
require("gsl")
include GSL

x, y = Vector.filescan("points")

#spline = Spline.alloc(Interp::CSPLINE, n)
#spline.init(x, y)
#spline = Spline.alloc(x, y, "cspline")
#spline = Spline.alloc("cspline", x, y)
spline = Spline.alloc(x, y)

x2 = Vector.linspace(x[0], x[-1], 100)
y2 = spline.eval(x2)

graph([x, y], [x2, y2], "-C -g 3 -S 2")

#p spline.name
#p spline.min_size

