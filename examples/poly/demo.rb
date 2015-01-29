#!/usr/bin/env ruby
require("gsl")
include GSL

# Polynomial p(x) = 1.5-1.25x-3.75x^2+x^4
poly = Poly[1.5, -1.25, -3.75, 0, 1]
# Solve the equation p(x) == 0
root = poly.solve    # Vector::Complex
# Extract only the real parts
# (imaginary parts are zero for this case)
re = root.real       # Vector::View

puts("p(x) = 1.5-1.25x-3.75x^2+x^4 == 0")
puts("Roots are found at #{re[0]}, #{re[1]}, #{re[2]}, #{re[3]}")

# Display the result
x = Vector.linspace(-2.5, 2.5, 20)
y = poly.eval(x)
zero = Vector.calloc(4)
graph([x, y], [re, zero], "-T X -C -g 3 -S 4 -X x -L 'p(x) = 1.5-1.25x-3.75x^2+x^4'")
