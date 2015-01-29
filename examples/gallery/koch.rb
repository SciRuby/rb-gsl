#!/usr/bin/env ruby
# Koch curve

require("gsl")
include Math

ONE_UNIT = 3

def koch(x, y, theta, size, order, file)
  if order == 0
    x += cos(theta)*size
    y += sin(theta)*size
    file.printf("%e %e\n", x, y)
  else
    x, y = koch(x, y, theta, size/3, order-1, file)
    theta += Math::PI/3
    x, y = koch(x, y, theta, size/3, order-1, file)
    theta -= 2.0*Math::PI/3
    x, y = koch(x, y, theta, size/3, order-1, file)
    theta += Math::PI/3
    x, y = koch(x, y, theta, size/3, order-1, file)
  end
  return [x, y]
end

SIZE = 243
ORDER = 4

x = 0.0
y = 0.0
theta = 0.0
IO.popen("graph -T X -C -N x -N y", "w") do |io|
  io.printf("%e %e\n", x, y)
  x, y = koch(x, y, theta, SIZE, ORDER, io)
  theta -= 2.0*Math::PI/3
  x, y = koch(x, y, theta, SIZE, ORDER, io)
  theta -= 2.0*Math::PI/3
  x, y = koch(x, y, theta, SIZE, ORDER, io)
  theta -= 2.0*Math::PI/3
end
