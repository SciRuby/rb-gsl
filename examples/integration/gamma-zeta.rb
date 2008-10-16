#!/usr/bin/env ruby
require("gsl")
include GSL
include Math

w = Integration::Workspace.alloc(1000)
xmin = 0.0

f = Function.alloc{ |x, n|
  pow(x, n-1)/(exp(x) - 1.0)
}

nv = Vector[2..5]
y1 = Sf::gamma(nv)*Sf::zeta(nv)

y2 = Vector[10]
begin
  file1 = File.open("y1.dat", "w")
  file2 = File.open("y2.dat", "w")
  nv.each do |n|
    f.set_params(n)
    y2[n.to_i-2] = f.qagiu(xmin, w)[0]
    file1.printf("%e %e\n", n, y1[n.to_i-2])
    file2.printf("%e %e\n", n, y2[n.to_i-2])
  end
ensure
  file1.close
  file2.close
end

system("graph -T X -C -g 3 -m 1 'y2.dat' -m -2 -S 4 'y1.dat' ")
File.delete("y1.dat")
File.delete("y2.dat")


