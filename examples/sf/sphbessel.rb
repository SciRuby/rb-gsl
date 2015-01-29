#!/usr/bin/env ruby
require("gsl")
include GSL::Sf

l = 2
x = 0.01
File.open("sphbessel.dat", "w") do |file|
  while x <= 12.0
    j0 = bessel_j0(x)
    j1 = bessel_j1(x)
    j2 = bessel_jl(2, x)
    j3 = bessel_jl(3,x)
    y0 = bessel_y0(x)
    y1 = bessel_y1(x)
    y2 = bessel_y2(x)
    i0 = bessel_i0_scaled(x)
    i1 = bessel_i1_scaled(x)
    i2 = bessel_i2_scaled(x)
    k0 = bessel_k0_scaled(x)
    k1 = bessel_k1_scaled(x)
    k2 = bessel_k2_scaled(x)
    file.printf("%e %e %e %e %e %e %e %e %e %e %e %e %e %e\n",
                x, j0, j1, j2, j3, y0, y1, y2, i0, i1, i2, k0, k1, k2)
    x += 0.1
  end
end

#puts("sphbessel.dat created. Try sphbessel.gp.")
system("gnuplot -persist sphbessel.gp")
File.delete("sphbessel.dat")
