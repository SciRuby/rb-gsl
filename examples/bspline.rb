#!/usr/bin/env ruby
require("gsl")

N = 200
NCOEFFS = 8
NBREAK = NCOEFFS - 2

GSL::Rng::env_setup()
r = GSL::Rng.alloc()

bw = GSL::BSpline.alloc(4, NBREAK)
B = GSL::Vector.alloc(NCOEFFS)
x = GSL::Vector.alloc(N)
y = GSL::Vector.alloc(N)
xx = GSL::Matrix.alloc(N, NCOEFFS)
w = GSL::Vector.alloc(N)

#printf("#m=0,S=0\n")
for i in 0...N do
  xi = (15.0/(N-1)/1)*i
  yi = Math::cos(xi)*Math::exp(-0.1*xi)

  sigma = 0.1
  dy = GSL::Ran.gaussian(r, sigma)
  yi += dy

  x[i] = xi
  y[i] = yi
  w[i] = sigma

#  printf("%f %f\n", xi, yi)
end

bw.knots_uniform(0.0, 15.0)

for i in 0...N do
  xi = x[i]
  bw.eval(xi, B)
  for j in 0...NCOEFFS do
    xx[i,j] = B[j]
  end
end

c, cov, chisq = GSL::MultiFit.wlinear(xx, w, y)

#printf("#m=0,S=0\n")

x2 = GSL::Vector.linspace(0, 15, 150)
y2 = GSL::Vector.alloc(150)
x2.each_index do |i|
  bw.eval(x2[i], B)
  yi, yerr = GSL::MultiFit::linear_est(B, c, cov)
  y2[i] = yi
#  printf("%f %f\n", xi, yi)
end

GSL::graph([x, y], [x2, y2], "-T X -C -X x -Y y -x 0 15 -y -1 1.3")
