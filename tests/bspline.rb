#!/usr/bin/env ruby
require("gsl")
require("./gsl_test.rb")
include GSL::Test

def test_bspline(bw)
  n = 100
  ncoeffs = bw.ncoeffs
  order = bw.order
  nbreak = bw.nbreak
  a = bw.breakpoint(0)
  b = bw.breakpoint(nbreak - 1)
  
  for i in 0...n do
    xi = a + (b - a)*(i/(n-1))
    sum = 0.0
    bb = bw.eval(xi)
    for j in 0...ncoeffs do
      bj = bb[j]
      s = bj < 0 or bj > 1
      GSL::Test::test(s,  "basis-spline coefficient #{j} is in range [0,1] for x=#{xi}")
      sum += bj
    end
    GSL::Test::test_rel(sum, 1.0, order*GSL::DBL_EPSILON, "basis-spline coefficient #{order} is in range [0,1] for x=#{xi}")
  end
end

NMAX = 10
BMAX = 100
for order in 1..NMAX do
  for breakpoints in 2..BMAX do
    a = -1.23*order
    b = 45.6*order
    bw = GSL::BSpline.alloc(order, breakpoints)
    knots = bw.knots_uniform(a, b)
    test_bspline(bw)
  end
end

for order in 1..NMAX do
  for breakpoints in 2..BMAX do
    a = -1.23*order
    b = 45.6*order
    bw = GSL::BSpline.alloc(order, breakpoints)
    k = GSL::Vector.alloc(breakpoints)
    for i in 0...breakpoints do
      f = Math.sqrt(i.to_f/(breakpoints - 1.0))
      k[i] = (1 - f)*a + f*b
    end
    knots = bw.knots(k)    
    test_bspline(bw)
  end
end
