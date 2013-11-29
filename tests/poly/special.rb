#!/usr/bin/env ruby
require("gsl")
require("../gsl_test2.rb")
include GSL::Test
include Math

_eps = 100.0 * GSL::DBL_EPSILON
GSL::IEEE::env_setup()

HermitPoly = Array[7]
HermitPoly[0] = GSL::Poly::Int[1]
HermitPoly[0][0] = 1
HermitPoly[1] = GSL::Poly::Int[0, 2]
HermitPoly[2] = GSL::Poly::Int[-2, 0, 4]
HermitPoly[3] = GSL::Poly::Int[0, -12, 0, 8]
HermitPoly[4] = GSL::Poly::Int[12, 0, -48, 0, 16]
HermitPoly[5] = GSL::Poly::Int[0, 120, 0, -160, 0, 32]
HermitPoly[6] = GSL::Poly::Int[-120, 0, 720, 0, -480, 0, 64]

for n in 0...6 do
  hn = GSL::Poly.hermite(n)
  GSL::Test::test2(hn == HermitPoly[n], "Hermite polynomial, n = #{n}")
end

LaguerrePoly = Array[7]
LaguerrePoly[0] = GSL::Poly::Int[1]
LaguerrePoly[0][0] = 1
LaguerrePoly[1] = GSL::Poly::Int[1, -1]
LaguerrePoly[2] = GSL::Poly::Int[2, -4, 1]
LaguerrePoly[3] = GSL::Poly::Int[6, -18, 9, -1]
LaguerrePoly[4] = GSL::Poly::Int[24, -96, 72, -16, 1]
LaguerrePoly[5] = GSL::Poly::Int[120, -600, 600, -200, 25, -1]
LaguerrePoly[6] = GSL::Poly::Int[720, -4320, 5400, -2400, 450, -36, 1]
for n in 0...7 do
  hn = GSL::Poly.laguerre(n)
  GSL::Test::test2(hn == LaguerrePoly[n], "Laguerre polynomial, n = #{n}")
end

ChebPoly = Array[7]
ChebPoly[0] = GSL::Poly::Int[1]
ChebPoly[0][0] = 1
ChebPoly[1] = GSL::Poly::Int[0, 1]
ChebPoly[2] = GSL::Poly::Int[-1, 0, 2]
ChebPoly[3] = GSL::Poly::Int[0, -3, 0, 4]
ChebPoly[4] = GSL::Poly::Int[1, 0, -8, 0, 8]
ChebPoly[5] = GSL::Poly::Int[0, 5, 0, -20, 0, 16]
ChebPoly[6] = GSL::Poly::Int[-1, 0, 18, 0, -48, 0, 32]
for n in 0...7 do
  hn = GSL::Poly.cheb(n)
  GSL::Test::test2(hn == ChebPoly[n], "Chebyshev polynomial, n = #{n}")
end

Cheb_IIPoly = Array[7]
Cheb_IIPoly[0] = GSL::Poly::Int[1]
Cheb_IIPoly[0][0] = 1
Cheb_IIPoly[1] = GSL::Poly::Int[0, 2]
Cheb_IIPoly[2] = GSL::Poly::Int[-1, 0, 4]
Cheb_IIPoly[3] = GSL::Poly::Int[0, -4, 0, 8]
Cheb_IIPoly[4] = GSL::Poly::Int[1, 0, -12, 0, 16]
Cheb_IIPoly[5] = GSL::Poly::Int[0, 6, 0, -32, 0, 32]
Cheb_IIPoly[6] = GSL::Poly::Int[-1, 0, 24, 0, -80, 0, 64]
for n in 0...7 do
  hn = GSL::Poly.cheb_II(n)
  GSL::Test::test2(hn == Cheb_IIPoly[n], "Chebyshev II polynomial, n = #{n}")
end
