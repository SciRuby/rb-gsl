#!/usr/bin/env ruby
# Ruby/GSL implementation of GSL "qrng/test.c"
require("gsl")
require("./gsl_test2.rb")
include GSL::Test
include Math

GSL::IEEE::env_setup()

def test_sobol()
  v = GSL::Vector.alloc(3)

  status = 0
  g = GSL::QRng.alloc(GSL::QRng::SOBOL, 2)
  g.get(v)
  g.get(v)
  g.get(v)
  status += (v[0] != 0.25 or v[1] != 0.75) ? 1 : 0
  g.get(v)
  status += (v[0] != 0.375 or v[1] != 0.375) ? 1 : 0
  GSL::Test::test(status, "Sobol d=2")

  status = 0
  g = GSL::QRng.alloc(GSL::QRng::SOBOL, 3)
  g.get(v)
  g.get(v)
  g.get(v)
  status += (v[0] != 0.25 or v[1] != 0.75) ? 1 : 0
  g.get(v)
  status += (v[0] != 0.375 or v[1] != 0.375) ? 1 : 0
  GSL::Test::test(status, "Sobol d=3")

  status = 0
  g.init
  g.get(v)
  g.get(v)
  g.get(v)
  status += (v[0] != 0.25 or v[1] != 0.75) ? 1 : 0
  g.get(v)
  status += (v[0] != 0.375 or v[1] != 0.375) ? 1 : 0
  GSL::Test::test(status, "Sobol d=3 (reinitialized)")
end

def test_nied2()
  v = GSL::Vector.alloc(3)

  status = 0
  g = GSL::QRng.alloc(GSL::QRng::NIEDERREITER_2, 2)
  g.get(v)
  g.get(v)
  g.get(v)
  status += (v[0] != 0.75 or v[1] != 0.25) ? 1 : 0
  g.get(v)
  status += (v[0] != 0.25 or v[1] != 0.75) ? 1 : 0
  g.get(v)
  g.get(v)
  g.get(v)
  status += (v[0] != 0.625 or v[1] != 0.125) ? 1 : 0
  GSL::Test::test(status, "Niederreiter d=2")

  status = 0
  g = GSL::QRng.alloc(GSL::QRng::NIEDERREITER_2, 3)
  g.get(v)
  g.get(v)
  g.get(v)
  status += (v[0] != 0.75 or v[1] != 0.25) ? 1 : 0
  g.get(v)
  status += (v[0] != 0.25 or v[1] != 0.75) ? 1 : 0
  g.get(v)
  g.get(v)
  g.get(v)
  status += (v[0] != 0.625 or v[1] != 0.125) ? 1 : 0
  GSL::Test::test(status, "Niederreiter d=3")

  status = 0
  g.init
  g.get(v)
  g.get(v)
  g.get(v)
  status += (v[0] != 0.75 or v[1] != 0.25) ? 1 : 0
  g.get(v)
  status += (v[0] != 0.25 or v[1] != 0.75) ? 1 : 0
  g.get(v)
  g.get(v)
  g.get(v)
  status += (v[0] != 0.625 or v[1] != 0.125) ? 1 : 0
  GSL::Test::test(status, "Niederreiter d=3 (reinitialized)")
end

test_sobol()
test_nied2()

def test_hdsobol()
  v = GSL::Vector.alloc(3)

  status = 0
  g = GSL::QRng.alloc(GSL::QRng::HDSOBOL, 2)
  g.get(v)
  g.get(v)
  g.get(v)
  status += (v[0] != 0.25 or v[1] != 0.75) ? 1 : 0
  g.get(v)
  status += (v[0] != 0.375 or v[1] != 0.375) ? 1 : 0
  GSL::Test::test(status, "HDSobol d=2")

  status = 0
  g = GSL::QRng.alloc(GSL::QRng::SOBOL, 3)
  g.get(v)
  g.get(v)
  g.get(v)
  status += (v[0] != 0.25 or v[1] != 0.75 or v[2] != 0.25) ? 1 : 0
  g.get(v)
  status += (v[0] != 0.375 or v[1] != 0.375 or v[2] != 0.625) ? 1 : 0
  GSL::Test::test(status, "HDSobol d=3")

  status = 0
  g.init
  g.get(v)
  g.get(v)
  g.get(v)
  status += (v[0] != 0.25 or v[1] != 0.75 or v[2] != 0.25) ? 1 : 0
  g.get(v)
  status += (v[0] != 0.375 or v[1] != 0.375 or v[2] != 0.625) ? 1 : 0
  GSL::Test::test(status, "HDSobol d=3 (reinitialized)")
end

test_hdsobol()
