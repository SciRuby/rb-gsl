#!/usr/bin/env ruby
# Ruby/GSL implementation of GSL "fit/test.c"
require("gsl")
require("./gsl_test.rb")
include GSL::Test
include Math

GSL::IEEE::env_setup()

norris_n = 36
norris_x = GSL::Vector.alloc(0.2, 337.4, 118.2, 884.6, 10.1, 226.5, 666.3, 996.3,
                      448.6, 777.0, 558.2, 0.4, 0.6, 775.5, 666.9, 338.0, 
                      447.5, 11.6, 556.0, 228.1, 995.8, 887.6, 120.2, 0.3, 
                      0.3, 556.8, 339.1, 887.2, 999.0, 779.0, 11.1, 118.3,
                      229.2, 669.1, 448.9, 0.5)
norris_y = GSL::Vector.alloc(0.1, 338.8, 118.1, 888.0, 9.2, 228.1, 668.5, 998.5,
                      449.1, 778.9, 559.2, 0.3, 0.1, 778.1, 668.8, 339.3, 
                      448.9, 10.8, 557.7, 228.3, 998.0, 888.8, 119.6, 0.3, 
                      0.6, 557.6, 339.3, 888.0, 998.5, 778.9, 10.2, 117.6,
                      228.9, 668.4, 449.2, 0.2)

noint1_n = 11
noint1_x = GSL::Vector.alloc(60, 61, 62, 63, 64, 65, 66, 67, 68, 69, 70)
noint1_y = GSL::Vector.alloc(130, 131, 132, 133, 134, 135, 136, 137, 138, 139, 140)

noint2_n = 3
noint2_x = GSL::Vector.alloc(4, 5, 6)
noint2_y = GSL::Vector.alloc(3, 4, 4)

x = GSL::Vector.alloc(1000)
y = GSL::Vector.alloc(1000)
w = GSL::Vector.alloc(1000)

xstride = 2
wstride = 3
ystride = 5

x.stride = xstride
w.stride = wstride
y.stride = ystride

for i in 0...norris_n
  x.set(i, norris_x[i])
  w.set(i, 1.0)
  y.set(i, norris_y[i])
end

expected_c0 = -0.262323073774029
expected_c1 =  1.00211681802045 
expected_cov00 = pow(0.232818234301152, 2.0)
expected_cov01 = -7.74327536339570e-05
expected_cov11 = pow(0.429796848199937E-03, 2.0)
expected_sumsq = 26.6173985294224

c0, c1, cov00, cov01, cov11, sumsq = GSL::Fit.linear(x, y, norris_n)
GSL::Test::test_rel(c0, expected_c0, 1e-10, "norris gsl_fit_linear c0")
GSL::Test::test_rel(c1, expected_c1, 1e-10, "norris gsl_fit_linear c1")
GSL::Test::test_rel(cov00, expected_cov00, 1e-10, "norris gsl_fit_linear cov00")
GSL::Test::test_rel(cov01, expected_cov01, 1e-10, "norris gsl_fit_linear cov01")
GSL::Test::test_rel(cov11, expected_cov11, 1e-10, "norris gsl_fit_linear cov11")
GSL::Test::test_rel(sumsq, expected_sumsq, 1e-10, "norris gsl_fit_linear sumsq")

expected_c0 = -0.262323073774029
expected_c1 =  1.00211681802045 
expected_cov00 = 6.92384428759429e-02 
expected_cov01 = -9.89095016390515e-05
expected_cov11 = 2.35960747164148e-07 
expected_sumsq = 26.6173985294224

c0, c1, cov00, cov01, cov11, sumsq = GSL::Fit.wlinear(x, w, y, norris_n)
GSL::Test::test_rel(c0, expected_c0, 1e-10, "norris gsl_fit_wlinear c0")
GSL::Test::test_rel(c1, expected_c1, 1e-10, "norris gsl_fit_wlinear c1")
GSL::Test::test_rel(cov00, expected_cov00, 1e-10, "norris gsl_fit_wlinear cov00")
GSL::Test::test_rel(cov01, expected_cov01, 1e-10, "norris gsl_fit_wlinear cov01")
GSL::Test::test_rel(cov11, expected_cov11, 1e-10, "norris gsl_fit_wlinear cov11")
GSL::Test::test_rel(sumsq, expected_sumsq, 1e-10, "norris gsl_fit_wlinear sumsq")

for i in 0...noint1_n
  x.set(i, noint1_x[i])
  w.set(i, 1.0)
  y.set(i, noint1_y[i])
end

expected_c1 = 2.07438016528926 
expected_cov11 = pow(0.165289256198347E-01, 2.0)  
expected_sumsq = 127.272727272727

c1, cov11, sumsq = GSL::Fit.mul(x, y, noint1_n)
GSL::Test::test_rel(c1, expected_c1, 1e-10, "noint1 gsl_fit_mul c1")
GSL::Test::test_rel(cov11, expected_cov11, 1e-10, "noint1 gsl_fit_mul cov11")
GSL::Test::test_rel(sumsq, expected_sumsq, 1e-10, "noint1 gsl_fit_mul sumsq")

expected_c1 = 2.07438016528926
expected_cov11 = 2.14661371686165e-05
expected_sumsq = 127.272727272727

c1, cov11, sumsq = GSL::Fit.wmul(x, w, y, noint1_n)
GSL::Test::test_rel(c1, expected_c1, 1e-10, "noint1 gsl_fit_wmul c1")
GSL::Test::test_rel(cov11, expected_cov11, 1e-10, "noint1 gsl_fit_wmul cov11")
GSL::Test::test_rel(sumsq, expected_sumsq, 1e-10, "noint1 gsl_fit_wmul sumsq")

for i in 0...noint2_n
  x.set(i, noint2_x[i])
  w.set(i, 1.0)
  y.set(i, noint2_y[i])
end

expected_c1 = 0.727272727272727 
expected_cov11 = pow(0.420827318078432E-01, 2.0)
expected_sumsq = 0.272727272727273

c1, cov11, sumsq = GSL::Fit.mul(x, y, noint2_n)
GSL::Test::test_rel(c1, expected_c1, 1e-10, "noint2 gsl_fit_mul c1")
GSL::Test::test_rel(cov11, expected_cov11, 1e-10, "noint2 gsl_fit_mul cov11")
GSL::Test::test_rel(sumsq, expected_sumsq, 1e-10, "noint2 gsl_fit_mul sumsq")

expected_c1 = 0.727272727272727 
expected_cov11 = 1.29870129870130e-02
expected_sumsq = 0.272727272727273

c1, cov11, sumsq = GSL::Fit.wmul(x, w, y, noint2_n)
GSL::Test::test_rel(c1, expected_c1, 1e-10, "noint2 gsl_fit_wmul c1")
GSL::Test::test_rel(cov11, expected_cov11, 1e-10, "noint2 gsl_fit_wmul cov11")
GSL::Test::test_rel(sumsq, expected_sumsq, 1e-10, "noint2 gsl_fit_wmul sumsq")
