#!/usr/bin/env ruby
require("gsl")
require("../gsl_test2.rb")
include GSL::Test
include Math

Longley_n = 16;
Longley_p = 7;


Longley_x = GSL::Vector.alloc(1,  83.0,   234289,   2356,     1590,    107608,  1947,
                       1,  88.5,   259426,   2325,     1456,    108632,  1948,
                       1,  88.2,   258054,   3682,     1616,    109773,  1949,
                       1,  89.5,   284599,   3351,     1650,    110929,  1950,
                       1,  96.2,   328975,   2099,     3099,    112075,  1951,
                       1,  98.1,   346999,   1932,     3594,    113270,  1952,
                       1,  99.0,   365385,   1870,     3547,    115094,  1953,
                       1, 100.0,   363112,   3578,     3350,    116219,  1954,
                       1, 101.2,   397469,   2904,     3048,    117388,  1955,
                       1, 104.6,   419180,   2822,     2857,    118734,  1956,
                       1, 108.4,   442769,   2936,     2798,    120445,  1957,
                       1, 110.8,   444546,   4681,     2637,    121950,  1958,
                       1, 112.6,   482704,   3813,     2552,    123366,  1959,
                       1, 114.2,   502601,   3931,     2514,    125368,  1960,
                       1, 115.7,   518173,   4806,     2572,    127852,  1961,
                       1, 116.9,   554894,   4007,     2827,    130081,  1962)

Longley_y = GSL::Vector.alloc(60323, 61122, 60171, 61187, 63221, 63639, 64989, 63761,
                       66019, 67857, 68169, 66513, 68655, 69564, 69331, 70551)

def test_longley()
  work = GSL::MultiFit::Workspace.alloc(Longley_n, Longley_p)
  x = GSL::Matrix.alloc(Longley_x, Longley_n, Longley_p).view
  y = Longley_y.view
  expected_c = GSL::Vector.alloc(-3482258.63459582,
                          15.0618722713733,
                          -0.358191792925910E-01,
                          -2.02022980381683,
                          -1.03322686717359,
                          -0.511041056535807E-01,
                          1829.15146461355)
  expected_sd = GSL::Vector.alloc(890420.383607373,      
                           84.9149257747669,      
                           0.334910077722432E-01, 
                           0.488399681651699,     
                           0.214274163161675,     
                           0.226073200069370,     
                           455.478499142212)
  expected_chisq = 836424.055505915

  c, cov, chisq, status = GSL::MultiFit.linear(x, y, work)
  for i in 0...7
    test_rel(c[i], expected_c[i], 1e-10, "longley gsl_fit_multilinear c#{i}") 
  end
  diag = cov.diagonal
  test_rel(diag[0], pow(expected_sd[0],2.0), 1e-10, "longley gsl_fit_multilinear cov00") ;
  test_rel(diag[1], pow(expected_sd[1],2.0), 1e-10, "longley gsl_fit_multilinear cov11") ;
  test_rel(diag[2], pow(expected_sd[2],2.0), 1e-10, "longley gsl_fit_multilinear cov22") ;
  test_rel(diag[3], pow(expected_sd[3],2.0), 1e-10, "longley gsl_fit_multilinear cov33") ;
  test_rel(diag[4], pow(expected_sd[4],2.0), 1e-10, "longley gsl_fit_multilinear cov44") ;
  test_rel(diag[5], pow(expected_sd[5],2.0), 1e-10, "longley gsl_fit_multilinear cov55") ;
  test_rel(diag[6], pow(expected_sd[6],2.0), 1e-10, "longley gsl_fit_multilinear cov66") ;
  test_rel(chisq, expected_chisq, 1e-10, "longley gsl_fit_multilinear chisq")

###################
  expected_cov = GSL::Matrix.alloc([8531122.56783558,
                             -166.727799925578, 0.261873708176346, 3.91188317230983,
                             1.1285582054705, -0.889550869422687, -4362.58709870581],
                            
                            [-166.727799925578, 0.0775861253030891, -1.98725210399982e-05,
                             -0.000247667096727256, -6.82911920718824e-05, 0.000136160797527761,
                             0.0775255245956248],
                            
                            [0.261873708176346, -1.98725210399982e-05, 1.20690316701888e-08,
                             1.66429546772984e-07, 3.61843600487847e-08, -6.78805814483582e-08,
                             -0.00013158719037715],
                            
                            [3.91188317230983, -0.000247667096727256, 1.66429546772984e-07,
                             2.56665052544717e-06, 6.96541409215597e-07, -9.00858307771567e-07,
                             -0.00197260370663974],
                            
                            [1.1285582054705, -6.82911920718824e-05, 3.61843600487847e-08,
                             6.96541409215597e-07, 4.94032602583969e-07, -9.8469143760973e-08,
                             -0.000576921112208274],
                            [-0.889550869422687, 0.000136160797527761, -6.78805814483582e-08,
-9.00858307771567e-07, -9.8469143760973e-08, 5.49938542664952e-07,
0.000430074434198215],

[-4362.58709870581, 0.0775255245956248, -0.00013158719037715,
-0.00197260370663974, -0.000576921112208274, 0.000430074434198215,
 2.23229587481535 ])

  expected_chisq = 836424.055505915
  w = GSL::Vector.alloc(Longley_n)
  w.set_all(1.0)
  c, cov, chisq, status = GSL::MultiFit.wlinear(x, w, y, work)
  for i in 0...7
    test_rel(c[i], expected_c[i], 1e-10, "longley gsl_fit_wmultilinear c#{i}") ;
  end

  for i in 0...Longley_p
    for j in 0...Longley_p
      test_rel(cov[i][j], expected_cov[i][j], 1e-7, 
               "longley gsl_fit_wmultilinear cov(#{i},#{j})")
    end
  end
  test_rel(chisq, expected_chisq, 1e-10, "longley gsl_fit_wmultilinear chisq")
end

test_longley()
