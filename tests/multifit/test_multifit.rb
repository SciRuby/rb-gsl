#!/usr/bin/env ruby
require("gsl")
require("../gsl_test2.rb")
include GSL::Test
include Math

def test_lmder(fdf, x, xx, f, cov)
  s = GSL::MultiFit::FdfSolver.alloc("lmsder", fdf.n, fdf.p)
  s.set(fdf, x)
  iter = 0
  begin
    _status = s.iterate
    for i in 0...fdf.p
      test_rel(s.x[i], xx[fdf.p*iter+i], 1e-5, "lmsder,  iter=#{iter}, x#{i}")
    end
    test_rel(Blas.dnrm2(s.f), f[iter], 1e-5, "lmsder, iter=#{iter}, f")
    iter += 1
  end while iter < 20
  covar = s.covar(0.0)
  for i in 0...fdf.p
    for j in 0...fdf.p
      test_rel(covar[i,j], cov[i*fdf.p+j], 1e-7, 
               "gsl_multifit_covar cov(#{i},#{j})")
    end
  end
end

def test_fdf(name, fdf, x, x_final, f_sumsq, sigma)
  s = GSL::MultiFit::FdfSolver.alloc("lmsder", fdf.n, fdf.p)
  s.set(fdf, x)
  iter = 0
  begin
    status = s.iterate
#    status = GSL::MultiFit::test_delta(s.dx, s.x, 0.0, 1e-7)
    status = s.test_delta(0.0, 1e-7)
    iter += 1
  end while status == GSL::CONTINUE and iter < 1000
  covar = s.covar(0.0)
  for i in 0...fdf.p
    test_rel(s.x[i], x_final[i], 1e-5, "#{name}, lmsder, x#{i}")
  end
  s2 = pow(Blas.dnrm2(s.f), 2.0)
  test_rel(s2, f_sumsq, 1e-5, "#{name}, lmsder, |f|^2")
  for i in 0...fdf.p
    ei = sqrt(s2/(fdf.n-fdf.p))*sqrt(covar[i,i]);
    test_rel(ei, sigma[i], 1e-4, "#{name}, sigma(#{i})")
  end
end

GSL::IEEE::env_setup()
  

