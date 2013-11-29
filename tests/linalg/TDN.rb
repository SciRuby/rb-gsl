#!/usr/bin/env ruby
require("gsl")
require("../gsl_test2.rb")
require("./linalg.rb")
include GSL::Test

def test_TDN_solve_dim(dim, d, a, b, actual, eps)
  s = 0
  
  abovediag = GSL::Vector.alloc(dim-1)
  belowdiag = GSL::Vector.alloc(dim-1)	
  
  diag = GSL::Vector.alloc(dim)
  diag.set_all(d)
  rhs = GSL::Vector.indgen(dim) + 1
  
  abovediag.set_all(a)
  belowdiag.set_all(b)
  
  x = GSL::Linalg::solve_tridiag(diag, abovediag, belowdiag, rhs)
  p x
  p actual
  for i in 0...dim do
    si = x[i]
    ai = actual[i]
    foo = check(si, ai, eps);
    if foo > 0
      printf("%3lu[%lu]: %22.18g   %22.18g\n", dim, i, x[i], actual[i]);
      s += foo;
    end
  end
  return s
end

def test_TDN_solve()
  s = 0
  
  actual = GSL::Vector.alloc(5)
  actual[0] =  -7.0/3.0;
  actual[1] =  5.0/3.0;
  actual[2] =  4.0/3.0;
  f = test_TDN_solve_dim(3, 1.0, 2.0, 1.0, actual, 2.0 * GSL::DBL_EPSILON);
  GSL::test(f, "  solve_TDN dim=2 A");
  s += f;

  actual[0] =  0.75;
  actual[1] =  0.75;
  actual[2] =  2.625;
  f = test_TDN_solve_dim(3, 1.0, 1.0/3.0, 1.0/2.0, actual, 2.0 * GSL::DBL_EPSILON);
  GSL::test(f, "  solve_TDN dim=2 B");
  s += f;

  actual[0] =  99.0/140.0;
  actual[1] =  41.0/35.0;
  actual[2] =  19.0/10.0;
  actual[3] =  72.0/35.0;
  actual[4] =  139.0/35.0;
  f = test_TDN_solve_dim(5, 1.0, 1.0/4.0, 1.0/2.0, actual, 35.0/8.0 * GSL::DBL_EPSILON);
  GSL::test(f, "  solve_TDN dim=5");
  s += f;

  return s;
end

def test_TDN_cyc_solve_dim(dim, d, a, b, actual, eps)
  s = 0
  
  abovediag = GSL::Vector.alloc(dim)
  belowdiag = GSL::Vector.alloc(dim)
  diag = GSL::Vector.alloc(dim)
  rhs = GSL::Vector.indgen(dim) + 1
  
  abovediag.set_all(a)
  belowdiag.set_all(b)
  diag.set_all(d)
  
  x = GSL::Linalg::solve_cyc_tridiag(diag, abovediag, belowdiag, rhs);
  p x
  p actual
  for i in 0...dim do
    si = x[i]
    ai = actual[i]
    foo = check(si, ai, eps);
    if foo > 0
      printf("%3lu[%lu]: %22.18g   %22.18g\n", dim, i, x[i], actual[i]);
    end
    s += foo;
  end
  return s
end

def test_TDN_cyc_solve()
  s = 0
  actual = GSL::Vector.alloc(5)
  actual[0] =  3.0/2.0;
  actual[1] = -1.0/2.0;
  actual[2] =  1.0/2.0;
  f = test_TDN_cyc_solve_dim(3, 1.0, 2.0, 1.0, actual, 32.0 * GSL::DBL_EPSILON);
  GSL::test(f, "  solve_TDN_cyc dim=2 A");
  s += f;
  
  actual[0] = -5.0/22.0;
  actual[1] = -3.0/22.0;
  actual[2] =  29.0/22.0;
  actual[3] = -9.0/22.0;
  actual[4] =  43.0/22.0;
  f = test_TDN_cyc_solve_dim(5, 3.0, 2.0, 1.0, actual, 66.0 * GSL::DBL_EPSILON);
  GSL::test(f, "  solve_TDN_cyc dim=5");
  s += f;
  
  return s;
end

test_TDN_solve()
test_TDN_cyc_solve()

