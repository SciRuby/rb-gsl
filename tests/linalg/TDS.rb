#!/usr/bin/env ruby
require("gsl")
require("../gsl_test2.rb")
require("./linalg.rb")
include GSL::Test

def test_TDS_solve_dim(dim, d, od, actual, eps)
	s = 0

	diag = GSL::Vector.alloc(dim)
	diag.set_all(d)
	rhs = GSL::Vector.indgen(dim) + 1
	offdiag = GSL::Vector.alloc(dim-1)
	offdiag.set_all(od)

	x = GSL::Linalg::solve_symm_tridiag(diag, offdiag, rhs)
	p x 
	p actual
	for i in 0...dim do
		si = x[i]
		ai = actual[i]
    foo = check(si, ai, eps);
    if foo > 0
      printf("%3lu[%lu]: %22.18g   %22.18g\n", dim, i,s[i], actual[i]);
    end
    s += foo;
  end
  return s
end

def test_TDS_solve()
	s = 0
	
  actual =  GSL::Vector[0.0, 2.0]
  f = test_TDS_solve_dim(2, 1.0, 0.5, actual, 8.0 * GSL::DBL_EPSILON)
  GSL::test(f, "  solve_TDS dim=2 A");
  s += f;

  actual =  GSL::Vector[3.0/8.0, 15.0/8.0]
  f = test_TDS_solve_dim(2, 1.0, 1.0/3.0, actual, 8.0 * GSL::DBL_EPSILON);
  GSL::test(f, "  solve_TDS dim=2 B");
  s += f

  actual =  GSL::Vector[5.0/8.0, 9.0/8.0, 2.0, 15.0/8.0, 35.0/8.0]
  f = test_TDS_solve_dim(5, 1.0, 1.0/3.0, actual, 8.0 * GSL::DBL_EPSILON);
  GSL::test(f, "  solve_TDS dim=5");
end

def test_TDS_cyc_solve_one(dim, d, od, r, actual, eps)
	s = 0
	
	diag = d.duplicate()
	offdiag = od.duplicate()
	rhs = r.duplicate()

	x = GSL::Linalg::solve_symm_cyc_tridiag(diag, offdiag, rhs)
	p x
	p actual
	for i in 0...dim do
		si = x[i]
		ai = actual[i]
    foo = check(si, ai, eps);
    if foo > 0
      printf("%3lu[%lu]: %22.18g   %22.18g\n", dim, i, si, actual[i]);
    end
    s += foo;
  end
  return s;
end

def test_TDS_cyc_solve()
	s = 0
=begin
	dim = 1
	diag = GSL::Vector.alloc(1)
	diag[0] = 2
	offdiag = GSL::Vector.alloc(1)
	offdiag[0] = 3
	rhs = GSL::Vector.alloc(1)
	rhs[0] = 7
	actual = GSL::Vector.alloc(1)
	actual[0] = 3.5
    
  f = test_TDS_cyc_solve_one(dim, diag, offdiag, rhs, actual, 28.0 * GSL::DBL_EPSILON);
  GSL::test(f, "  solve_TDS_cyc dim=%lu A", dim);
  s += f;

	dim = 2
	diag = GSL::Vector[1, 2]
	offdiag = GSL::Vector[3, 4]
	rhs = GSL::Vector[7, -7]
	actual = GSL::Vector[-5, 4]

  f = test_TDS_cyc_solve_one(dim, diag, offdiag, rhs, actual, 28.0 * GSL::DBL_EPSILON);
  GSL::test(f, "  solve_TDS_cyc dim=%lu A", dim);
  s += f;
=end
	dim = 3
	diag = GSL::Vector[1, 1, 1]
	offdiag = GSL::Vector[3, 3, 3]
	rhs = GSL::Vector[7, -7, 7]
	actual = GSL::Vector[-2, 5, -2]

  f = test_TDS_cyc_solve_one(dim, diag, offdiag, rhs, actual, 28.0 * GSL::DBL_EPSILON);
  GSL::test(f, "  solve_TDS_cyc dim=#{dim} A");
  s += f;

	dim = 5
	diag = GSL::Vector[4, 2, 1, 2, 4]
	offdiag = GSL::Vector[1, 1, 1, 1, 1]
	rhs = GSL::Vector[30, -24, 3, 21, -30]
	actual = GSL::Vector[ 12, 3, -42, 42, -21 ]

  f = test_TDS_cyc_solve_one(dim, diag, offdiag, rhs, actual, 35.0 * GSL::DBL_EPSILON);
  GSL::test(f, "  solve_TDS_cyc dim=#{dim} B");
  s += f;
	return s
end

test_TDS_solve()
test_TDS_cyc_solve()

