#!/usr/bin/env ruby
require("gsl")
require("../gsl_test2.rb")
require("./linalg.rb")
include GSL::Test

def test_HH_solve_dim(m, actual, eps)
  dim = m.size1
  s = 0
  _perm = GSL::Permutation.alloc(dim)
  x = GSL::Vector.indgen(dim) + 1
  hh = m.duplicate
  GSL::Linalg::HH.svx(hh, x)
  for i in 0...dim do
    foo = check(x[i],actual[i],eps)
    if foo > 0
      printf("%3lu[%lu]: %22.18g   %22.18g\n", dim, i, x[i], actual[i])
    end
    s += foo;
  end
  return s
end

def test_HH_solve()
	s = 0
  f = test_HH_solve_dim(Hilb2, Hilb2_solution, 8.0 * GSL::DBL_EPSILON);

  GSL::test(f, "  HH_solve Hilbert(2)");
  s += f;

  f = test_HH_solve_dim(Hilb3, Hilb3_solution, 128.0 * GSL::DBL_EPSILON);
  GSL::test(f, "  HH_solve Hilbert(3)");
  s += f;

  f = test_HH_solve_dim(Hilb4, Hilb4_solution, 2.0 * 1024.0 * GSL::DBL_EPSILON);
  GSL::test(f, "  HH_solve Hilbert(4)");
  s += f;

  f = test_HH_solve_dim(Hilb12, Hilb12_solution, 0.5);
  GSL::test(f, "  HH_solve Hilbert(12)");
  s += f;

  f = test_HH_solve_dim(Vander2, Vander2_solution, 8.0 * GSL::DBL_EPSILON);
  GSL::test(f, "  HH_solve Vander(2)");
  s += f;

  f = test_HH_solve_dim(Vander3, Vander3_solution, 64.0 * GSL::DBL_EPSILON);
  GSL::test(f, "  HH_solve Vander(3)");
  s += f;

  f = test_HH_solve_dim(Vander4, Vander4_solution, 1024.0 * GSL::DBL_EPSILON);
  GSL::test(f, "  HH_solve Vander(4)");
  s += f;

  f = test_HH_solve_dim(Vander12, Vander12_solution, 0.05);
  GSL::test(f, "  HH_solve Vander(12)");
  s += f;

  return s;
end

test_HH_solve()


