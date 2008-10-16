#!/usr/bin/env ruby
require("gsl")
require("../gsl_test2.rb")
require("./linalg.rb")
include GSL::Test


def test_bidiag_decomp_dim(m, eps)
	s = 0
	mm = m.size1
	nn = m.size2

	a = m.duplicate
	b = GSL::Matrix.calloc(nn, nn)
	aa, tau1, tau2 = GSL::Linalg::Bidiag::decomp(a)
	u, v, d, sd = GSL::Linalg::Bidiag::unpack(aa, tau1, tau2)

	b.set_diagonal(d)
	for i in 0...(nn-1) do
		b[i][i+1] = sd[i]
	end

	a = u*b*v.trans  

	for i in 0...mm do
		for j in 0...nn do
			aij = a[i][j]
			mij = m[i][j]
      foo = check(aij, mij, eps);
      if foo > 0
        printf("(%3lu,%3lu)[%lu,%lu]: %22.18g   %22.18g\n", M, N, i,j, aij, mij);
      end
      s += foo;
    end
  end
  return s
end

def test_bidiag_decomp()
	s = 0
	
	m53 = create_general_matrix(5,3)
	m97 = create_general_matrix(9,7)
		
  f = test_bidiag_decomp_dim(m53, 2 * 64.0 * GSL::DBL_EPSILON);
  GSL::test(f, "  bidiag_decomp m(5,3)");
  s += f;

  f = test_bidiag_decomp_dim(m97, 2 * 64.0 * GSL::DBL_EPSILON);
  GSL::test(f, "  bidiag_decomp m(9,7)");
  s += f;

  f = test_bidiag_decomp_dim(Hilb2, 2 * 8.0 * GSL::DBL_EPSILON);
  GSL::test(f, "  bidiag_decomp hilbert(2)");
  s += f;

  f = test_bidiag_decomp_dim(Hilb3, 2 * 64.0 * GSL::DBL_EPSILON);
  GSL::test(f, "  bidiag_decomp hilbert(3)");
  s += f;

  f = test_bidiag_decomp_dim(Hilb4, 2 * 1024.0 * GSL::DBL_EPSILON);
  GSL::test(f, "  bidiag_decomp hilbert(4)");
  s += f;

  f = test_bidiag_decomp_dim(Hilb12, 2 * 1024.0 * GSL::DBL_EPSILON);
  GSL::test(f, "  bidiag_decomp hilbert(12)");
  s += f;

  return s;
end

test_bidiag_decomp()

