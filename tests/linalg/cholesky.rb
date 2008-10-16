#!/usr/bin/env ruby
require("gsl")
require("../gsl_test2.rb")
include GSL::Test
include Linalg

m = GSL::Matrix.pascal(6)

p m
c_exp = GSL::Matrix[[1, 0, 0, 0, 0, 0],
               [1, 1, 0, 0, 0, 0],
               [1, 2, 1, 0, 0, 0],
               [1, 3, 3, 1, 0, 0],
               [1, 4, 6, 4, 1, 0],
               [1, 5, 10, 10, 5, 1]]

c = m.cholesky_decomp
a = c.lower
test2(a == c_exp, "#{m.class}#cholesky_decomp") 
test2((a*a.trans) == m, "#{m.class}#cholesky_decomp")

require("test/unit")
class CholeskyTest < Test::Unit::TestCase
	Data7 = GSL::Vector.alloc(66,0, 0,64, 126,63, 124,-62, 61,-61, 60,60, 0,-59,
                     0,-64, 65,0, 62,-124, -61,-122, -60,-60, 59,-59, -58,0,
                     126,-63, 62,124, 308,0, 180,-240, 59,-177, 174,58, -57,
                     -114,
                     124,62, -61,122, 180,240, 299,0, 174,-58, 57,171, 56,-112,
                     61,61, -60,60, 59,177, 174,58, 119,0, 0,112, 55,-55,
                     60,-60, 59,59, 174,-58, 57,-171, 0,-112, 116,0, -54,-54,
                     0,59, -58,0, -57,114, 56,112, 55,55, -54,54, 60,0).to_complex2.matrix_view(7, 7)
	Data7_sol = GSL::Vector.alloc(-0.524944196428570,0.209123883928571,
                         1.052873883928572,0.712444196428571,
                         0.117568824404762,0.443191964285714,
                         0.412862723214286,-0.356696428571429,
                         0.815931919642858,-0.265820312500000,
                         0.777929687500000,0.119484747023810,
                         1.058733258928571,-0.132087053571429).to_complex2
	RHS = GSL::Vector.alloc(1, 2, 3, 4, 5, 6, 7).to_complex                         
	def test_linalg_complex_cholesky_solve
		c = GSL::Linalg::Complex::Cholesky::decomp(Data7)
		x = GSL::Linalg::Complex::Cholesky::solve(c, RHS)
		assert_equal(x, Data7_sol)
	end
	def test_linalg_complex_cholesky_svx
		c = GSL::Linalg::Complex::Cholesky::decomp(Data7)
		x = RHS.clone
		GSL::Linalg::Complex::Cholesky::svx(c, x)
		assert_equal(x, Data7_sol)
	end	
end
