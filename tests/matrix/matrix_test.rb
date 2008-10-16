#!/usr/bin/enm ruby

require("gsl")
require("test/unit")

class MatrixTest < Test::Unit::TestCase

	def test_matrix_ispos_neg
		m = GSL::Matrix::Int.alloc([0, 1, 2, 3, 4, 5, 6, 7, 8, 9], 3, 3)
		assert_equal(m.ispos, 0)
		assert_equal(m.ispos?, false)		
		assert_equal(m.isneg, 0)
		assert_equal(m.isneg?, false)
		
		m += 1
		assert_equal(m.ispos, 1)
		assert_equal(m.ispos?, true)		
		assert_equal(m.isneg, 0)
		assert_equal(m.isneg?, false)		
		
		m -= 100
		assert_equal(m.ispos, 0)
		assert_equal(m.ispos?, false)		
		assert_equal(m.isneg, 1)
		assert_equal(m.isneg?, true)				
	end
		
	def test_matrix_isnonneg
		m = GSL::Matrix::Int.alloc([0, 1, 2, 3, 4, 5, 6, 7, 8, 9], 3, 3)
		assert_equal(m.isnonneg, 1)
		assert_equal(m.isnonneg?, true)		
		assert_equal(m.isneg, 0)
		assert_equal(m.isneg?, false)
		
		m -= 100
		assert_equal(m.isnonneg, 0)
		assert_equal(m.isnonneg?, false)		
		assert_equal(m.isneg, 1)
		assert_equal(m.isneg?, true)		
		
		m += 200
		assert_equal(m.isnonneg, 1)
		assert_equal(m.isnonneg?, true)		
		assert_equal(m.ispos, 1)
		assert_equal(m.ispos?, true)				
	end
end

