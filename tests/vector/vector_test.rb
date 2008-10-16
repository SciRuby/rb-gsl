#!/usr/bin/env ruby

require("gsl")
require("test/unit")

class VectorTest < Test::Unit::TestCase
	def test_vector_get
		v = GSL::Vector::Int.indgen(5)
		assert_equal(GSL::Vector::Int[3, 1, 2], v.get([3, 1, 2]))
	end
	
	def test_vector_addsub
		a = GSL::Vector::Int[2, 5, 4]
		b = GSL::Vector::Int[10, 30, 20]
		c = GSL::Vector::Int[12, 35, 24]
		d = GSL::Vector::Int[8, 25, 16]
		assert_equal(c, a+b)
		assert_equal(d, b-a)		
	end
	
	def test_vector_collect
		v = GSL::Vector::Int.indgen(5)
		u = GSL::Vector::Int[0, 1, 4, 9, 16]
		w = v.collect { |val| val*val }
		assert_equal(u, w)
	end

	def test_vector_ispos_neg
		v = GSL::Vector::Int.indgen(5)
		assert_equal(v.ispos, 0)
		assert_equal(v.ispos?, false)		
		assert_equal(v.isneg, 0)
		assert_equal(v.isneg?, false)
		
		v += 1
		assert_equal(v.ispos, 1)
		assert_equal(v.ispos?, true)		
		assert_equal(v.isneg, 0)
		assert_equal(v.isneg?, false)		
		
		v -= 100
		assert_equal(v.ispos, 0)
		assert_equal(v.ispos?, false)		
		assert_equal(v.isneg, 1)
		assert_equal(v.isneg?, true)				
	end
		
	def test_vector_isnonneg
		v = GSL::Vector::Int.indgen(5)
		assert_equal(v.isnonneg, 1)
		assert_equal(v.isnonneg?, true)		
		assert_equal(v.isneg, 0)
		assert_equal(v.isneg?, false)
		
		v -= 100
		assert_equal(v.isnonneg, 0)
		assert_equal(v.isnonneg?, false)		
		assert_equal(v.isneg, 1)
		assert_equal(v.isneg?, true)		
		
		v += 200
		assert_equal(v.isnonneg, 1)
		assert_equal(v.isnonneg?, true)		
		assert_equal(v.ispos, 1)
		assert_equal(v.ispos?, true)				
	end
end

