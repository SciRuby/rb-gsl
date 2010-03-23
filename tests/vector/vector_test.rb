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

  def test_vector_subvector
    v = GSL::Vector::Int.indgen(12)

    # args = []
    vv = v.subvector
    assert_not_equal(v.object_id, vv.object_id)
    assert_equal(v.subvector, v)

    # args = [Fixnum]
    vv = v.subvector(3)
    assert_equal([0, 1, 2], vv.to_a)
#    assert_raise(ArgumentError) {v.subvector(-1)}

    # args = [Fixnum, Fixnum]
    vv = v.subvector(2, 3)
    assert_equal([2, 3, 4], vv.to_a)

    vv = v.subvector(-4, 3)
    assert_equal([8, 9, 10], vv.to_a)
#    assert_raise(ArgumentError) {v.subvector(2, -1)}

    # args = [Fixnum, Fixnum, Fixnum]
    vv = v.subvector(1, 3, 4)
    assert_equal([1, 4, 7, 10], vv.to_a)

    # args = [Range]
    tests = {
    # ( range ) => [0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11]
      ( 1..  3) => [   1, 2, 3                          ],                                   
      ( 1... 3) => [   1, 2                             ],
      ( 3..  1) => [   3, 2, 1                          ],
      ( 3... 1) => [      3, 2                          ],
      (-7..  9) => [               5, 6, 7, 8, 9        ],
      (-7... 9) => [               5, 6, 7, 8           ],
      ( 4.. -3) => [            4, 5, 6, 7, 8, 9        ],
      ( 4...-3) => [            4, 5, 6, 7, 8           ],
      ( 2.. -2) => [      2, 3, 4, 5, 6, 7, 8, 9, 10    ],
      ( 2...-2) => [      2, 3, 4, 5, 6, 7, 8, 9        ],
      (-2..  2) => [     10, 9, 8, 7, 6, 5, 4, 3,  2    ],
      (-2... 2) => [     10, 9, 8, 7, 6, 5, 4, 3        ],
      (-3.. -1) => [                           9, 10, 11],
      (-3...-1) => [                           9, 10    ],
      (-1.. -3) => [                          11, 10,  9],
      (-1...-3) => [                          11, 10    ],
      # Add more test cases here...
    }
    tests.each do |r, x|
#      assert_nothing_raised("subvector(#{r})") {v.subvector(r)}
#      assert_equal(x, v.subvector(r).to_a, "subvector(#{r})")
    end

    # args = [Range, Fixnum]
    tests = {
    # [( range ), s] => [0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11]
      [( 1..  6), 2] => [   1,    3,    5                    ],
      [( 1... 6), 2] => [   1,    3,    5                    ],
      [( 0..  6), 3] => [0,       3,      6                  ],
      [( 0... 6), 3] => [0,    3                             ],
      # Add more test cases here...
    }
    tests.each do |(r,s), x|
#      assert_nothing_raised("subvector(#{r},#{s})") {v.subvector(r)}
#      assert_equal(x, v.subvector(r,s).to_a, "subvector(#{r},#{s})")
    end
  end
end
