#!/usr/bin/env ruby

require("gsl")
require("test/unit")

class VectorComplexTest < Test::Unit::TestCase
	def test_vector_complex_get
		v = GSL::Vector::Complex.indgen(5)
		assert_equal(GSL::Vector::Complex[[3,0],[1,0],[2,0]], v.get([3, 1, 2]))
	end
	
	def test_vector_complex_addsub
		a = GSL::Vector::Complex[[-2,  5], [ 4, -1]]
		b = GSL::Vector::Complex[[10, 30], [20, 40]]
		c = GSL::Vector::Complex[[ 8, 35], [24, 39]]
		d = GSL::Vector::Complex[[12, 25], [16, 41]]
		assert_equal(c, a+b)
		assert_equal(d, b-a)		
	end
	
	def test_vector_complex_collect
		v = GSL::Vector::Complex.indgen(5)
		u = GSL::Vector::Complex[[0,0], [1,0], [4,0], [9,0], [16,0]]
		w = v.collect { |val| val*val }
		assert_equal(u, w)
	end

  def test_vector_complex_subvector
    v = GSL::Vector::Complex.indgen(12)

    # args = []
    vv = v.subvector
    assert_not_equal(v.object_id, vv.object_id)
    assert_equal(v.subvector, v)

    # args = [Fixnum]
    vv = v.subvector(3)
    assert_equal([0, 0, 1, 0, 2, 0], vv.to_a)
    assert_raise(ArgumentError) {v.subvector(-1)}

    # args = [Fixnum, Fixnum]
    vv = v.subvector(2, 3)
    assert_equal([2, 0, 3, 0, 4, 0], vv.to_a)

    vv = v.subvector(-4, 3)
    assert_equal([8, 0, 9, 0, 10, 0], vv.to_a)
    assert_raise(ArgumentError) {v.subvector(2, -1)}

    # args = [Fixnum, Fixnum, Fixnum]
    vv = v.subvector(1, 3, 4)
    assert_equal([1, 0, 4, 0, 7, 0, 10, 0], vv.to_a)

    # args = [Range]
    tests = {
    # ( range ) => [0, 0, 1, 0, 2, 0, 3, 0, 4, 0, 5, 0, 6, 0, 7, 0, 8, 0, 9, 0, 10, 0, 11, 0]
      ( 1..  3) => [      1, 0, 2, 0, 3, 0                                                  ],                                   
      ( 1... 3) => [      1, 0, 2, 0,                                                       ],
      ( 3..  1) => [      3, 0, 2, 0, 1, 0                                                  ],
      ( 3... 1) => [            3, 0, 2, 0                                                  ],
      (-7..  9) => [                              5, 0, 6, 0, 7, 0, 8, 0, 9, 0              ],
      (-7... 9) => [                              5, 0, 6, 0, 7, 0, 8, 0                    ],
      ( 4.. -3) => [                        4, 0, 5, 0, 6, 0, 7, 0, 8, 0, 9, 0              ],
      ( 4...-3) => [                        4, 0, 5, 0, 6, 0, 7, 0, 8, 0                    ],
      ( 2.. -2) => [            2, 0, 3, 0, 4, 0, 5, 0, 6, 0, 7, 0, 8, 0, 9, 0, 10, 0       ],
      ( 2...-2) => [            2, 0, 3, 0, 4, 0, 5, 0, 6, 0, 7, 0, 8, 0, 9, 0              ],
      (-2..  2) => [           10, 0, 9, 0, 8, 0, 7, 0, 6, 0, 5, 0, 4, 0, 3, 0,  2, 0       ],
      (-2... 2) => [           10, 0, 9, 0, 8, 0, 7, 0, 6, 0, 5, 0, 4, 0, 3, 0              ],
      (-3.. -1) => [                                                      9, 0, 10, 0, 11, 0],
      (-3...-1) => [                                                      9, 0, 10, 0       ],
      (-1.. -3) => [                                                     11, 0, 10, 0,  9, 0],
      (-1...-3) => [                                                     11, 0, 10, 0       ],
      # Add more test cases here...
    }
    tests.each do |r, x|
      assert_nothing_raised("subvector(#{r})") {v.subvector(r)}
      assert_equal(x, v.subvector(r).to_a, "subvector(#{r})")
    end

    # args = [Range, Fixnum]
    tests = {
    # [( range ), s] => [0, 0, 1, 0, 2, 0, 3, 0, 4, 0, 5, 0, 6, 0, 7, 0, 8, 0, 9, 0, 10, 0, 11, 0]
      [( 1..  6), 2] => [      1, 0,       3, 0,       5, 0                                      ],
      [( 1... 6), 2] => [      1, 0,       3, 0,       5, 0                                      ],
      [( 0..  6), 3] => [0, 0,             3, 0,             6, 0                                ],
      [( 0... 6), 3] => [0, 0,             3, 0                                                  ],
      # Add more test cases here...
    }
    tests.each do |(r,s), x|
      assert_nothing_raised("subvector(#{r},#{s})") {v.subvector(r)}
      assert_equal(x, v.subvector(r,s).to_a, "subvector(#{r},#{s})")
    end
  end
end
