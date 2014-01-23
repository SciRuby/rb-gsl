require 'test_helper'

class VectorTest < GSL::TestCase

  def test_get
    v = GSL::Vector::Int.indgen(5)
    assert_equal GSL::Vector::Int[3, 1, 2], v.get([3, 1, 2])
  end

  def test_addsub
    a = GSL::Vector::Int[2, 5, 4]
    b = GSL::Vector::Int[10, 30, 20]
    c = GSL::Vector::Int[12, 35, 24]
    d = GSL::Vector::Int[8, 25, 16]

    assert_equal c, a + b
    assert_equal d, b - a
  end

  def test_collect
    v = GSL::Vector::Int.indgen(5)
    u = GSL::Vector::Int[0, 1, 4, 9, 16]

    assert_equal u, v.collect { |val| val * val }
  end

  def test_ispos_neg
    v = GSL::Vector::Int.indgen(5)
    assert_equal 0,     v.ispos
    assert_equal false, v.ispos?
    assert_equal 0,     v.isneg
    assert_equal false, v.isneg?

    v += 1
    assert_equal 1,     v.ispos
    assert_equal true,  v.ispos?
    assert_equal 0,     v.isneg
    assert_equal false, v.isneg?

    v -= 100
    assert_equal 0,     v.ispos
    assert_equal false, v.ispos?
    assert_equal 1,     v.isneg
    assert_equal true,  v.isneg?
  end

  def test_isnonneg
    v = GSL::Vector::Int.indgen(5)
    assert_equal 1,     v.isnonneg
    assert_equal true,  v.isnonneg?
    assert_equal 0,     v.isneg
    assert_equal false, v.isneg?

    v -= 100
    assert_equal 0,     v.isnonneg
    assert_equal false, v.isnonneg?
    assert_equal 1,     v.isneg
    assert_equal true,  v.isneg?

    v += 200
    assert_equal 1,     v.isnonneg
    assert_equal true,  v.isnonneg?
    assert_equal 1,     v.ispos
    assert_equal true,  v.ispos?
  end

  def test_subvector
    v = GSL::Vector::Int.indgen(12)

    vv = v.subvector
    assert_not_equal v.object_id, vv.object_id
    assert_equal     v.subvector, v

    vv = v.subvector(3)
    assert_equal [0, 1, 2], vv.to_a
    assert_nothing_raised('subvector(-1)') { v.subvector(-1) }

    vv = v.subvector(-1)
    assert_equal [11], vv.to_a

    vv = v.subvector(-2)
    assert_equal [10, 11], vv.to_a
    assert_raises(RangeError) { v.subvector(-13) }

    vv = v.subvector(2, 3)
    assert_equal [2, 3, 4], vv.to_a

    vv = v.subvector(-4, 3)
    assert_equal [8, 9, 10], vv.to_a
    assert_nothing_raised('subvector(-4, -3)') { v.subvector(-4, -3) }

    vv = v.subvector(-4, -3)
    assert_equal [8, 7, 6], vv.to_a
    assert_raises(GSL::ERROR::EINVAL) { v.subvector(-11, -3) }

    vv = v.subvector(1, 3, 4)
    assert_equal [1, 4, 7, 10], vv.to_a

    {
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
      (-1...-3) => [                          11, 10    ]
    }.each { |r, x|
      assert_nothing_raised("subvector(#{r})") { v.subvector(r) }
      assert_equal x, v.subvector(r).to_a, "subvector(#{r})"
    }

    {
    # [( range ), s] => [0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11]
      [( 1..  6), 2] => [   1,    3,    5                    ],
      [( 1... 6), 2] => [   1,    3,    5                    ],
      [( 0..  6), 3] => [0,       3,      6                  ],
      [( 0... 6), 3] => [0,    3                             ]
    }.each { |(r, s), x|
      assert_nothing_raised("subvector(#{r},#{s})") { v.subvector(r) }
      assert_equal x, v.subvector(r,s).to_a, "subvector(#{r},#{s})"
    }
  end

  def test_complex_get
    v = GSL::Vector::Complex.indgen(5)
    assert_equal GSL::Vector::Complex[[3, 0], [1, 0], [2, 0]], v.get([3, 1, 2])
  end

  def test_complex_addsub
    a = GSL::Vector::Complex[[-2,  5], [ 4, -1]]
    b = GSL::Vector::Complex[[10, 30], [20, 40]]
    c = GSL::Vector::Complex[[ 8, 35], [24, 39]]
    d = GSL::Vector::Complex[[12, 25], [16, 41]]

    assert_equal c, a + b
    assert_equal d, b - a
  end

  def test_complex_collect
    v = GSL::Vector::Complex.indgen(5)
    u = GSL::Vector::Complex[[0, 0], [1, 0], [4, 0], [9, 0], [16, 0]]

    assert_equal u, v.collect { |val| val * val }
  end

  def test_complex_subvector
    v = GSL::Vector::Complex.indgen(12)

    vv = v.subvector
    assert_not_equal v.object_id, vv.object_id
    assert_equal     v.subvector, v

    vv = v.subvector(3)
    assert_equal [0, 0, 1, 0, 2, 0], vv.to_a
    assert_nothing_raised('subvector(-1)') { v.subvector(-1) }

    vv = v.subvector(-1)
    assert_equal [11, 0], vv.to_a

    vv = v.subvector(-2)
    assert_equal [10, 0, 11, 0], vv.to_a
    assert_raises(RangeError) { v.subvector(-13) }

    vv = v.subvector(2, 3)
    assert_equal [2, 0, 3, 0, 4, 0], vv.to_a

    vv = v.subvector(-4, 3)
    assert_equal [8, 0, 9, 0, 10, 0], vv.to_a
    assert_nothing_raised('subvector(-4, -3)') { v.subvector(-4, -3) }

    vv = v.subvector(-4, -3)
    assert_equal [8, 0, 7, 0, 6, 0], vv.to_a
    assert_raises(GSL::ERROR::EINVAL) { v.subvector(-11, -3) }

    vv = v.subvector(1, 3, 4)
    assert_equal [1, 0, 4, 0, 7, 0, 10, 0], vv.to_a

    {
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
      (-1...-3) => [                                                     11, 0, 10, 0       ]
    }.each { |r, x|
      assert_nothing_raised("subvector(#{r})") { v.subvector(r) }
      assert_equal x, v.subvector(r).to_a, "subvector(#{r})"
    }

    {
    # [( range ), s] => [0, 0, 1, 0, 2, 0, 3, 0, 4, 0, 5, 0, 6, 0, 7, 0, 8, 0, 9, 0, 10, 0, 11, 0]
      [( 1..  6), 2] => [      1, 0,       3, 0,       5, 0                                      ],
      [( 1... 6), 2] => [      1, 0,       3, 0,       5, 0                                      ],
      [( 0..  6), 3] => [0, 0,             3, 0,             6, 0                                ],
      [( 0... 6), 3] => [0, 0,             3, 0                                                  ]
    }.each { |(r, s), x|
      assert_nothing_raised("subvector(#{r},#{s})") { v.subvector(r) }
      assert_equal x, v.subvector(r,s).to_a, "subvector(#{r},#{s})"
    }
  end

end
