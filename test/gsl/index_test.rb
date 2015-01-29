require 'test_helper'

class IndexTest < GSL::TestCase

  # helper(s)

  def _create_index array
    i = GSL::Index.alloc(array.size)
    array.each_with_index { |e,idx| i[idx] = e }
    i
  end

  # tests

  def test_get_int
    i = GSL::Index.alloc(5)
    assert_equal 2, i.get(2)
    assert_equal 3, i.get(-2)
  end

  def test_get_array
    i = GSL::Index.alloc(5)
    assert_equal _create_index([2, 3]), i.get([2, 3])
    assert_equal _create_index([4, 2]), i.get([-1, 2])
    assert_equal _create_index([4, 3, 1]), i.get([4, -2, 1])
  end

  def test_get_range
    i = GSL::Index.alloc(5)
    assert_equal _create_index([2, 3]), i.get(2..3)
    assert_equal _create_index([1, 2, 3]), i.get(-4...-1) # note the exclusive range operator!
  end

  def test_get_failure
    i = GSL::Index.alloc(5)
    assert_nothing_raised('get(4)') { i.get(4) }
    assert_raises(RangeError) { i.get(5) }
    assert_raises(RangeError) { i.get(1_000_000) }
    assert_raises(ArgumentError) { i.get(10**100) }

    assert_nothing_raised('get(-5)') { i.get(-5) }
    assert_raises(RangeError) { i.get(-6) }

    assert_nothing_raised('get([0, 4, -1, -5])') { i.get([0, 4, -1, -5]) }
    assert_raises(RangeError) { i.get([5]) }
    assert_raises(RangeError) { i.get([-6]) }
    assert_raises(RangeError) { i.get([-6, 0, 5]) }

    assert_nothing_raised('get(0..4)') { i.get(0..4) }
    assert_nothing_raised('get(-5..-1)') { i.get(-5..-1) }
    assert_raises(RangeError) { i.get(0..5) }
    assert_raises(RangeError) { i.get(-6..-1) }
    #assert_raises(RangeError) { i.get(-5..0) } # ???
    #assert_raises(RangeError) { i.get(3..-5) } # ??? (3..-5).to_a == nil anyway
  end

end
