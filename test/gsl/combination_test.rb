require 'test_helper'

class CombinationTest < GSL::TestCase

  def setup
    @c63 = GSL::Matrix.alloc(
      [0, 1, 2], [0, 1, 3], [0, 1, 4], [0, 1, 5],
      [0, 2, 3], [0, 2, 4], [0, 2, 5],
      [0, 3, 4], [0, 3, 5],
      [0, 4, 5],
      [1, 2, 3], [1, 2, 4], [1, 2, 5],
      [1, 3, 4], [1, 3, 5],
      [1, 4, 5],
      [2, 3, 4], [2, 3, 5], [2, 4, 5],
      [3, 4, 5]
    )
  end

  def test_6_3
    c = GSL::Combination.alloc(6, 3)
    c.init_first

    status, i = false, 0

    begin
      if i >= 20
        status = true
        break
      end

      3.times { |j| status |= c.data[j] != @c63[i, j] }
      refute c.valid?, 'GSL::Combination#valid(%u)' % i

      i += 1
    end while c.next == GSL::SUCCESS

    refute status, 'GSL::Combination#next, 6 choose 3 combination, 20 steps'

    3.times { c.next }
    3.times { |j| status |= c.data[j] != @c63[19, j] }

    refute status, 'GSL::Combination#prev on the first combination'
    refute c.valid?, 'GSL::Combination#valid on the last combination'

    d = GSL::Combination.alloc(6, 3)
    GSL::Combination.memcpy(d, c)

    status = false
    3.times { |j| status |= d.data[j] != c.data[j] }
    refute status, 'GSL::Combination.memcpy, 6 choose 3 combination'

    c.init_last
    i = 20

    begin
      if i == 0
        status = true
        break
      end

      i -= 1

      3.times { |j| status |= c.data[j] != @c63[i, j] }
      refute c.valid?, 'GSL::Combination#valid(%u)' % i
    end while c.prev == GSL::SUCCESS

    refute status, 'GSL::Combination#prev, 6 choose 3 combination, 20 steps'

    3.times { c.prev }
    3.times { |j| status |= c.data[j] != @c63[0, j] }

    refute status, 'GSL::Combination#prev on the first combination'
    refute c.valid?, 'GSL::Combination#valid on the first combination'

    d = GSL::Combination.alloc(6, 3)
    GSL::Combination.memcpy(d, c)

    status = false
    3.times { |j| status |= d.data[j] != c.data[j] }
    refute status, 'GSL::Combination.memcpy, 6 choose 3 combination'
  end

  def test_7_0
    c, desc = GSL::Combination.calloc(7, 0), 'GSL::Combination 7 choose 0'
    2.times { assert c.next == GSL::FAILURE, desc }
    2.times { assert c.prev == GSL::FAILURE, desc }
  end

  def test_7_7
    c, desc = GSL::Combination.calloc(7, 7), 'GSL::Combination 7 choose 7'

    3.times {
      7.times { |j| assert c.get(j) == j, desc }
      assert c.next == GSL::FAILURE, desc
    }

    7.times { |j| assert c.get(j) == j, desc }
  end

end
