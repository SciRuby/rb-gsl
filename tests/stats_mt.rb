require 'minitest/autorun'

require 'gsl'

class StatsTest < MiniTest::Unit::TestCase

  def test_variance_with_fixed_mean
    v = GSL::Vector[1..8]
    assert_raises(ArgumentError, 'check for no args') do
      # This exposes a segfault(!)
      v.variance_with_fixed_mean
    end
  end

end

