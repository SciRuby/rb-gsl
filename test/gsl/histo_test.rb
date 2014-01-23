require 'test_helper'

class HistoTest < GSL::TestCase

  def test_histo
    h = GSL::Histogram.alloc(10, [0, 10])

    assert h
    assert h.get_range(2)
    assert h.range
    assert h.bin
  end

end
