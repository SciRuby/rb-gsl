require 'test_helper'

class InterpTest < GSL::TestCase

  def test_bsearch
    x_array = GSL::Vector.alloc(0.0, 1.0, 2.0, 3.0, 4.0)

    res = GSL::Interp.bsearch(x_array, 1.5, 0, 4)
    refute res != 1, 'simple bsearch'

    res = x_array.bsearch(4.0, 0, 4)
    refute res != 3, 'upper endpoint bsearch'

    res = GSL::Interp.bsearch(x_array, 0.0, 0, 4)
    refute res != 0, 'lower endpoint bsearch'

    res = GSL::Interp.bsearch(x_array, 2.0, 0, 4)
    refute res != 2, 'degenerate bsearch'

    res = GSL::Interp.bsearch(x_array, 10.0, 0, 4)
    refute res != 3, 'out of bounds bsearch +'

    res = GSL::Interp.bsearch(x_array, -10.0, 0, 4)
    refute res != 0, 'out of bounds bsearch -'
  end

end
