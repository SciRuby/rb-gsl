require_relative '../../test_helper'

class NMatrixChebTest < GSL::TestCase
  def test_eval
    f = GSL::Function.alloc do |x|
      if x < 0.5
        0.25
      else
        0.75
      end
    end

    n = 1000
    order = 40
    cs = GSL::Cheb.alloc(order)
    x = NMatrix.new([n], GSL::Vector.linspace(0, 1, n).to_a, dtype: :float64)

    ff = f.eval(x)

    assert ff.class == NMatrix

    cs.init(f, 0, 1)
    r10 = cs.eval_n(10, x)
    r40 = cs.eval(x)

    assert r10.class == NMatrix
    assert_rel r10.to_a.last, 0.758879, 0.001, 'test r10.last'
    assert_rel r10[5]  , 0.247816, 0.001, 'test r10[5]'

    assert r40.class == NMatrix
    assert_rel r40[5]   , 0.255682, 0.001, 'test r40[5]'
    assert_rel r40.first, 0.25633 , 0.001, 'test r40.first'
  end
end