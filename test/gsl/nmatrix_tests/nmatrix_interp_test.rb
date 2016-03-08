require 'test_helper'

class NMatrixInterpTest < GSL::TestCase
  def test_interp_init_eval
    n = 10
    x = NMatrix.new([n], dtype: :float64)
    y = x.clone_structure

    0.upto(n-1) do |i|
      a = i.to_f
      x[i] = i + 0.5*Math::sin(a)
      y[i] = i + Math::cos(a*a)
    end

    interp = GSL::Interp.alloc("akima", n)
    interp.init(x, y)

    yi = []
    xi = x[0]
    r = []
    while xi < x[9]
      r << xi
      xi += 0.01
    end

    yi = interp.eval(x,y,NMatrix.new([r.size], r, dtype: :float64))

    assert_rel yi[1] , 1.0066, 0.001, 'yi[1]'
    assert_rel yi.to_a.last, 9.7618, 0.001, 'yi.last'
  end

  def test_spline_init_eval
    n = 10
    x = NMatrix.new([n], (1..10).to_a, dtype: :float64)
    y = x.dup

    spline = GSL::Spline.alloc(x, y)
    
    xi = NMatrix.new([n], (1..9).map { |a| a += 0.5 }, dtype: :float64)
    yi = spline.eval(xi)

    assert yi.class == NMatrix
    assert_rel yi[0], 1.5, 0.0001, 'yi[0]'
  end
end