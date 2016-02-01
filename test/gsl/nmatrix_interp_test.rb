require 'test_helper'

class NMatrixEigenTest < GSL::TestCase
  def test_interp_init_eval
    n = 10
    x = NMatrix.new([n], dtype: :float64)
    y = x.clone_structure

    0.upto(n-1) do |i|
      a = i.to_f
      x[i] = i + 0.5*Math::sin(a)
      y[i] = i + Math::cos(a*a)
    end

    interp = Interp.alloc("akima", n)
    interp.init(x, y)

    xi = []
    yi = []

    xi << x[0]
    i = 0
    while xi < x[9]
      yi << interp.eval(x, y, xi)
      xi[i] = xi[i + 1] + 0.5
      i += 1
    end

    assert_rel yi[1] , 1.2685, 0.001
    assert_rel yi[19], 9.2816, 0.001
  end

  def test_spline_init_eval
    n = 10
    x = NMatrix.new([n], (1..10).to_a, dtype: :float64)
    y = x.dup

    spline = GSL::Spline.alloc(x, y)
    
    xi = NMatrix.new([n], (1..9).map { |a| a += 0.5 }, dtype: :float64)
    yi = spline.eval(xi)

    assert yi.class, NMatrix
    assert_rel yi[0], 1.5, 0.0001
  end
end