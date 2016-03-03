require 'test_helper'

class BsplineTest < GSL::TestCase

  N = 100

  NMAX = 10
  BMAX = 10  # 100

  def _test_bspline(bw)
    ncoeffs, order = bw.ncoeffs, bw.order

    a = bw.breakpoint(0)
    b = bw.breakpoint(bw.nbreak - 1)

    N.times { |i|
      xi, sum = a + (b - a) * (i / (N - 1)), 0.0

      bb = bw.eval(xi)

      ncoeffs.times { |j|
        bj = bb[j]
        refute bj < 0 || bj > 1,
          "basis-spline coefficient #{j} is in range [0,1] for x=#{xi}"

        sum += bj
      }

      assert_rel sum, 1.0, order * GSL::DBL_EPSILON,
        "basis-spline coefficient #{order} is in range [0,1] for x=#{xi}"
    }
  end

  def test_bspline_knots_uniform
    1.upto(NMAX) { |order|
      2.upto(BMAX) { |breakpoints|
        bw = GSL::BSpline.alloc(order, breakpoints)
        bw.knots_uniform(-1.23 * order, 45.6 * order)
        _test_bspline(bw)
      }
    }
  end

  def test_bspline_knots
    1.upto(NMAX) do |order|
      2.upto(BMAX) do |breakpoints|
        a, b = -1.23 * order, 45.6 * order

        bw = GSL::BSpline.alloc(order, breakpoints)
        test_data = [GSL::Vector.alloc(breakpoints)]
        test_data << NMatrix.new([breakpoints], dtype: :float64) if ENV['NMATRIX']
        test_data.each do |k|
          breakpoints.times do |i|
            f = GSL.sqrt(i.to_f / (breakpoints - 1.0))
            k[i] = (1 - f) * a + f * b
          end

          bw.knots(k)
          _test_bspline(bw)
        end
      end
    end
  end
end
