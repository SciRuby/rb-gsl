require 'test_helper'

class Spline2dTest < GSL::TestCase
  def saddle(x,y)
    x*x - y*y
  end

  def setup
    x_samples =(-20..20).map(&:to_f).to_a
    y_samples = (-20..20).map(&:to_f).to_a
    z_samples = []
    x_samples.each do |x|
      y_samples.each do |y|
        z_samples << saddle(x,y)
      end
    end

    @x_array = GSL::Vector.alloc(x_samples)
    @y_array = GSL::Vector.alloc(y_samples)
    @z_array = GSL::Vector.alloc(z_samples)

    @i2d_bicubic = GSL::Spline2d.alloc(GSL::Interp2d::BICUBIC,
      @x_array, @y_array, @z_array)

    @i2d_bilinear = GSL::Spline2d.alloc(GSL::Interp2d::BILINEAR,
      @x_array, @y_array, @z_array)
  end

  def test_info
    str = <<-EOF
Class:      GSL::Spline2d
SuperClass: GSL::Object
Type:       bilinear
xmin:       -20.000000
xmax:       20.000000
ymin:       -20.000000
ymax:       20.000000
xSize:       41
ySize:       41
EOF
    assert_equal str, @i2d_bilinear.info
  end

  def test_use_case_saddle_interpolation
    tolerance = 0.05 # 5% inaccuracy is tolerated in tests below
    interpolator_type = GSL::Interp2d::BICUBIC

    i2d = GSL::Spline2d.alloc(interpolator_type, @x_array, @y_array, @z_array)

    # confirm that the fit passes very close to the sampled data points
    @x_array.each do |x|
      @y_array.each do |y|
        expected_z = saddle(x,y)
        z = i2d.eval(x, y)
        error = (z - expected_z).abs
        max_error = (expected_z.abs)*tolerance
        if max_error == 0
          max_error = tolerance
        end
        refute error > max_error, "error @ sample #{x},#{y} z: #{z} expected_z: #{expected_z}"
      end
    end

    interstitial_x_values = @x_array.to_a.first(@x_array.size-1).map {|v| v+ 0.5}
    interstitial_y_values = @y_array.to_a.first(@y_array.size-1).map {|v| v+ 0.3}

    # confirm that interstitial values are interpolated accurately 
    interstitial_x_values.each do |x|
      interstitial_y_values.each do |y|
        expected_z = saddle(x,y)
        z = i2d.eval(x,y)
        error = (z - expected_z).abs
        max_error = (expected_z.abs)*tolerance
        if max_error == 0
          max_error = tolerance
        end
        refute error > max_error, "error @ interstitial #{x},#{y}"
      end
    end
  end
end if GSL::GSL_VERSION >= '2.0'