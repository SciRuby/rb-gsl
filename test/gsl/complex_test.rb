require 'test_helper'

class ComplexTest < GSL::TestCase

  def test_complex
    10.times { |i|
      r = (i - 5.0) * 0.3
      t = 2.0 * GSL::M_PI * i / 5.0

      z = GSL::Complex.polar(r, t)

      assert_rel z.real, r * Math.cos(t), 10 * GSL::DBL_EPSILON,
        'gsl_complex_polar real part at (r=%g,t=%g)' % [r, t]

      assert_rel z.imag, r * Math.sin(t), 10 * GSL::DBL_EPSILON,
        'gsl_complex_polar imag part at (r=%g,t=%g)' % [r, t]
    }
  end
  
  # Test if it is possible to create a GSL::Complex from ::Complex
  def test_rb_complex_creation
    rb_comp = Complex(rand, rand)
    
    z = GSL::Complex.alloc(rb_comp)
    
    assert_rel z.real, rb_comp.real, GSL::DBL_EPSILON,
      "gsl_complex real part.  Re(#{rb_comp}) = #{z.real}"
    assert_rel z.imag, rb_comp.imag, GSL::DBL_EPSILON,
      "gsl_complex imag part.  Im(#{rb_comp}) = #{z.imag}"
  end
end
