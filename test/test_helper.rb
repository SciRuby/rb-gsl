require 'test/unit'
begin
  require 'narray'
rescue LoadError => e
end

begin
  require 'nmatrix/nmatrix'
rescue LoadError => e
end
require 'gsl'

GSL::IEEE.env_setup
GSL::Rng.env_setup

class GSL::TestCase < Test::Unit::TestCase

  def assert_factor(result, expected, factor, desc)
    refute result == expected ? false : (
      expected.zero? ? result != expected :
      (u = result / expected; u > factor || u < 1.0 / factor)),
      '%s (%.18g observed vs %.18g expected)' % [desc, result, expected]
  end

  def assert_rel(result, expected, relerr, desc)
    refute((GSL.isnan?(result) || GSL.isnan?(expected)) ?
       GSL.isnan?(result) != GSL.isnan?(expected) :
      (GSL.isinf?(result) || GSL.isinf?(expected)) ?
       GSL.isinf?(result) != GSL.isinf?(expected) :
      expected.zero? ? result.abs > relerr :
      (result - expected).abs / expected.abs > relerr,
      '%s (%.18g observed vs %.18g expected)' % [desc, result, expected])
  end

  def assert_abs(result, expected, abserr, desc)
    refute((GSL.isnan?(result) || GSL.isnan?(expected)) ?
       GSL.isnan?(result) != GSL.isnan?(expected) :
      (GSL.isinf?(result) || GSL.isinf?(expected)) ?
       GSL.isinf?(result) != GSL.isinf?(expected) :
      (result - expected).abs > abserr,
      '%s (%.18g observed vs %.18g expected)' % [desc, result, expected])
  end

  def assert_int(result, expected, desc)
    assert(result == expected, '%s (%d observed vs %d expected)' % [desc, result, expected])
  end

  def assert_tol(a, b, msg)
    assert((a - b).abs < (self.class::EPSREL * GSL.MIN(a.abs, b.abs) + self.class::EPSABS), msg)
  end

  # Assert each element in an enumerable with absolute error.
  def assert_enum_abs(result, expected, abserr, desc)
    assert result.size == expected.size, 'size mismatch.'

    r_enum = result.each
    e_enum = expected.each

    while true
      begin
        res_value = r_enum.next
        exp_value = e_enum.next
        assert_abs res_value, exp_value, abserr, "value mismatch. #{res_value} should be #{exp_value}"
      rescue StopIteration
        break
      end
    end
  end

end
