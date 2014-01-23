require 'test_helper'

class GSLTest < GSL::TestCase

  def test_expm1
    y = GSL.expm1(0.0)
    y_expected = 0.0
    assert_rel y, y_expected, 1e-15, 'GSL.expm1(0.0)'

    y = GSL.expm1(1e-10)
    y_expected = 1.000000000050000000002e-10
    assert_rel y, y_expected, 1e-15, 'GSL.expm1(1e-10)'

    y = GSL.expm1(-1e-10)
    y_expected = -9.999999999500000000017e-11
    assert_rel y, y_expected, 1e-15, 'GSL.expm1(-1e-10)'

    y = GSL.expm1(0.1)
    y_expected = 0.1051709180756476248117078264902
    assert_rel y, y_expected, 1e-15, 'GSL.expm1(0.1)'

    y = GSL.expm1(-0.1)
    y_expected = -0.09516258196404042683575094055356
    assert_rel y, y_expected, 1e-15, 'GSL.expm1(-0.1)'

    y = GSL.expm1(10.0)
    y_expected = 22025.465794806716516957900645284
    assert_rel y, y_expected, 1e-15, 'GSL.expm1(10.0)'

    y = GSL.expm1(-10.0)
    y_expected = -0.99995460007023751514846440848444
    assert_rel y, y_expected, 1e-15, 'GSL.expm1(-10.0)'
  end

  def test_log1p
    y = GSL.log1p(0.0)
    y_expected = 0.0
    assert_rel y, y_expected, 1e-15, 'GSL.log1p(0.0)'

    y = GSL.log1p(1e-10)
    y_expected = 9.9999999995000000000333333333308e-11
    assert_rel y, y_expected, 1e-15, 'GSL.log1p(1e-10)'

    y = GSL.log1p(0.1)
    y_expected = 0.095310179804324860043952123280765
    assert_rel y, y_expected, 1e-15, 'GSL.log1p(0.1)'

    y = GSL.log1p(10.0)
    y_expected = 2.3978952727983705440619435779651
    assert_rel y, y_expected, 1e-15, 'GSL.log1p(10.0)'
  end

  def test_hypot
    y = GSL.hypot(0.0, 0.0)
    y_expected = 0.0
    assert_rel y, y_expected, 1e-15, 'GSL.hypot(0.0, 0.0)'

    y = GSL.hypot(1e-10, 1e-10)
    y_expected = 1.414213562373095048801688e-10
    assert_rel y, y_expected, 1e-15, 'GSL.hypot(1e-10, 1e-10)'

    y = GSL.hypot(1e-38, 1e-38)
    y_expected = 1.414213562373095048801688e-38
    assert_rel y, y_expected, 1e-15, 'GSL.hypot(1e-38, 1e-38)'

    y = GSL.hypot(1e-10, -1.0)
    y_expected = 1.000000000000000000005
    assert_rel y, y_expected, 1e-15, 'GSL.hypot(1e-10, -1)'

    y = GSL.hypot(-1.0, 1e-10)
    y_expected = 1.000000000000000000005
    assert_rel y, y_expected, 1e-15, 'GSL.hypot(-1, 1e-10)'

    #y = GSL.hypot(1e307, 1e301)
    #y_expected = 1.000000000000499999999999e307
    #assert_rel y, y_expected, 1e-15, 'GSL.hypot(1e307, 1e301)'

    #y = GSL.hypot(1e301, 1e307)
    #y_expected = 1.000000000000499999999999e307
    #assert_rel y, y_expected, 1e-15, 'GSL.hypot(1e301, 1e307)'

    #y = GSL.hypot(1e307, 1e307)
    #y_expected = 1.414213562373095048801688e307
    #assert_rel y, y_expected, 1e-15, 'GSL.hypot(1e307, 1e307)'
  end

  def test_acosh
    y = GSL.acosh(1.0)
    y_expected = 0.0
    assert_rel y, y_expected, 1e-15, 'GSL.acosh(1.0)'

    y = GSL.acosh(1.1)
    y_expected = 4.435682543851151891329110663525e-1
    assert_rel y, y_expected, 1e-15, 'GSL.acosh(1.1)'

    y = GSL.acosh(10.0)
    y_expected = 2.9932228461263808979126677137742e0
    assert_rel y, y_expected, 1e-15, 'GSL.acosh(10.0)'

    y = GSL.acosh(1e10)
    y_expected = 2.3718998110500402149594646668302e1
    assert_rel y, y_expected, 1e-15, 'GSL.acosh(1e10)'
  end

  def test_asinh
    y = GSL.asinh(0.0)
    y_expected = 0.0
    assert_rel y, y_expected, 1e-15, 'GSL.asinh(0.0)'

    y = GSL.asinh(1e-10)
    y_expected = 9.9999999999999999999833333333346e-11
    assert_rel y, y_expected, 1e-15, 'GSL.asinh(1e-10)'

    y = GSL.asinh(-1e-10)
    y_expected = -9.9999999999999999999833333333346e-11
    assert_rel y, y_expected, 1e-15, 'GSL.asinh(1e-10)'

    y = GSL.asinh(0.1)
    y_expected = 9.983407889920756332730312470477e-2
    assert_rel y, y_expected, 1e-15, 'GSL.asinh(0.1)'

    y = GSL.asinh(-0.1)
    y_expected = -9.983407889920756332730312470477e-2
    assert_rel y, y_expected, 1e-15, 'GSL.asinh(-0.1)'

    y = GSL.asinh(1.0)
    y_expected = 8.8137358701954302523260932497979e-1
    assert_rel y, y_expected, 1e-15, 'GSL.asinh(1.0)'

    y = GSL.asinh(-1.0)
    y_expected = -8.8137358701954302523260932497979e-1
    assert_rel y, y_expected, 1e-15, 'GSL.asinh(-1.0)'

    y = GSL.asinh(10.0)
    y_expected = 2.9982229502979697388465955375965e0
    assert_rel y, y_expected, 1e-15, 'GSL.asinh(10)'

    y = GSL.asinh(-10.0)
    y_expected = -2.9982229502979697388465955375965e0
    assert_rel y, y_expected, 1e-15, 'GSL.asinh(-10)'

    y = GSL.asinh(1e10)
    y_expected = 2.3718998110500402149599646668302e1
    assert_rel y, y_expected, 1e-15, 'GSL.asinh(1e10)'

    y = GSL.asinh(-1e10)
    y_expected = -2.3718998110500402149599646668302e1
    assert_rel y, y_expected, 1e-15, 'GSL.asinh(-1e10)'
  end

  def test_atanh
    y = GSL.atanh(0.0)
    y_expected = 0.0
    assert_rel y, y_expected, 1e-15, 'GSL.atanh(0.0)'

    y = GSL.atanh(1e-20)
    y_expected = 1e-20
    assert_rel y, y_expected, 1e-15, 'GSL.atanh(1e-20)'

    y = GSL.atanh(-1e-20)
    y_expected = -1e-20
    assert_rel y, y_expected, 1e-15, 'GSL.atanh(-1e-20)'

    y = GSL.atanh(0.1)
    y_expected = 1.0033534773107558063572655206004e-1
    assert_rel y, y_expected, 1e-15, 'GSL.atanh(0.1)'

    y = GSL.atanh(-0.1)
    y_expected = -1.0033534773107558063572655206004e-1
    assert_rel y, y_expected, 1e-15, 'GSL.atanh(-0.1)'

    y = GSL.atanh(0.9)
    y_expected = 1.4722194895832202300045137159439e0
    assert_rel y, y_expected, 1e-15, 'GSL.atanh(0.9)'

    y = GSL.atanh(-0.9)
    y_expected = -1.4722194895832202300045137159439e0
    assert_rel y, y_expected, 1e-15, 'GSL.atanh(0.9)'
  end

  def test_pow_int
    y = GSL.pow_2(-3.14)
    y_expected = GSL.pow(-3.14, 2.0)
    assert_rel y, y_expected, 1e-15, 'GSL.pow_2(-3.14)'

    y = GSL.pow_3(-3.14)
    y_expected = GSL.pow(-3.14, 3.0)
    assert_rel y, y_expected, 1e-15, 'GSL.pow_3(-3.14)'

    y = GSL.pow_4(-3.14)
    y_expected = GSL.pow(-3.14, 4.0)
    assert_rel y, y_expected, 1e-15, 'GSL.pow_4(-3.14)'

    y = GSL.pow_5(-3.14)
    y_expected = GSL.pow(-3.14, 5.0)
    assert_rel y, y_expected, 1e-15, 'GSL.pow_5(-3.14)'

    y = GSL.pow_6(-3.14)
    y_expected = GSL.pow(-3.14, 6.0)
    assert_rel y, y_expected, 1e-15, 'GSL.pow_6(-3.14)'

    y = GSL.pow_7(-3.14)
    y_expected = GSL.pow(-3.14, 7.0)
    assert_rel y, y_expected, 1e-15, 'GSL.pow_7(-3.14)'

    y = GSL.pow_8(-3.14)
    y_expected = GSL.pow(-3.14, 8.0)
    assert_rel y, y_expected, 1e-15, 'GSL.pow_8(-3.14)'

    y = GSL.pow_9(-3.14)
    y_expected = GSL.pow(-3.14, 9.0)
    assert_rel y, y_expected, 1e-15, 'GSL.pow_9(-3.14)'

    -9.upto(9) { |n|
      y = GSL.pow_int(-3.14, n)
      y_expected = GSL.pow(-3.14, n)
      assert_rel y, y_expected, 1e-15, "GSL.pow_n(-3.14,#{n})"
    }
  end

  def test_ldexp
    y = GSL.ldexp(GSL::M_PI, -2)
    y_expected = GSL::M_PI_4
    assert_rel y, y_expected, 1e-15, 'GSL.ldexp(pi,-2)'

    y = GSL.ldexp(1.0, 2)
    y_expected = 4.000000
    assert_rel y, y_expected, 1e-15, 'GSL.ldexp(1.0,2)'

    y = GSL.ldexp(0.0, 2)
    y_expected = 0.0
    assert_rel y, y_expected, 1e-15, 'GSL.ldexp(0.0,2)'
  end

  def test_frexp
    y, e = GSL.frexp(GSL::M_PI)
    y_expected = GSL::M_PI_4
    e_expected = 2
    assert_rel y, y_expected, 1e-15, 'GSL.frexp(pi) fraction'
    assert_int e, e_expected, 'GSL.frexp(pi) exponent'

    y, e = GSL.frexp(2.0)
    y_expected = 0.5
    e_expected = 2
    assert_rel y, y_expected, 1e-15, 'GSL.frexp(2.0) fraction'
    assert_int e, e_expected, 'GSL.frexp(2.0) exponent'

    y, e = GSL.frexp(1.0 / 4.0)
    y_expected = 0.5
    e_expected = -1
    assert_rel y, y_expected, 1e-15, 'GSL.frexp(0.25) fraction'
    assert_int e, e_expected, 'GSL.frexp(0.25) exponent'

    y, e = GSL.frexp(1.0 / 4.0 - 4.0 * GSL::DBL_EPSILON)
    y_expected = 0.999999999999996447
    e_expected = -2
    assert_rel y, y_expected, 1e-15, 'GSL.frexp(0.25-eps) fraction'
    assert_int e, e_expected, 'GSL.frexp(0.25-eps) exponent'
  end

  def test_gsl
    x = GSL::M_PI
    y = 22.0 / 7.0

    10.times { |i|
      tol = GSL.pow(10, -i)
      res = GSL.fcmp(x, y, tol)
      assert_int res, i >= 4 ? -1 : 0, "GSL.fcmp(#{x},#{y},#{tol})"

      res = GSL.fcmp(y, x, tol)
      assert_int res, i >= 4 ? 1 : 0, "GSL.fcmp(#{y},#{x},#{tol})"
    }

    zero = 0.0
    one = 1.0
    inf = Math.exp(1.0e10)
    nan = inf / inf

    s = GSL.isinf(zero)
    assert_int s, 0, 'GSL.isinf(0)'

    s = GSL.isinf(one)
    assert_int s, 0, 'GSL.isinf(1)'

    s = GSL.isinf(inf)
    assert_int s, 1, 'GSL.isinf(inf)'

    # Commented out 2008/Oct/17 by YT
    # This test fails in (Darwin 9.5.0, gcc4.0.1):
    #  gsl_isinf() returns 1 for -inf
    #s = GSL.isinf(-inf)
    #assert_int s, -1, 'GSL.isinf(-inf)'

    s = GSL.isinf(nan)
    assert_int s, 0, 'GSL.isinf(nan)'

    s = GSL.isnan(zero)
    assert_int s, 0, 'GSL.isnan(0)'

    s = GSL.isnan(one)
    assert_int s, 0, 'GSL.isnan(1)'
    s = GSL.isnan(inf)
    assert_int s, 0, 'GSL.isnan(inf)'

    s = GSL.isnan(nan)
    assert_int s, 1, 'GSL.isnan(nan)'

    s = GSL.finite(zero)
    assert_int s, 1, 'GSL.finite(0)'

    s = GSL.finite(one)
    assert_int s, 1, 'GSL.finite(1)'

    s = GSL.finite(inf)
    assert_int s, 0, 'GSL.finite(inf)'

    s = GSL.finite(nan)
    assert_int s, 0, 'GSL.finite(nan)'
  end

end
