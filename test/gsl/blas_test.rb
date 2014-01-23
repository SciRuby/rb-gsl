require 'test_helper'

class BlasTest < GSL::TestCase

  DBLEPS = 1e-6

  def test_amax
    v = GSL::Vector.alloc(0.537, 0.826)
    assert_int v.idamax, 1, 'damax'

    vz = GSL::Vector::Complex.alloc([0.913, -0.436], [-0.134, 0.129])
    assert_int vz.izamax, 0, 'zmax'
  end

  def test_asum
    v = GSL::Vector.alloc(0.271, -0.012)
    assert_rel v.dasum, 0.283, DBLEPS, 'dasum'

    vz = GSL::Vector::Complex.alloc([-0.046, -0.671], [-0.323, 0.785])
    assert_rel vz.dzasum, 1.825, DBLEPS, 'dzasum'
  end

  def test_axpy
    x = GSL::Vector.alloc(0.029)
    y = GSL::Vector.alloc(-0.992)
    e = GSL::Vector.alloc(-1.0007)

    y2 = GSL::Blas.daxpy(-0.3, x, y)

    assert_rel y2[0], e[0], DBLEPS, 'daxpy'

    x = GSL::Vector::Complex.alloc([[0.776, -0.671]])
    y = GSL::Vector::Complex.alloc([[0.39, 0.404]])
    e = GSL::Vector::Complex.alloc([[1.061, 1.18]])

    y2 = GSL::Blas.zaxpy(GSL::Complex.alloc(0, 1), x, y)

    assert_rel y2[0].re, e[0].re, DBLEPS, 'zaxpy real'
    assert_rel y2[0].im, e[0].im, DBLEPS, 'zaxpy imag'
  end

  def test_copy
    x = GSL::Vector.alloc(0.002)
    y = GSL::Vector.alloc(-0.921)
    e = GSL::Vector.alloc(0.002)

    GSL::Blas.dcopy(x, y)

    assert_rel y[0], e[0], DBLEPS, 'dcopy'

    x = GSL::Vector::Complex.alloc([[ 0.315, -0.324]])
    y = GSL::Vector::Complex.alloc([[-0.312, -0.748]])
    e = GSL::Vector::Complex.alloc([[0.315, -0.324]])

    GSL::Blas.zcopy(x, y)

    assert_rel y[0].re, e[0].re, DBLEPS, 'zcopy real'
    assert_rel y[0].im, e[0].im, DBLEPS, 'zcopy imag'
  end

  def test_dnrm2
    require 'narray'

    e = Math.sqrt((0..4).inject { |m, x| m + x * x })

    v = GSL::Vector.indgen(5)
    v_dnrm2 = GSL::Blas.dnrm2(v)

    assert_rel v_dnrm2, e, DBLEPS, 'GSL::Blas.dnrm2(GSL::Vector)'

    na = NArray.float(5).indgen!
    na_dnrm2 = GSL::Blas.dnrm2(na)

    assert_rel na_dnrm2, e, DBLEPS, 'GSL::Blas.dnrm2(NArray)'

    assert_rel na_dnrm2, v_dnrm2, 0, 'GSL::Blas.dnrm2(NArray) == GSL::Blas.dnrm2(GSL::Vector)'
  end

end
