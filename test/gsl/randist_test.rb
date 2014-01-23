require 'test_helper'

class RandistTest < GSL::TestCase

  N = 100000
  MULTI_DIM = 10

  BINS = 100
  STEPS = 100

  R_GLOBAL = GSL::Rng.alloc

  def test_shuffle
    n = 10

    count = GSL::Matrix.calloc(n, n)
    x = GSL::Permutation.alloc(n)

    N.times { |i|
      n.times { |j| x[j] = j }
      GSL::Ran.shuffle(R_GLOBAL, x)
      n.times { |j| count.set(x[j], j, count[x[j], j] + 1) }
    }

    expected = N / 10.0

    n.times { |i|
      n.times { |j|
        d = (count[i, j] - expected).abs

        refute d > 1 && d / Math.sqrt(expected) > 5,
          "gsl_ran_shuffle #{i},#{j} (#{count[i, j] / N} observed vs 0.1 expected)"
      }
    }
  end

  def _test_moments(name, arg, a, b, pp)
    count, expected = 0, pp * N

    N.times {
      r = R_GLOBAL.send(*[name, arg].compact)
      count += 1 if r < b && r > a
    }

    refute((count - expected).abs / Math.sqrt(expected) > 3,
      "#{name}(#{arg}) [#{a},#{b}] (#{count.to_f / N} observed vs #{pp} expected)")
  end

  def _test_pdf(name, *args)
    pdf = "#{name}_pdf"

    a, b = -5.0, 5.0
    dx = (b - a) / BINS

    status = status_i = 0

    count = GSL::Vector.calloc(BINS)
    pp = GSL::Vector.calloc(BINS)

    N.times { |i|
      r = R_GLOBAL.send(name, *args)

      if r < b && r > a
        j = ((r - a) / dx).to_i
        count[j] = count[j] + 1
      end
    }

    BINS.times { |i|
      x = a + i * dx
      x = 0.0 if x.abs < 1e-10

      sum = 0.0

      STEPS.times { |j| sum += GSL::Ran.send(pdf, x + j * dx / STEPS, *args) }

      pp[i] = 0.5 * (GSL::Ran.send(pdf, x, *args) + 2 * sum + GSL::Ran.send(pdf, x + dx - 1e-7, *args)) * dx / STEPS
    }

    BINS.times { |i|
      x = a + i * dx
      d = count[i] - N * pp[i]

      status_i = (pp[i] == 0 ? count[i] != 0 : d < 1 && d / Math.sqrt(N * pp[i]) > 5) ? 1 : 0
      status |= status_i

      refute status_i == 1,
        "#{name} [#{x},#{x + dx}) (#{count[i]}/#{N}=#{count.to_f / N} observed vs #{pp[i]} expected)"
    }

    assert status.zero?, "#{name}, sampling against pdf over range [#{a},#{b})"
  end

  def test_randist
    _test_moments(:ugaussian,      nil,  0.0,   100.0, 0.5)
    _test_moments(:ugaussian,      nil, -1.0,     1.0, 0.6826895)
    _test_moments(:ugaussian,      nil,  3.0,     3.5, 0.0011172689)
    _test_moments(:ugaussian_tail, 3.0,  3.0,     3.5, 0.0011172689 / 0.0013498981)
    _test_moments(:exponential,    2.0,  0.0,     1.0, 1 - Math.exp(-0.5))
    _test_moments(:cauchy,         2.0,  0.0, 10000.0, 0.5)

    _test_moments(:discrete, GSL::Ran::Discrete.alloc(GSL::Vector.alloc(0.59, 0.4, 0.01)), -0.5, 0.5, 0.59)
    _test_moments(:discrete, GSL::Ran::Discrete.alloc(GSL::Vector.alloc(0.59, 0.4, 0.01)),  0.5, 1.5, 0.40)
    _test_moments(:discrete, GSL::Ran::Discrete.alloc(GSL::Vector.alloc(0.59, 0.4, 0.01)),  1.5, 3.5, 0.01)

    _test_moments(:discrete, GSL::Ran::Discrete.alloc(GSL::Vector.alloc(1, 9, 3, 4, 5, 8, 6, 7, 2, 0)), -0.5,  0.5, 1.0 / 45.0)
    _test_moments(:discrete, GSL::Ran::Discrete.alloc(GSL::Vector.alloc(1, 9, 3, 4, 5, 8, 6, 7, 2, 0)),  8.5,  9.5, 0)

    _test_pdf(:beta,   2.0, 3.0)
    _test_pdf(:cauchy, 2.0)
    _test_pdf(:chisq,  2.0)

    _test_pdf(:exponential,    2.0)
    _test_pdf(:exppow,         3.7, 0.3)
    _test_pdf(:fdist,          3.0, 4.0)
    _test_pdf(:flat,           3.0, 4.0)
    _test_pdf(:gamma,          2.5, 2.17)
    _test_pdf(:gaussian,       3.0)
    _test_pdf(:ugaussian_tail, 0.1, 2.0)
  end

end
