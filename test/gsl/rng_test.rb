require 'test_helper'

class RngTest < GSL::TestCase

  N = 1000
  N2 = 20000

  BINS = 17
  EXTRA = 10

  {
    'rand'       => [[1, 10000, 1910041713]],
    'randu'      => [[1, 10000, 1623524161]],
    'cmrg'       => [[1, 10000, 719452880]],
    'minstd'     => [[1, 10000, 1043618065]],
    'mrg'        => [[1, 10000, 2064828650]],
    'taus'       => [[1, 10000, 2733957125]],
    'taus113'    => [[1, 1000, 1925420673]],
    'transputer' => [[1, 10000, 1244127297]],
    'vax'        => [[1, 10000, 3051034865]],

    'borosh13'  => [[1, 10000, 2513433025]],
    'fishman18' => [[1, 10000, 330402013]],
    'fishman2x' => [[1, 10000, 540133597]],
    'knuthran2' => [[1, 10000, 1084477620]],
    'knuthran'  => [[310952, 1009 * 2009 + 1, 461390032]],
    'lecuyer21' => [[1, 10000, 2006618587]],

    'waterman14' => [[1, 10000, 3776680385]],
    'coveyou'    => [[6, 10000, 1416754246]],
    'fishman20'  => [[6, 10000, 248127575]],
    'ranlux'     => [[314159265, 10000, 12077992]],
    'ranlux389'  => [[314159265, 10000, 165942]],
    'ranlxs0'    => [[1, 10000, 11904320]],

    'ranlxs1' => [[1, 10000, 8734328]],
    'ranlxs2' => [[1, 10000, 6843140]],
    'ranlxd1' => [[1, 10000, 1998227290]],
    'ranlxd2' => [[1, 10000, 3949287736]],
    'slatec'  => [[1, 10000, 45776]],
    'uni'     => [[1, 10000, 9214]],
    'uni32'   => [[1, 10000, 1155229825]],
    'zuf'     => [[1, 10000, 3970]],

    'r250'         => [[1, 10000, 1100653588]],
    'mt19937'      => [[4357, 1000, 1186927261]],
    'mt19937_1999' => [[4357, 1000, 1030650439]],
    'mt19937_1998' => [[4357, 1000, 1309179303]],
    'tt800'        => [[0, 10000, 2856609219]],

    'ran0' => [[0, 10000, 1115320064]],
    'ran1' => [[0, 10000, 1491066076]],
    'ran2' => [[0, 10000, 1701364455]],
    'ran3' => [[0, 10000, 186340785]],

    'ranmar' => [[1, 10000, 14428370]],

    'rand48' => [[0, 10000, 0xDE095043],
                 [1, 10000, 0xEDA54977]],

    'random_glibc2'    => [[0, 10000, 1908609430]],
    'random8_glibc2'   => [[0, 10000, 1910041713]],
    'random32_glibc2'  => [[0, 10000, 1587395585]],
    'random64_glibc2'  => [[0, 10000, 52848624]],
    'random128_glibc2' => [[0, 10000, 1908609430]],
    'random256_glibc2' => [[0, 10000, 179943260]],

    'random_bsd'    => [[0, 10000, 1457025928]],
    'random8_bsd'   => [[0, 10000, 1910041713]],
    'random32_bsd'  => [[0, 10000, 1663114331]],
    'random64_bsd'  => [[0, 10000, 864469165]],
    'random128_bsd' => [[0, 10000, 1457025928]],
    'random256_bsd' => [[0, 10000, 1216357476]],

    'random_libc5'    => [[0, 10000, 428084942]],
    'random8_libc5'   => [[0, 10000, 1910041713]],
    'random32_libc5'  => [[0, 10000, 1967452027]],
    'random64_libc5'  => [[0, 10000, 2106639801]],
    'random128_libc5' => [[0, 10000, 428084942]],
    'random256_libc5' => [[0, 10000, 116367984]],

    'ranf' => [[0, 10000, 2152890433],
               [2, 10000, 339327233]]
  }.each { |type, args|
    args.each_with_index { |(seed, n, result), i|
      define_method("test_#{type}_#{i}") {
        r, k = GSL::Rng.alloc(type), nil
        r.set(seed) if seed != 0

        n.times { k = r.get }

        assert k == result, "#{r.name}, #{n} steps (#{k} observed vs #{result} expected)"
      }
    }
  }

  GSL::Rng.types.each { |type|
    define_method("test_float_#{type}")          { _rng_float_test(type) }
    define_method("test_state_#{type}")          { _rng_state_test(type) }
    define_method("test_parallel_state_#{type}") { _rng_parallel_state_test(type) }
    define_method("test_read_write_#{type}")     { _rng_read_write_test(type) }
    define_method("test_generic_#{type}")        { _generic_rng_test(type) }
  }

  def _rng_float_test(type)
    ri = GSL::Rng.alloc(type)
    rf = GSL::Rng.alloc(type)

    status = k = 0

    begin
      k = ri.get
      u = rf.get
    end while k == 0

    c = k / u

    N2.times {
      k = ri.get
      u = rf.get

      if c * k != u
        status = 1
        break
      end
    }

    assert status.zero?, "#{ri.name}, ratio of int to double (#{c} observed vs #{k/u} expected)"
  end

  def _rng_state_test(type)
    r = GSL::Rng.alloc(type)
    r_save = GSL::Rng.alloc(type)

    N.times { r.get }

    GSL::Rng.memcpy(r_save, r)

    test_a = GSL::Vector.alloc(N)
    test_b = GSL::Vector.alloc(N)

    N.times { |i| test_a[i] = r.get }

    GSL::Rng.memcpy(r, r_save)
    N.times { |i| test_b[i] = r.get }

    assert((0...N).all? { |i| test_b[i] == test_a[i] }, "#{r.name}, random number state consistency")
  end

  def _rng_parallel_state_test(type)
    r1 = GSL::Rng.alloc(type)
    r2 = GSL::Rng.alloc(type)

    test_a = GSL::Vector.alloc(N)
    test_b = GSL::Vector.alloc(N)
    test_c = GSL::Vector.alloc(N)
    test_d = GSL::Vector.alloc(N)
    test_e = GSL::Vector.alloc(N)
    test_f = GSL::Vector.alloc(N)

    N.times { r1.get }

    GSL::Rng.memcpy(r2, r1)

    N.times { |i|
      test_a[i] = r1.get
      test_b[i] = r2.get
      test_c[i] = r1.uniform_int(1234)
      test_d[i] = r2.uniform_int(1234)
      test_e[i] = r1.uniform
      test_f[i] = r2.uniform
    }

    assert((0...N).all? { |i| test_a[i] == test_b[i] }, "#{r1.name}, parallel random number state consistency (a/b)")
    assert((0...N).all? { |i| test_c[i] == test_d[i] }, "#{r1.name}, parallel random number state consistency (c/d)")
    assert((0...N).all? { |i| test_e[i] == test_f[i] }, "#{r1.name}, parallel random number state consistency (e/f)")
  end

  def _rng_read_write_test(type)
    r = GSL::Rng.alloc(type)

    test_a = GSL::Vector.alloc(N)
    test_b = GSL::Vector.alloc(N)

    N.times { r.get }

    r.fwrite('test.dat')
    N.times { |i| test_a[i] = r.get }

    r.fread('test.dat')
    N.times { |i| test_b[i] = r.get }

    assert((0...N).all? { |i| test_b[i] == test_a[i] }, "#{r.name}, random number generator read and write")
  ensure
    File.delete('test.dat') if FileTest.exist?('test.dat')
  end

  def _generic_rng_test(type)
    r = GSL::Rng.alloc(type)

    name, ran_max, ran_min = r.name, r.max, r.min
    kmax, kmin, sigma = 0, 1000, 0.0

    kmax, status = _rng_max_test(r, ran_max)
    assert status.zero?, "#{name}, observed vs theoretical maximum (#{kmax} vs #{ran_max})"

    kmin, status = _rng_min_test(r, ran_min, ran_max)
    assert status.zero?, "#{name}, observed vs theoretical minimum (#{kmin} vs #{ran_min})"

    sigma, status = _rng_sum_test(r)
    assert status.zero?, "#{r.name}, sum test within acceptable sigma (observed #{sigma} sigma)"

    sigma, status = _rng_bin_test(r)
    assert status.zero?, "#{r.name}, bin test within acceptable chisq (observed #{sigma} sigma)"

    r.set(1)
    kmax, status = _rng_max_test(r, ran_max)

    r.set(1)
    kmin, s = _rng_min_test(r, ran_min, ran_max)
    status |= s

    r.set(1)
    sigma, s = _rng_sum_test(r)
    status |= s

    r.set(12345)
    kmax, s = _rng_max_test(r, ran_max)
    status |= s

    r.set(12345)
    kmin, s = _rng_min_test(r, ran_min, ran_max)
    status |= s

    r.set(12345)
    sigma, s = _rng_sum_test(r)
    status |= s

    assert status.zero?, "#{r.name}, maximum and sum tests for non-default seeds"
  end

  def _rng_max_test(r, ran_max)
    max = 0

    N2.times {
      k = r.get
      max = k if k > max
    }

    actual_uncovered = ran_max - max
    expect_uncovered = ran_max.to_f / N2.to_f

    [max, max > ran_max || actual_uncovered > 7 * expect_uncovered ? 1 : 0]
  end

  def _rng_min_test(r, ran_min, ran_max)
    min = 1000000000

    N2.times {
      k = r.get
      min = k if k < min
    }

    actual_uncovered = min - ran_min
    expect_uncovered = ran_max.to_f / N2.to_f

    [min, min < ran_min || actual_uncovered > 7 * expect_uncovered ? 1 : 0]
  end

  def _rng_sum_test(r)
    sum = 0.0

    N2.times {
      x = r.uniform - 0.5
      sum += x
    }

    sum /= N2.to_f
    sigma = sum * Math.sqrt(12.0 * N2)

    [sigma, sigma.abs > 3 || sigma.abs < 0.003 ? 1 : 0]
  end

  def _rng_bin_test(r)
    count = GSL::Vector.calloc(BINS + EXTRA)
    chisq = 0.0

    N2.times { count[r.uniform_int(BINS)] += 1 }

    BINS.times { |i|
      x = N2.to_f / BINS
      d = count[i] - x
      chisq += d * d / x
    }

    BINS.upto(EXTRA - 1) { |i|
      assert count[i].zero?, "#{r.name}, wrote outside range in bin test (#{i} observed vs #{BINS - 1} expected)"
    }

    [sigma = Math.sqrt(chisq / BINS), sigma.abs > 3 || sigma.abs < 0.003 ? 1 : 0]
  end

end
