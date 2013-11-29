#!/usr/bin/env ruby
# Ruby/GSL implementation of GSL "rng/test.c"
require("gsl")
require("./gsl_test2.rb")
include GSL::Test
include Math

N = 1000
N2 = 20000

GSL::IEEE::env_setup()
GSL::Rng::env_setup()

def rng_test(type, seed, n, result)
  r, k = GSL::Rng.alloc(type), nil
  if seed != 0; r.set(seed); end

  n.times do
    k = r.get
  end

  status = (k != result) ? 1 : 0
  GSL::Test::test(status, "#{r.name}, #{n} steps (#{k} observed vs #{result} expected)")
end

def rng_float_test(type)
  ri = GSL::Rng.alloc(type)
  rf = GSL::Rng.alloc(type)
  k = 0
  status = 0
  begin
    k = ri.get
    u = rf.get
  end while k == 0

  c = k/u
  N2.times do
    k = ri.get
    u = rf.get
    if (c*k != u) 
      status = 1
      break
    end
  end
  GSL::Test::test(status, "#{ri.name}, ratio of int to double (#{c} observed vs #{k/u} expected)")
end

def rng_state_test(type)
  r = GSL::Rng.alloc(type)
  r_save = GSL::Rng.alloc(type)
  for i in 0...N; r.get; end
  GSL::Rng.memcpy(r_save, r)
  test_a = GSL::Vector.alloc(N)
  test_b = GSL::Vector.alloc(N)
  for i in 0...N; test_a[i] = r.get; end
  GSL::Rng.memcpy(r, r_save)
  for i in 0...N; test_b[i] = r.get; end
  status = 0
  for i in 0...N
    status |= (test_b[i] != test_a[i]) ? 1 : 0
  end
  GSL::Test::test(status, "#{r.name}, random number state consistency")
end

def rng_parallel_state_test(type)
  r1 = GSL::Rng.alloc(type)
  r2 = GSL::Rng.alloc(type)
  test_a = GSL::Vector.alloc(N)
  test_b = GSL::Vector.alloc(N)
  test_c = GSL::Vector.alloc(N)
  test_d = GSL::Vector.alloc(N)
  test_e = GSL::Vector.alloc(N)
  test_f = GSL::Vector.alloc(N)
  for i in 0...N; r1.get; end
  GSL::Rng.memcpy(r2, r1)
  for i in 0...N
    test_a[i] = r1.get
    test_b[i] = r2.get
    test_c[i] = r1.uniform_int(1234)
    test_d[i] = r2.uniform_int(1234)
    test_e[i] = r1.uniform
    test_f[i] = r2.uniform
  end
  status = 0
  for i in 0...N
    status |= (test_b[i] != test_a[i]) ? 1 : 0
    status |= (test_c[i] != test_d[i]) ? 1 : 0
    status |= (test_e[i] != test_f[i]) ? 1 : 0
  end
  GSL::Test::test(status, "#{r1.name}, parallel random number state consistency")
end

def rng_read_write_test(type)
  r = GSL::Rng.alloc(type)
  test_a = GSL::Vector.alloc(N)
  test_b = GSL::Vector.alloc(N)
  for i in 0...N; r.get; end
  r.fwrite("test.dat")
  for i in 0...N; test_a[i] = r.get; end
  r.fread("test.dat")
  for i in 0...N; test_b[i] = r.get; end
  status = 0
  for i in 0...N
    status |= (test_b[i] != test_a[i]) ? 1 : 0
  end
  GSL::Test::test(status, "#{r.name}, random number generator read and write")
  if FileTest.exist?("test.dat")
    File.delete("test.dat")
  end
end

def rng_max_test(r, ran_max)
  max = 0
  N2.times do
    k = r.get
    if k > max; max = k; end
  end
  actual_uncovered = ran_max - max;
  expect_uncovered = ran_max.to_f / (N2.to_f);
  status = ((max > ran_max) || (actual_uncovered > 7 * expect_uncovered)) ? 1 : 0
  return max, status;
end

def rng_min_test(r, ran_min, ran_max)
  min = 1000000000
  N2.times do
    k = r.get
    if k < min; min = k; end
  end
  actual_uncovered = min - ran_min
  expect_uncovered = ran_max.to_f/(N2.to_f)
  status = ((min < ran_min) || (actual_uncovered > 7 * expect_uncovered)) ? 1 : 0
  return min, status
end

def rng_sum_test(r)
  sum = 0.0
  N2.times do
    x = r.uniform - 0.5
    sum += x
  end
  sum = sum / (N2.to_f)
  sigma = sum*sqrt(12.0*N2)
  status = (sigma.abs > 3 || sigma.abs < 0.003) ? 1 : 0
  if status == 1
    fprintf(STDERR, "sum=%g, sigma=%g\n", sum. sigma)
  end
  return sigma, status
end

BINS = 17
EXTRA = 10

def rng_bin_test(r)
  count = GSL::Vector.calloc(BINS+EXTRA)
  chisq = 0.0
  for i in 0...N2
    j = r.uniform_int(BINS)
    count[j] = count[j] + 1
  end
  for i in 0...BINS
    x = (N2.to_f)/(BINS)
    d = count[i] - x
    chisq += (d*d)/(x)
  end
  sigma = sqrt(chisq/(BINS))
  status = (sigma.abs > 3 || sigma.abs < 0.003) ? 1 : 0
  for i in BINS...EXTRA
    if count[i] != 0
      status = 1
      GSL::Test::test(status, "#{r.name}, wrote outside range in bin test (#{i} observed vs #{BINS-1} expected)")
    end
  end
  return sigma, status
end

def generic_rng_test(type)
  r = GSL::Rng.alloc(type)
  name = r.name
  kmax = 0
  kmin = 1000
  sigma = 0.0
  ran_max = r.max
  ran_min = r.min
  kmax, status = rng_max_test(r, ran_max)
  GSL::Test::test(status, "#{name}, observed vs theoretical maximum (#{kmax} vs #{ran_max})")
  kmin, status = rng_min_test(r, ran_min, ran_max)
  GSL::Test::test(status, "#{name}, observed vs theoretical minimum (#{kmin} vs #{ran_min})")
  sigma, status = rng_sum_test(r)
  GSL::Test::test(status, "#{r.name}, sum test within acceptable sigma (observed #{sigma} sigma)")
  sigma, status = rng_bin_test(r)
  GSL::Test::test(status, "#{r.name}, bin test within acceptable chisq (observed #{sigma} sigma)")
  r.set(1)
  kmax, status = rng_max_test(r, ran_max)
  r.set(1)
  kmin, s = rng_min_test(r, ran_min, ran_max)
  status |= s
  r.set(1)
  sigma, s = rng_sum_test(r)
  status |= s
  r.set(12345)
  kmax, s = rng_max_test(r, ran_max)
  status |= s
  r.set(12345)
  kmin, s = rng_min_test(r, ran_min, ran_max)
  status |= s
  r.set(12345)
  sigma, s = rng_sum_test(r)
  status |= s
  GSL::Test::test(status, "#{r.name}, maximum and sum tests for non-default seeds")
end

rng_test("rand", 1, 10000, 1910041713);
rng_test("randu", 1, 10000, 1623524161);
rng_test("cmrg", 1, 10000, 719452880);
rng_test("minstd", 1, 10000, 1043618065);
rng_test("mrg", 1, 10000, 2064828650);
rng_test("taus", 1, 10000, 2733957125);
rng_test("taus113", 1, 1000, 1925420673);
rng_test("transputer", 1, 10000, 1244127297);
rng_test("vax", 1, 10000, 3051034865);

rng_test("borosh13", 1, 10000, 2513433025);
rng_test("fishman18", 1, 10000, 330402013);
rng_test("fishman2x", 1, 10000, 540133597);
rng_test("knuthran2", 1, 10000, 1084477620);
rng_test("knuthran", 310952, 1009 * 2009 + 1, 461390032);
rng_test("lecuyer21", 1, 10000, 2006618587);

rng_test("waterman14", 1, 10000, 3776680385);
rng_test("coveyou", 6, 10000, 1416754246);
rng_test("fishman20", 6, 10000, 248127575);
rng_test("ranlux", 314159265, 10000, 12077992);
rng_test("ranlux389", 314159265, 10000, 165942);
rng_test("ranlxs0", 1, 10000, 11904320);

rng_test("ranlxs1", 1, 10000, 8734328);
rng_test("ranlxs2", 1, 10000, 6843140); 
rng_test("ranlxd1", 1, 10000, 1998227290);
rng_test("ranlxd2", 1, 10000, 3949287736);
rng_test("slatec", 1, 10000, 45776);
rng_test("uni", 1, 10000, 9214);
rng_test("uni32", 1, 10000, 1155229825);
rng_test("zuf", 1, 10000, 3970);

rng_test("r250", 1, 10000, 1100653588);
rng_test("mt19937", 4357, 1000, 1186927261);
rng_test("mt19937_1999", 4357, 1000, 1030650439);
rng_test("mt19937_1998", 4357, 1000, 1309179303);
rng_test("tt800", 0, 10000, 2856609219);

rng_test("ran0", 0, 10000, 1115320064);
rng_test("ran1", 0, 10000, 1491066076);
rng_test("ran2", 0, 10000, 1701364455);
rng_test("ran3", 0, 10000, 186340785);

rng_test("ranmar", 1, 10000, 14428370);

rng_test("rand48", 0, 10000, 0xDE095043);
rng_test("rand48", 1, 10000, 0xEDA54977);

rng_test("random_glibc2", 0, 10000, 1908609430);
rng_test("random8_glibc2", 0, 10000, 1910041713);
rng_test("random32_glibc2", 0, 10000, 1587395585);
rng_test("random64_glibc2", 0, 10000, 52848624);
rng_test("random128_glibc2", 0, 10000, 1908609430);
rng_test("random256_glibc2", 0, 10000, 179943260);

rng_test("random_bsd", 0, 10000, 1457025928);
rng_test("random8_bsd", 0, 10000, 1910041713);
rng_test("random32_bsd", 0, 10000, 1663114331);
rng_test("random64_bsd", 0, 10000, 864469165);
rng_test("random128_bsd", 0, 10000, 1457025928);
rng_test("random256_bsd", 0, 10000, 1216357476);

rng_test("random_libc5", 0, 10000, 428084942);
rng_test("random8_libc5", 0, 10000, 1910041713);
rng_test("random32_libc5", 0, 10000, 1967452027);
rng_test("random64_libc5", 0, 10000, 2106639801);
rng_test("random128_libc5", 0, 10000, 428084942);
rng_test("random256_libc5", 0, 10000, 116367984);

rng_test("ranf", 0, 10000, 2152890433);
rng_test("ranf", 2, 10000, 339327233);

Rng.types.each do |type|
  rng_float_test(type)
end


Rng.types.each do |type|
  rng_state_test(type)
end

Rng.types.each do |type|
  rng_parallel_state_test(type)
end

Rng.types.each do |type|
  rng_read_write_test(type)
end

Rng.types.each do |type|
  generic_rng_test(type)
end
