#!/usr/bin/env ruby
# Ruby/GSL implementation of GSL "randist/test.c"
require("gsl")
require("./gsl_test2.rb")
include GSL::Test
include Math
include GSL::Ran

N = 100000
MULTI_DIM = 10

GSL::IEEE::env_setup()
GSL::Rng::env_setup()

R_global = GSL::Rng.alloc()

def test_shuffle()
  n = 10
  status = 0
  count = GSL::Matrix.calloc(n, n)

  x = GSL::Permutation.alloc(n)
  for i in 0...N
    for j in 0...n
      x[j] = j
    end
    GSL::Ran.shuffle(R_global, x)
    for j in 0...n
      count.set(x[j], j, count[x[j],j]+1)
    end
  end

  for i in 0...n
    for j in 0...n
      expected = N/10.0
      d = (count[i,j] - expected).abs
      sigma = d/sqrt(expected)
      if sigma > 5 and d > 1
        status = 1
        GSL::Test::test(status, 
                        "gsl_ran_shuffle #{i},#{j} (#{count[i,j]/N} observed vs 0.1 expected)")
      end
    end
  end
  test(status, "gsl_ran_shuffle on {0, 1, 2, 3, 4, 5, 6, 7, 8, 9}")
end

def testMoments(name, a, b, pp)
  count = 0
  N.times do
    r = eval("R_global.#{name}")
    if r < b and r > a; count += 1; end
  end
  expected = pp*N
  sigma = (count - expected).abs/sqrt(expected)
  status = (sigma > 3) ? 1 : 0
  GSL::Test::test(status, "#{name} [#{a},#{b}] (#{count.to_f/N} observed vs #{pp} expected)")
end

BINS = 100
STEPS = 100
def testPDF(name, args)
  a = -5.0
  b = +5.0
  dx = (b - a)/BINS
  status = 0
  status_i = 0
  count = GSL::Vector.calloc(BINS)
  pp = GSL::Vector.calloc(BINS)

  for i in 0...N
    r = eval("R_global.#{name}(#{args})")
    if r < b and r > a
      j = ((r - a)/dx).to_i
      count[j] = count[j] + 1;
    end
  end
  for i in 0...BINS
    x = a + i*dx
    sum = 0.0
    if x.abs < 1e-10; x = 0.0; end
    for j in 1...STEPS
      sum += eval("GSL::Ran::#{name}_pdf(#{x+j*dx/STEPS}, #{args})")
    end
    pp[i] = 0.5*(eval("GSL::Ran::#{name}_pdf(#{x}, #{args})") + 2*sum + eval("GSL::Ran::#{name}_pdf(#{x+dx-1e-7}, #{args})"))*dx/STEPS
  end
  for i in 0...BINS
    x = a + i*dx
    d = (count[i] - N*pp[i])
    if pp[i] != 0
      s = d/sqrt(N*pp[i])
      status_i = ((s > 5) && (d > 1)) ? 1 : 0
    else
      status_i = (count[i] != 0) ? 1 : 0
    end
    status |= status_i
    if status_i == 1
      GSL::Test::test(status_i, "#{name} [#{x},#{x+dx}) (#{count[i]}/#{N}=#{count.to_f/N} observed vs #{pp[i]} expected)")
    end
  end
  if status == 0
    GSL::Test::test(status, "#{name}, sampling against pdf over range [#{a},#{b})")
  end
end

testMoments("ugaussian", 0.0, 100.0, 0.5)
testMoments("ugaussian", -1.0, 1.0, 0.6826895);
testMoments("ugaussian", 3.0, 3.5, 0.0011172689);
testMoments("ugaussian_tail(3.0)", 3.0, 3.5, 0.0011172689 / 0.0013498981);
testMoments("exponential(2.0)", 0.0, 1.0, 1 - exp(-0.5));
testMoments("cauchy(2.0)", 0.0, 10000.0, 0.5);

testMoments("discrete(GSL::Ran::Discrete.alloc(GSL::Vector.alloc(0.59, 0.4, 0.01)))", -0.5, 0.5, 0.59);
testMoments("discrete(GSL::Ran::Discrete.alloc(GSL::Vector.alloc(0.59, 0.4, 0.01)))", 0.5, 1.5, 0.40);
testMoments("discrete(GSL::Ran::Discrete.alloc(GSL::Vector.alloc(0.59, 0.4, 0.01)))", 1.5, 3.5, 0.01);

testMoments("discrete(GSL::Ran::Discrete.alloc(GSL::Vector.alloc(1, 9, 3, 4, 5, 8, 6, 7, 2, 0)))", -0.5,  0.5, 1.0/45.0 );
testMoments("discrete(GSL::Ran::Discrete.alloc(GSL::Vector.alloc(1, 9, 3, 4, 5, 8, 6, 7, 2, 0)))",  8.5,  9.5, 0 );

testPDF("beta", "2.0, 3.0")
testPDF("cauchy", "2.0")
testPDF("chisq", "2.0")

testPDF("exponential", "2.0")
testPDF("exppow", "3.7, 0.3")
testPDF("fdist", "3.0, 4.0")
testPDF("flat", "3.0, 4.0")
testPDF("gamma", "2.5, 2.17")
testPDF("gaussian", "3.0")
testPDF("ugaussian_tail", "0.1, 2.0")




