#!/usr/bin/env ruby
require("gsl")

GSL::Rng.env_setup()

r = GSL::Rng.alloc(GSL::Rng::DEFAULT)

n = 19
dim = 3
X = GSL::Matrix.alloc(n, dim)
y = GSL::Vector.alloc(n)
w = GSL::Vector.alloc(n)

file0 = "data.dat"
file1 = "fit.dat"

a = 0.1
File.open(file0, "w") do |f|
  for i in 0...n
    y0 = Math::exp(a)
    sigma = 0.1*y0
    val = r.gaussian(sigma)
    X[i,0] = 1.0
    X[i,1] = a
    X[i,2] = a*a
    y[i] = y0 + val
    w[i] = 1.0/(sigma*sigma)
    f.printf("%g %g %g\n", a, y[i], sigma)
    a += 0.1
  end
end

c, cov, chisq, status = GSL::MultiFit.wlinear(X, w, y)

printf("# best fit: Y = %g + %g X + %g X^2\n", c[0], c[1], c[2])
printf("# covariance matrix:\n")
printf("[ %+.5e, %+.5e, %+.5e\n", cov[0,0], cov[0,1], cov[0,2])
printf("  %+.5e, %+.5e, %+.5e\n", cov[1,0], cov[1,1], cov[1,2])
printf("  %+.5e, %+.5e, %+.5e ]\n", cov[2,0], cov[2,1], cov[2,2])
printf("# chisq = %g\n", chisq)

str = sprintf("%4.3f", c[0])
if c[1] > 0.0
  str += sprintf("+ %4.3f*x", c[1].abs)
else
  str += sprintf("- %4.3f*x", c[1].abs)
end
if c[2] > 0.0
  str += sprintf("+ %4.3f*x*x", c[2].abs)
else
  str += sprintf("- %4.3f*x*x", c[2].abs)
end

func = GSL::Function.alloc { |x, params|
  c0 = params[0]; c1 = params[1]; c2 = params[2]
  c0 + x*(c1 + x*c2)
}

func.set_params(c)

File.open(file1, "w") do |f|
  x = 0
  while x < 2
    f.printf("%e %e\n", x, func.eval(x))
    x += 0.01
  end
end

system("graph -T X -C -g 3 -y 0 7 -L '#{str}' -I e -m -1 -S 4 data.dat -I a -m 2 -S 0 fit.dat")

File.delete(file0)
File.delete(file1)
