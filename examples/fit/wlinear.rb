#!/usr/bin/env ruby
require("gsl")

n = 4
x = GSL::Vector[1970.0, 1980, 1990, 2000]
y = GSL::Vector[12.0, 11, 14, 13]
w = GSL::Vector[0.1, 0.2, 0.3, 0.4]

c0, c1, cov00, cov01, cov11, chisq = GSL::Fit::wlinear(x, w, y)

printf("# best fit: Y = %g + %g X\n", c0, c1);
printf("# covariance matrix:\n");
printf("# [ %g, %g\n#   %g, %g]\n",
        cov00, cov01, cov01, cov11);
printf("# chisq = %g\n", chisq);

File.open("data.dat", "w") do |f|
  for i in 0...n do
    f.printf("%e %e %e\n", x[i], y[i], 1.0/Math::sqrt(w[i]))
  end
end

begin
  ffit = File.open("fit.dat", "w")
  fhi = File.open("hi.dat", "w")
  flo = File.open("lo.dat", "w")
  for i in -30...130 do
    xf = x[0] + (i/100.0) * (x[n-1] - x[0])

    yf, yf_err = GSL::Fit::linear_est(xf, c0, c1, cov00, cov01, cov11)

    ffit.printf("%g %g\n", xf, yf)
    fhi.printf("%g %g\n", xf, yf + yf_err)
    flo.printf("%g %g\n", xf, yf - yf_err)
  end
ensure
  ffit.close
  fhi.close
  flo.close
end

system("graph -T X -C -g 3 -X x -Y y -x 1960 2010 -y 0 20 -m 0 -S 2 -Ie data.dat -S 0 -I a -m 1 fit.dat -m 2 hi.dat -m 2 lo.dat")
File.delete("data.dat")
File.delete("fit.dat")
File.delete("hi.dat")
File.delete("lo.dat")
