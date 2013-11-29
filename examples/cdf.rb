#!/usr/bin/env ruby
require("gsl")

x = 2.0

P = GSL::Cdf::ugaussian_P(x);
printf("prob(x < %f) = %f\n", x, P);

Q = GSL::Cdf::ugaussian_Q(x);
printf("prob(x > %f) = %f\n", x, Q);

x = GSL::Cdf::ugaussian_Pinv(P);
printf("Pinv(%f) = %f\n", P, x);

x = GSL::Cdf::ugaussian_Qinv(Q);
printf("Qinv(%f) = %f\n", Q, x);
