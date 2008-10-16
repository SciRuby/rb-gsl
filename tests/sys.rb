#!/usr/bin/env ruby
# Ruby/GSL implementation of GSL "sum/test.c"
require("gsl")
require("./gsl_test2.rb")
include GSL::Test

GSL::IEEE::env_setup()

y = GSL::expm1(0.0);
y_expected = 0.0;
GSL::Test::test_rel(y, y_expected, 1e-15, "GSL::expm1(0.0)");

y = GSL::expm1(1e-10);
y_expected = 1.000000000050000000002e-10;
GSL::Test::test_rel(y, y_expected, 1e-15, "GSL::expm1(1e-10)");

y = GSL::expm1(-1e-10);
y_expected = -9.999999999500000000017e-11;
GSL::Test::test_rel(y, y_expected, 1e-15, "GSL::expm1(-1e-10)");

y = GSL::expm1(0.1);
y_expected = 0.1051709180756476248117078264902;
GSL::Test::test_rel(y, y_expected, 1e-15, "GSL::expm1(0.1)");

y = GSL::expm1(-0.1);
y_expected = -0.09516258196404042683575094055356;
GSL::Test::test_rel(y, y_expected, 1e-15, "GSL::expm1(-0.1)");

y = GSL::expm1(10.0);
y_expected = 22025.465794806716516957900645284;
GSL::Test::test_rel(y, y_expected, 1e-15, "GSL::expm1(10.0)");

y = GSL::expm1(-10.0);
y_expected = -0.99995460007023751514846440848444;
GSL::Test::test_rel(y, y_expected, 1e-15, "GSL::expm1(-10.0)");

# Test for log1p 

y = GSL::log1p(0.0);
y_expected = 0.0;
GSL::Test::test_rel(y, y_expected, 1e-15, "GSL::log1p(0.0)");

y = GSL::log1p(1e-10);
y_expected = 9.9999999995000000000333333333308e-11;
GSL::Test::test_rel(y, y_expected, 1e-15, "GSL::log1p(1e-10)");

y = GSL::log1p(0.1);
y_expected = 0.095310179804324860043952123280765;
GSL::Test::test_rel(y, y_expected, 1e-15, "GSL::log1p(0.1)");

y = GSL::log1p(10.0);
y_expected = 2.3978952727983705440619435779651;
GSL::Test::test_rel(y, y_expected, 1e-15, "GSL::log1p(10.0)");

# Test for GSL::hypot 

y = GSL::hypot(0.0, 0.0);
y_expected = 0.0;
GSL::Test::test_rel(y, y_expected, 1e-15, "GSL::hypot(0.0, 0.0)");

y = GSL::hypot(1e-10, 1e-10);
y_expected = 1.414213562373095048801688e-10;
GSL::Test::test_rel(y, y_expected, 1e-15, "GSL::hypot(1e-10, 1e-10)");

y = GSL::hypot(1e-38, 1e-38);
y_expected = 1.414213562373095048801688e-38;
GSL::Test::test_rel(y, y_expected, 1e-15, "GSL::hypot(1e-38, 1e-38)");

y = GSL::hypot(1e-10, -1.0);
y_expected = 1.000000000000000000005;
GSL::Test::test_rel(y, y_expected, 1e-15, "GSL::hypot(1e-10, -1)");

y = GSL::hypot(-1.0, 1e-10);
y_expected = 1.000000000000000000005;
GSL::Test::test_rel(y, y_expected, 1e-15, "GSL::hypot(-1, 1e-10)");

#y = GSL::hypot(1e307, 1e301);
#y_expected = 1.000000000000499999999999e307;
#GSL::Test::test_rel(y, y_expected, 1e-15, "GSL::hypot(1e307, 1e301)");

#y = GSL::hypot(1e301, 1e307);
#y_expected = 1.000000000000499999999999e307;
#GSL::Test::test_rel(y, y_expected, 1e-15, "GSL::hypot(1e301, 1e307)");

#y = GSL::hypot(1e307, 1e307);
#y_expected = 1.414213562373095048801688e307;
#GSL::Test::test_rel(y, y_expected, 1e-15, "GSL::hypot(1e307, 1e307)");

y = GSL::acosh(1.0);
y_expected = 0.0;
GSL::Test::test_rel(y, y_expected, 1e-15, "GSL::acosh(1.0)");

y = GSL::acosh(1.1);
y_expected = 4.435682543851151891329110663525e-1;
GSL::Test::test_rel(y, y_expected, 1e-15, "GSL::acosh(1.1)");

y = GSL::acosh(10.0);
y_expected = 2.9932228461263808979126677137742e0;
GSL::Test::test_rel(y, y_expected, 1e-15, "GSL::acosh(10.0)");

y = GSL::acosh(1e10);
y_expected = 2.3718998110500402149594646668302e1;
GSL::Test::test_rel(y, y_expected, 1e-15, "GSL::acosh(1e10)");

# Test for asinh 

y = GSL::asinh(0.0);
y_expected = 0.0;
GSL::Test::test_rel(y, y_expected, 1e-15, "GSL::asinh(0.0)");

y = GSL::asinh(1e-10);
y_expected = 9.9999999999999999999833333333346e-11;
GSL::Test::test_rel(y, y_expected, 1e-15, "GSL::asinh(1e-10)");

y = GSL::asinh(-1e-10);
y_expected = -9.9999999999999999999833333333346e-11;
GSL::Test::test_rel(y, y_expected, 1e-15, "GSL::asinh(1e-10)");

y = GSL::asinh(0.1);
y_expected = 9.983407889920756332730312470477e-2;
GSL::Test::test_rel(y, y_expected, 1e-15, "GSL::asinh(0.1)");

y = GSL::asinh(-0.1);
y_expected = -9.983407889920756332730312470477e-2;
GSL::Test::test_rel(y, y_expected, 1e-15, "GSL::asinh(-0.1)");

y = GSL::asinh(1.0);
y_expected = 8.8137358701954302523260932497979e-1;
GSL::Test::test_rel(y, y_expected, 1e-15, "GSL::asinh(1.0)");

y = GSL::asinh(-1.0);
y_expected = -8.8137358701954302523260932497979e-1;
GSL::Test::test_rel(y, y_expected, 1e-15, "GSL::asinh(-1.0)");

y = GSL::asinh(10.0);
y_expected = 2.9982229502979697388465955375965e0;
GSL::Test::test_rel(y, y_expected, 1e-15, "GSL::asinh(10)");

y = GSL::asinh(-10.0);
y_expected = -2.9982229502979697388465955375965e0;
GSL::Test::test_rel(y, y_expected, 1e-15, "GSL::asinh(-10)");

y = GSL::asinh(1e10);
y_expected = 2.3718998110500402149599646668302e1;
GSL::Test::test_rel(y, y_expected, 1e-15, "GSL::asinh(1e10)");

y = GSL::asinh(-1e10);
y_expected = -2.3718998110500402149599646668302e1;
GSL::Test::test_rel(y, y_expected, 1e-15, "GSL::asinh(-1e10)");

y = GSL::atanh(0.0);
y_expected = 0.0;
GSL::Test::test_rel(y, y_expected, 1e-15, "GSL::atanh(0.0)");

y = GSL::atanh(1e-20);
y_expected = 1e-20;
GSL::Test::test_rel(y, y_expected, 1e-15, "GSL::atanh(1e-20)");

y = GSL::atanh(-1e-20);
y_expected = -1e-20;
GSL::Test::test_rel(y, y_expected, 1e-15, "GSL::atanh(-1e-20)");

y = GSL::atanh(0.1);
y_expected = 1.0033534773107558063572655206004e-1;
GSL::Test::test_rel(y, y_expected, 1e-15, "GSL::atanh(0.1)");

y = GSL::atanh(-0.1);
y_expected = -1.0033534773107558063572655206004e-1;
GSL::Test::test_rel(y, y_expected, 1e-15, "GSL::atanh(-0.1)");

y = GSL::atanh(0.9);
y_expected = 1.4722194895832202300045137159439e0;
GSL::Test::test_rel(y, y_expected, 1e-15, "GSL::atanh(0.9)");

y = GSL::atanh(-0.9);
y_expected = -1.4722194895832202300045137159439e0;
GSL::Test::test_rel(y, y_expected, 1e-15, "GSL::atanh(0.9)");

# Test for pow_int 

y = GSL::pow_2(-3.14);
y_expected = pow(-3.14, 2.0);
GSL::Test::test_rel(y, y_expected, 1e-15, "GSL::pow_2(-3.14)");

y = GSL::pow_3(-3.14);
y_expected = pow(-3.14, 3.0);
GSL::Test::test_rel(y, y_expected, 1e-15, "GSL::pow_3(-3.14)");

y = GSL::pow_4(-3.14);
y_expected = pow(-3.14, 4.0);
GSL::Test::test_rel(y, y_expected, 1e-15, "GSL::pow_4(-3.14)");

y = GSL::pow_5(-3.14);
y_expected = pow(-3.14, 5.0);
GSL::Test::test_rel(y, y_expected, 1e-15, "GSL::pow_5(-3.14)");

y = GSL::pow_6(-3.14);
y_expected = pow(-3.14, 6.0);
GSL::Test::test_rel(y, y_expected, 1e-15, "GSL::pow_6(-3.14)");

y = GSL::pow_7(-3.14);
y_expected = pow(-3.14, 7.0);
GSL::Test::test_rel(y, y_expected, 1e-15, "GSL::pow_7(-3.14)");

y = GSL::pow_8(-3.14);
y_expected = pow(-3.14, 8.0);
GSL::Test::test_rel(y, y_expected, 1e-15, "GSL::pow_8(-3.14)");

y = GSL::pow_9(-3.14);
y_expected = pow(-3.14, 9.0);
GSL::Test::test_rel(y, y_expected, 1e-15, "GSL::pow_9(-3.14)");

for n in -9...10
  y = GSL::pow_int(-3.14, n);
  y_expected = pow(-3.14, n);
  GSL::Test::test_rel(y, y_expected, 1e-15, "GSL::pow_n(-3.14,#{n})")
end
  
# Test for ldexp 

y = GSL::ldexp(M_PI, -2);
y_expected = M_PI_4;
GSL::Test::test_rel(y, y_expected, 1e-15, "GSL::ldexp(pi,-2)");

y = GSL::ldexp(1.0, 2);
y_expected = 4.000000;
GSL::Test::test_rel(y, y_expected, 1e-15, "GSL::ldexp(1.0,2)");

y = GSL::ldexp(0.0, 2);
y_expected = 0.0;
GSL::Test::test_rel(y, y_expected, 1e-15, "GSL::ldexp(0.0,2)");

# Test for frexp 

y, e = GSL::frexp(M_PI);
y_expected = M_PI_4;
e_expected = 2;
GSL::Test::test_rel(y, y_expected, 1e-15, "GSL::frexp(pi) fraction");
GSL::Test::test_int(e, e_expected, "GSL::frexp(pi) exponent");

y, e = GSL::frexp(2.0);
y_expected = 0.5;
e_expected = 2;
GSL::Test::test_rel(y, y_expected, 1e-15, "GSL::frexp(2.0) fraction");
GSL::Test::test_int(e, e_expected, "GSL::frexp(2.0) exponent");

y, e = GSL::frexp(1.0 / 4.0);
y_expected = 0.5;
e_expected = -1;
GSL::Test::test_rel(y, y_expected, 1e-15, "GSL::frexp(0.25) fraction");
GSL::Test::test_int(e, e_expected, "GSL::frexp(0.25) exponent");

y, e = GSL::frexp(1.0 / 4.0 - 4.0 * GSL::DBL_EPSILON);
y_expected = 0.999999999999996447;
e_expected = -2;
GSL::Test::test_rel(y, y_expected, 1e-15, "GSL::frexp(0.25-eps) fraction");
GSL::Test::test_int(e, e_expected, "GSL::frexp(0.25-eps) exponent");


x = M_PI;
y = 22.0 / 7.0;

# test the basic function 

for i in 0...10
  tol = pow(10, -i);
  res = GSL::fcmp(x, y, tol);
  GSL::Test::test_int(res, -((i >= 4) ? 1 : 0), "GSL::fcmp(#{x},#{y},#{tol})")
end

for i in 0...10
  tol = pow(10, -i);
  res = GSL::fcmp(y, x, tol);
  GSL::Test::test_int(res,(i >= 4) ? 1 : 0, "GSL::fcmp(#{y},#{x},#{tol})")
end

zero = 0.0;
one = 1.0;
inf = Math::exp(1.0e10)
nan = inf / inf;

s = GSL::isinf(zero);
GSL::Test::test_int(s, 0, "GSL::isinf(0)");

s = GSL::isinf(one);
GSL::Test::test_int(s, 0, "GSL::isinf(1)");

s = GSL::isinf(inf);
GSL::Test::test_int(s, 1, "GSL::isinf(inf)");

s = GSL::isinf(-inf);
GSL::Test::test_int(s, -1, "GSL::isinf(-inf)");

s = GSL::isinf(nan);
GSL::Test::test_int(s, 0, "GSL::isinf(nan)");


s = GSL::isnan(zero);
GSL::Test::test_int(s, 0, "GSL::isnan(0)");

s = GSL::isnan(one);
GSL::Test::test_int(s, 0, "GSL::isnan(1)");
s = GSL::isnan(inf);
GSL::Test::test_int(s, 0, "GSL::isnan(inf)");

s = GSL::isnan(nan);
GSL::Test::test_int(s, 1, "GSL::isnan(nan)");


s = GSL::finite(zero);
GSL::Test::test_int(s, 1, "GSL::finite(0)");

s = GSL::finite(one);
GSL::Test::test_int(s, 1, "GSL::finite(1)");

s = GSL::finite(inf);
GSL::Test::test_int(s, 0, "GSL::finite(inf)");

s = GSL::finite(nan);
GSL::Test::test_int(s, 0, "GSL::finite(nan)");
