#!/usr/bin/env ruby
require("gsl")
require("../gsl_test2.rb")
include GSL::Test

exit unless GSL::Monte::Miser.method_defined?("params_get")

dim = 1
miser = GSL::Monte::Miser.alloc(dim)

params = miser.params_get

params.estimate_frac = 99
miser.params_set(params)
test_abs(miser.estimate_frac, 99, 1e-5, "miser_estimate_frac")

params.min_calls = 9
miser.params_set(params)
test_int(miser.min_calls, 9, "miser_min_calls")

params.min_calls_per_bisection = 7
miser.params_set(params)
test_int(miser.min_calls_per_bisection, 7, "miser_min_calls_per_bisection")

params.alpha = 3
miser.params_set(params)
test_abs(miser.alpha, 3, 1e-5, "miser_alpha")

params.dither = 4
miser.params_set(params)
test_abs(miser.dither, 4, 1e-5, "miser_dither")
