#!/usr/bin/env ruby
require("gsl")
require("../gsl_test2.rb")
include GSL::Test

exit unless GSL::Monte::Vegas.method_defined?("params_get")

dim = 1
vegas = GSL::Monte::Vegas.alloc(dim)

params = vegas.params_get

params.alpha = 1
vegas.params_set(params)
test_abs(vegas.alpha, 1, 1e-5, "vegas_alpha")

params.iterations = 4
vegas.params_set(params)
test_int(vegas.iterations, 4, "vegas_iterations")

params.stage = 3
vegas.params_set(params)
test_int(vegas.stage, 3, "vegas_stage")

params.mode = GSL::Monte::Vegas::MODE_IMPORTANCE
vegas.params_set(params)
test_int(vegas.mode, GSL::Monte::Vegas::MODE_IMPORTANCE, "vegas_mode MODE_IMPORTANCE")

params.mode = GSL::Monte::Vegas::MODE_IMPORTANCE_ONLY
vegas.params_set(params)
test_int(vegas.mode, GSL::Monte::Vegas::MODE_IMPORTANCE_ONLY, "vegas_mode MODE_IMPORTANCE_ONLY")

params.mode = GSL::Monte::Vegas::MODE_STRATIFIED
vegas.params_set(params)
test_int(vegas.mode, GSL::Monte::Vegas::MODE_STRATIFIED, "vegas_mode MODE_STRATIFIED")

params.verbose = 0
vegas.params_set(params)
test_int(vegas.verbose, 0, "vegas_verbose 0")
params.verbose = 1
vegas.params_set(params)
test_int(vegas.verbose, 1, "vegas_verbose 1")
params.verbose = -1
vegas.params_set(params)
test_int(vegas.verbose, -1, "vegas_verbose -1")
