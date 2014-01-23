require 'test_helper'

class MonteTest < GSL::TestCase

  DIM = 1

  def test_miser
    return unless GSL::Monte::Miser.method_defined?(:params_get)

    miser = GSL::Monte::Miser.alloc(DIM)
    params = miser.params_get

    params.estimate_frac = 99
    miser.params_set(params)
    assert_abs miser.estimate_frac, 99, 1e-5, 'miser_estimate_frac'

    params.min_calls = 9
    miser.params_set(params)
    assert_int miser.min_calls, 9, 'miser_min_calls'

    params.min_calls_per_bisection = 7
    miser.params_set(params)
    assert_int miser.min_calls_per_bisection, 7, 'miser_min_calls_per_bisection'

    params.alpha = 3
    miser.params_set(params)
    assert_abs miser.alpha, 3, 1e-5, 'miser_alpha'

    params.dither = 4
    miser.params_set(params)
    assert_abs miser.dither, 4, 1e-5, 'miser_dither'
  end

  def test_vegas
    return unless GSL::Monte::Vegas.method_defined?(:params_get)

    vegas = GSL::Monte::Vegas.alloc(DIM)
    params = vegas.params_get

    params.alpha = 1
    vegas.params_set(params)
    assert_abs vegas.alpha, 1, 1e-5, 'vegas_alpha'

    params.iterations = 4
    vegas.params_set(params)
    assert_int vegas.iterations, 4, 'vegas_iterations'

    params.stage = 3
    vegas.params_set(params)
    assert_int vegas.stage, 3, 'vegas_stage'

    params.mode = GSL::Monte::Vegas::MODE_IMPORTANCE
    vegas.params_set(params)
    assert_int vegas.mode, GSL::Monte::Vegas::MODE_IMPORTANCE, 'vegas_mode MODE_IMPORTANCE'

    params.mode = GSL::Monte::Vegas::MODE_IMPORTANCE_ONLY
    vegas.params_set(params)
    assert_int vegas.mode, GSL::Monte::Vegas::MODE_IMPORTANCE_ONLY, 'vegas_mode MODE_IMPORTANCE_ONLY'

    params.mode = GSL::Monte::Vegas::MODE_STRATIFIED
    vegas.params_set(params)
    assert_int vegas.mode, GSL::Monte::Vegas::MODE_STRATIFIED, 'vegas_mode MODE_STRATIFIED'

    params.verbose = 0
    vegas.params_set(params)
    assert_int vegas.verbose, 0, 'vegas_verbose 0'

    params.verbose = 1
    vegas.params_set(params)
    assert_int vegas.verbose, 1, 'vegas_verbose 1'

    params.verbose = -1
    vegas.params_set(params)
    assert_int vegas.verbose, -1, 'vegas_verbose -1'
  end

end
