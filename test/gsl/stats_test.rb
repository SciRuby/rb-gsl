require 'test_helper'

class StatsTest < GSL::TestCase

  def test_stats
    lew = GSL::Vector.alloc(
      -213, -564,  -35,  -15,  141,  115, -420, -360,  203, -338, -431,  194,
      -220, -513,  154, -125, -559,   92,  -21, -579,  -52,   99, -543, -175,
       162, -457, -346,  204, -300, -474,  164, -107, -572,   -8,   83, -541,
      -224,  180, -420, -374,  201, -236, -531,   83,   27, -564, -112,  131,
      -507, -254,  199, -311, -495,  143,  -46, -579,  -90,  136, -472, -338,
       202, -287, -477,  169, -124, -568,   17,   48, -568, -135,  162, -430,
      -422,  172,  -74, -577,  -13,   92, -534, -243,  194, -355, -465,  156,
       -81, -578,  -64,  139, -449, -384,  193, -198, -538,  110,  -44, -577,
        -6,   66, -552, -164,  161, -460, -344,  205, -281, -504,  134,  -28,
      -576, -118,  156, -437, -381,  200, -220, -540,   83,   11, -568, -160,
       172, -414, -408,  188, -125, -572,  -32,  139, -492, -321,  205, -262,
      -504,  142,  -83, -574,    0,   48, -571, -106,  137, -501, -266,  190,
      -391, -406,  194, -186, -553,   83,  -13, -577,  -49,  103, -515, -280,
       201,  300, -506,  131,  -45, -578,  -80,  138, -462, -361,  201, -211,
      -554,   32,   74, -533, -235,  187, -372, -442,  182, -147, -566,   25,
        68, -535, -244,  194, -351, -463,  174, -125, -570,   15,   72, -550,
      -190,  172, -424, -385,  198, -218, -536,   96
    )

    mean = lew.mean
    sd   = lew.sd
    lag1 = lew.lag1_autocorrelation

    assert_rel mean, -177.435000000000,    1e-15, 'lew gsl_stats_mean'
    assert_rel sd,    277.332168044316,    1e-15, 'lew gsl_stats_sd'
    assert_rel lag1,   -0.307304800605679, 1e-14, 'lew autocorrelation'

    rel = 1e-10

    rawa = GSL::Vector.alloc(
      0.0421, 0.0941, 0.1064, 0.0242, 0.1331,
      0.0773, 0.0243, 0.0815, 0.1186, 0.0356,
      0.0728, 0.0999, 0.0614, 0.0479
    )

    rawb = GSL::Vector.alloc(
      0.1081, 0.0986, 0.1566, 0.1961, 0.1125,
      0.1942, 0.1079, 0.1021, 0.1583, 0.1673,
      0.1675, 0.1856, 0.1688, 0.1512
    )

    raww = GSL::Vector.alloc(
      0.0000, 0.0000, 0.0000, 3.000, 0.0000,
      1.0000, 1.0000, 1.0000, 0.000, 0.5000,
      7.0000, 5.0000, 4.0000, 0.123
    )

    assert_rel mean = rawa.mean, 0.0728, rel, 'gsl_stats_mean'

    assert_rel rawa.variance_with_fixed_mean(mean), 0.00113837428571429,
      rel, 'gsl_stats_variance_with_fixed_mean'

    assert_rel rawa.sd_with_fixed_mean(mean), 0.0337398026922845,
      rel, 'gsl_stats_sd_with_fixed_mean'

    assert_rel rawb.variance,         0.00124956615384615,  rel, 'gsl_stats_variance'
    assert_rel rawa.sd,               0.0350134479659107,   rel, 'gsl_stats_sd'
    assert_rel rawa.absdev,           0.0287571428571429,   rel, 'gsl_stats_absdev'
    assert_rel rawa.skew,             0.0954642051479004,   rel, 'gsl_stats_skew'
    assert_rel rawa.kurtosis,        -1.38583851548909,     rel, 'gsl_stats_kurtosis'
    assert_rel rawa.wmean(raww),      0.0678111523670601,   rel, 'gsl_stats_wmean'
    assert_rel rawa.wvariance(raww),  0.000769562962860317, rel, 'gsl_stats_wvariance'
    assert_rel rawa.wsd(raww),        0.0277409978706664,   rel, 'gsl_stats_wsd'
    assert_rel rawa.wabsdev(raww),    0.0193205027504008,   rel, 'gsl_stats_wabsdev'
    assert_rel rawa.wskew(raww),     -0.373631000307076,    rel, 'gsl_stats_wskew'
    assert_rel rawa.wkurtosis(raww), -1.48114233353963,     rel, 'gsl_stats_wkurtosis'

    assert_rel GSL::Stats.covariance(rawa, rawb), -0.000139021538461539, rel, 'gsl_stats_covariance'

    if GSL::GSL_VERSION >= '1.9.90'
      assert_rel GSL::Stats.correlation(rawa, rawb), -0.112322712666074171, rel, 'gsl_stats_correlation'
      assert_rel GSL::Stats.pvariance(rawa, rawb),    0.00123775384615385,  rel, 'gsl_stats_pvariance'
      assert_rel GSL::Stats.ttest(rawa, rawb),       -5.67026326985851,     rel, 'gsl_stats_ttest'
    end

    assert_rel rawa.max, 0.1331, rel, 'gsl_stats_max'
    assert_rel rawa.min, 0.0242, rel, 'gsl_stats_min'

    min, max = rawa.minmax
    assert_rel max, 0.1331, rel, 'gsl_stats_minmax: max'
    assert_rel min, 0.0242, rel, 'gsl_stats_minmax: min'

    assert rawa.max_index == 4, 'gsl_stats_max_index'
    assert rawa.min_index == 3, 'gsl_stats_min_index'

    min_index, max_index = rawa.minmax_index
    assert max_index == 4, 'gsl_stats_minmax_index: max'
    assert min_index == 3, 'gsl_stats_minmax_index: min'

    rawa.sort!

    assert_rel rawa.median_from_sorted_data, 0.07505,
      rel, 'gsl_stats_median_from_sorted_data'

    assert_rel rawa.subvector(0, rawa.size - 1).median_from_sorted_data, 0.0728,
      rel, 'gsl_stats_median_from_sorted_data'

    assert_rel rawa.quantile_from_sorted_data(0.0), 0.0242,
      rel, 'gsl_stats_quantile_from_sorted_data'

    assert_rel rawa.quantile_from_sorted_data(1.0), 0.1331,
      rel, 'gsl_stats_quantile_from_sorted_data (100)'

    assert_rel rawa.quantile_from_sorted_data(0.5), 0.07505,
      rel, 'gsl_stats_quantile_from_sorted_data (50even)'

    assert_rel rawa.subvector(0, rawa.size-1).quantile_from_sorted_data(0.5), 0.0728,
      rel, 'gsl_stats_quantile_from_sorted_data (50odd)'
  end

  def test_variance_with_fixed_mean
    v = GSL::Vector[1..8]
    assert_raises(ArgumentError, 'check for no args') { v.variance_with_fixed_mean }
  end

end
