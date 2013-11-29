#!/usr/bin/env ruby
require("gsl")
require("./gsl_test.rb")
include GSL::Test

lew = GSL::Vector.alloc(    -213, -564,  -35,  -15,  141,  115, -420, -360,  203, -338, -431,  194,
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
    -190,  172, -424, -385,  198, -218, -536,   96 )

expected_mean = -177.435000000000;
expected_sd = 277.332168044316;
expected_lag1 = -0.307304800605679;

mean = lew.mean()
sd = lew.sd()
lag1 = lew.lag1_autocorrelation()

GSL::Test::test_rel(mean, expected_mean, 1e-15, "lew gsl_stats_mean") ;
GSL::Test::test_rel(sd, expected_sd, 1e-15, "lew gsl_stats_sd") ;
GSL::Test::test_rel(lag1, expected_lag1, 1e-14, "lew autocorrelation") ;

rel = 1e-10

rawa = GSL::Vector.alloc(0.0421, 0.0941, 0.1064, 0.0242, 0.1331,0.0773, 0.0243, 0.0815, 0.1186, 0.0356,0.0728, 0.0999, 0.0614, 0.0479)

rawb = GSL::Vector.alloc(0.1081, 0.0986, 0.1566, 0.1961, 0.1125,
 0.1942, 0.1079, 0.1021, 0.1583, 0.1673,
 0.1675, 0.1856, 0.1688, 0.1512)
 
raww = GSL::Vector.alloc(0.0000, 0.0000, 0.0000, 3.000, 0.0000,
 1.000, 1.000, 1.000, 0.000, 0.5000,
 7.000, 5.000, 4.000, 0.123)
 
mean = rawa.mean
expected = 0.0728
GSL::Test::test_rel(mean, expected, rel, "gsl_stats_mean")

var = rawa.variance_with_fixed_mean(mean)
expected = 0.00113837428571429
GSL::Test::test_rel(var, expected, rel, "gsl_stats_variance_with_fixed_mean")

sd = rawa.sd_with_fixed_mean(mean)
expected = 0.0337398026922845
GSL::Test::test_rel(sd, expected, rel, "gsl_stats_sd_with_fixed_mean")


 
var = rawb.variance()
expected = 0.00124956615384615
GSL::Test::test_rel(var, expected, rel, "gsl_stats_variance")

sd = rawa.sd()
expected = 0.0350134479659107
GSL::Test::test_rel(sd, expected, rel, "gsl_stats_sd")

absdev = rawa.absdev()
expected = 0.0287571428571429
GSL::Test::test_rel(absdev, expected, rel, "gsl_stats_absdev")

skew = rawa.skew()
expected = 0.0954642051479004
GSL::Test::test_rel(skew, expected, rel, "gsl_stats_skew")

kurtosis = rawa.kurtosis()
expected = -1.38583851548909
GSL::Test::test_rel(kurtosis, expected, rel, "gsl_stats_kurtosis")

wmean = rawa.wmean(raww)
expected = 0.0678111523670601
GSL::Test::test_rel(wmean, expected, rel, "gsl_stats_wmean")

wvariance = rawa.wvariance(raww)
expected = 0.000769562962860317
GSL::Test::test_rel(wvariance, expected, rel, "gsl_stats_wvariance")

wsd = rawa.wsd(raww)
expected = 0.0277409978706664
GSL::Test::test_rel(wsd, expected, rel, "gsl_stats_wsd")

wabsdev = rawa.wabsdev(raww)
expected = 0.0193205027504008
GSL::Test::test_rel(wabsdev, expected, rel, "gsl_stats_wabsdev")

wskew = rawa.wskew(raww)
expected = -0.373631000307076
GSL::Test::test_rel(wskew, expected, rel, "gsl_stats_wskew")

wkurtosis = rawa.wkurtosis(raww)
expected =  -1.48114233353963
GSL::Test::test_rel(wkurtosis, expected, rel, "gsl_stats_wkurtosis")

c = GSL::Stats::covariance(rawa, rawb)
expected = -0.000139021538461539
GSL::Test::test_rel(c, expected, rel, "gsl_stats_covariance")

if GSL_VERSION >= "1.9.90"
	r = GSL::Stats::correlation(rawa, rawb)
	expected = -0.112322712666074171
	GSL::Test::test_rel(r, expected, rel, "gsl_stats_correlation")
	
	pv = GSL::Stats::pvariance(rawa, rawb)
	expected = 0.00123775384615385
	GSL::Test::test_rel(pv, expected, rel, "gsl_stats_pvariance")	
	
	t = GSL::Stats::ttest(rawa, rawb)
	expected =  -5.67026326985851
	GSL::Test::test_rel(t, expected, rel, "gsl_stats_ttest")		
end

expected = 0.1331
GSL::Test::test_rel(rawa.max, expected, rel, "gsl_stats_max")

expected = 0.0242
GSL::Test::test_rel(rawa.min, expected, rel, "gsl_stats_min")
min, max = rawa.minmax
expected = 0.1331
GSL::Test::test_rel(max, expected, rel, "gsl_stats_minmax: max")
expected = 0.0242
GSL::Test::test_rel(min, expected, rel, "gsl_stats_minmax: min")

expected = 4
max_index = rawa.max_index
GSL::Test::test(max_index != expected, "gsl_stats_max_index")

expected = 3
min_index = rawa.min_index
GSL::Test::test(min_index != expected, "gsl_stats_min_index")

min_index, max_index = rawa.minmax_index
expected = 4
GSL::Test::test(max_index != expected, "gsl_stats_minmax_index: max")

expected = 3
GSL::Test::test(min_index != expected, "gsl_stats_minmax_index: min")

rawa.sort!

median = rawa.median_from_sorted_data()
expected = 0.07505
GSL::Test::test_rel(median, expected, rel, "gsl_stats_median_from_sorted_data")

median = rawa.subvector(0, rawa.size-1).median_from_sorted_data()
expected = 0.0728
GSL::Test::test_rel(median, expected, rel, "gsl_stats_median_from_sorted_data")

quantile = rawa.quantile_from_sorted_data(0.0)
expected = 0.0242
GSL::Test::test_rel(quantile, expected, rel, "gsl_stats_quantile_from_sorted_data")

quantile = rawa.quantile_from_sorted_data(1.0)
expected = 0.1331
GSL::Test::test_rel(quantile, expected, rel, "gsl_stats_quantile_from_sorted_data (100)")

quantile = rawa.quantile_from_sorted_data(0.5)
expected = 0.07505
GSL::Test::test_rel(quantile, expected, rel, "gsl_stats_quantile_from_sorted_data (50even)")

quantile = rawa.subvector(0, rawa.size-1).quantile_from_sorted_data(0.5)
expected = 0.0728
GSL::Test::test_rel(quantile, expected, rel, "gsl_stats_quantile_from_sorted_data (50odd)")
