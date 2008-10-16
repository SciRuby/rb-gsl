#!/usr/bin/env ruby
require("gsl")
include GSL

module GSL
  module Test
    Verbose = true
    $tests = 0
    $passed = 0
    $failed = 0

    def test(status, desc)
      $tests += 1
      if !status
        $passed += 1
        printf("PASS: #{desc}\n")
      else
        $failed += 1
        printf("FAIL: #{desc}\n")
      end
    end

    def test_factor(result, expected, factor, desc)
      status = nil
      if result == expected
        status = false
      elsif expected == 0.0
        status = (result > expected or result < expected)
      else
        u = result/expected
        status = (u > factor or u < 1.0/factor)
      end
      $tests += 1
      if !status
        $passed += 1
        printf("PASS: #{desc} (%g observed vs %g expected)\n", result, expected) 
      else
        $failed += 1
        printf("FAIL: #{desc} (%.18g observed vs %.18g expected)\n", result, expected) 
      end
    end

   def test_factor2(result, expected, factor, desc)
      status = nil
      if result == expected
        status = false
      elsif expected == 0.0
        status = (result > expected or result < expected)
      else
        u = result/expected
        status = (u > factor or u < 1.0/factor)
      end
      $tests += 1
      if !status
        $passed += 1
        printf("PASS: #{desc} (%.18g observed vs %g expected)\n", result, expected) 
      else
        $failed += 1
        printf("FAIL: #{desc} (%.18g observed vs %.18g expected)\n", result, expected) 
      end
    end

    def test_rel(result, expected, relerr, desc)
      status = nil
      if isnan?(result) or isnan?(expected)
        status = isnan?(result) != isnan?(expected)
      elsif isinf?(result) or isinf?(expected)
        status = isinf?(result) != isinf?(expected)
      elsif expected.to_f != 0.0
        status = (result - expected).abs/expected.abs > relerr
      else
        status = result.abs > relerr
      end
      $tests += 1
      if !status
        $passed += 1
        printf("PASS: #{desc} (%.18g observed vs %g expected)\n", result, expected) 
      else
        $failed += 1
        printf("FAIL: #{desc} (%.18g observed vs %.18g expected)\n", result, expected) 
      end
    end

    def test_abs(result, expected, abserr, desc)
      status = nil
      if isnan?(result) or isnan?(expected)
        status = isnan?(result) != isnan?(expected)
      elsif isinf?(result) or isinf?(expected)
        status = isinf?(result) != isinf?(expected)
      else
        status = (result - expected).abs > abserr
      end
      $tests += 1
      if !status
        $passed += 1
#        printf("PASS: #{desc} (%g observed vs %g expected)\n", result, expected) 
        printf("PASS: #{desc} (%.18g observed vs %g expected)\n", result, expected) 
      else
        $failed += 1
        printf("FAIL: #{desc} (%.18g observed vs %.18g expected)\n", result, expected) 
      end      
    end

    def test_int(result, expected, desc)
      status = (result != expected)
      $tests += 1
      if !status
        $passed += 1
        printf("PASS: #{desc} (%d observed vs %d expected)\n", result, expected) 
      else
        $failed += 1
        printf("FAIL: #{desc} (%d observed vs %d expected)\n", result, expected) 
      end
    end

  end
end

