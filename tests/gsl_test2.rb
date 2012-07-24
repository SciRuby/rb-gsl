#!/usr/bin/env ruby
require("gsl")
include GSL

module GSL
  module Test
    Verbose = true
    $tests = 0
    $passed = 0
    $failed = 0

# PASS if status == 0, FAIL otherwise
    def test(status, desc)
      $tests += 1
      if status == 0
        $passed += 1
        printf("PASS: #{desc}\n")
      else
        $failed += 1
        printf("FAIL: #{desc}\n")
      end
    end

# PASS if status == true, FAIL otherwise
    def test2(status, desc)
      if status == true
        printf("PASS: #{desc}\n")
      else
        printf("FAIL: #{desc}\n")
      end
    end

    def test_factor(result, expected, factor, desc)
      status = nil
      if result == expected
        status = 0
      elsif expected == 0.0
        status = (result > expected or result < expected) ? 1 : 0
      else
        u = result/expected
        status = (u > factor or u < 1.0/factor) ? 1 : 0
      end
      $tests += 1
      if status == 0
        $passed += 1
        printf("PASS: #{desc} (%g observed vs %g expected)\n", result, expected) 
      else
        $failed += 1
        printf("FAIL: #{desc} (%.18g observed vs %.18g expected)\n", result, expected) 
      end
    end

    def test_rel(result, expected, relerr, desc)
      status = nil
      if isnan?(result) or isnan?(expected)
        status = isnan?(result) != isnan?(expected) ? 1 : 0
      elsif isinf?(result) or isinf?(expected)
        status = isinf?(result) != isinf?(expected) ? 1 : 0
      elsif expected.to_f != 0.0
        status = (result - expected).abs/expected.abs > relerr ? 1 : 0
      else
        status = result.abs > relerr ? 1 : 0
      end
      $tests += 1
      if status == 0
        $passed += 1
        printf("PASS: #{desc} (%g observed vs %g expected)\n", result, expected) 
      else
        $failed += 1
        printf("FAIL: #{desc} (%.18g observed vs %.18g expected)\n", result, expected) 
      end
    end

    def test_abs(result, expected, abserr, desc)
      status = nil
      if isnan?(result) or isnan?(expected)
        status = isnan?(result) != isnan?(expected) ? 1 : 0
      elsif isinf?(result) or isinf?(expected)
        status = isinf?(result) != isinf?(expected) ? 1 : 0
      else
        status = (result - expected).abs > abserr ? 1 : 0
      end
      $tests += 1
      if status == 0
        $passed += 1
        printf("PASS: #{desc} (%g observed vs %g expected)\n", result, expected) 
      else
        $failed += 1
        printf("FAIL: #{desc} (%.18g observed vs %.18g expected)\n", result, expected) 
      end      
    end

    def test_int(result, expected, desc)
      status = (result != expected) ? 1 : 0
      $tests += 1
      if status == 0
        $passed += 1
        printf("PASS: #{desc} (%d observed vs %d expected)\n", result, expected) 
      else
        $failed += 1
        printf("FAIL: #{desc} (%d observed vs %d expected)\n", result, expected) 
      end
    end

  end
end

