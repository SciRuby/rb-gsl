#!/usr/bin/env ruby
require("gsl")
require("../gsl_test2.rb")
include GSL::Test
include Math

module GSL
  module Test
    module Sf
      TEST_TOL0 = 2.0*GSL::DBL_EPSILON
      TEST_TOL1 = 16.0*GSL::DBL_EPSILON
      TEST_TOL2 = 256.0*GSL::DBL_EPSILON
      TEST_TOL3 = 2048.0*GSL::DBL_EPSILON
      TEST_TOL4 = 16384.0*GSL::DBL_EPSILON
      TEST_TOL5 = 131072.0*GSL::DBL_EPSILON
      TEST_TOL6 = 1048576.0*GSL::DBL_EPSILON
      TEST_SQRT_TOL0 = 2.0*GSL::SQRT_DBL_EPSILON
      TEST_SNGL = 1.0e-06
      
      TEST_SF_INCONS = 1
      TEST_SF_ERRNEG = 2
      TEST_SF_TOLBAD = 4
      TEST_SF_RETBAD = 8
      
      TEST_FACTOR = 100.0
      TEST_SIGMA = 1.5
      
      def test_sf_frac_diff(x1, x2)
        if x1.zero? and x2.zero?
          return 0.0
        elsif x1 <= DBL_MAX and x2 < DBL_MAX and (x1 + x1 != 0.0)
          return ((x1 - x2)/(x1 + x2)).abs
        else
          return 1.0
        end
      end
      
      def test_sf_check_result(mbuf, r, val, tol)
        s = 0
        f = 0.0
        if isnan?(r.val) or isnan?(val)
          s = isnan(r.val) != isnan(val) ? TEST_SF_INCONS : s
        elsif isinf?(r.val) or isinf?(val)
          s = isinf(r.val) != isinf(val) ? TEST_SF_INCONS : s
        else
          f = test_sf_frac_diff(val, r.val)
          if (val - r.val).abs > 2.0*TEST_SIGMA*r.err
            s |= TEST_SF_INCONS
          end
          if r.err < 0.0
            s |= TEST_SF_ERRNEG
          end
          if f > TEST_FACTOR * tol
            s |= TEST_SF_TOLBAD
          end
        end
        if s != 0
          message = sprintf("  expected: %20.16g\n", val)
          mbuf += message
          message = sprintf("  obtained: %20.16g   %20.16g  %g\n", 
                            r.val, r.err, r.err/(r.val.abs + r.err))
          mbuf += message
          message = sprintf("  fracdiff: %20.16g\n", f)
          mbuf += message
        end
        
        if s & TEST_SF_INCONS
          mbuf += "  value/expected not consistent within reported error\n"
        end
        if s &  TEST_SF_ERRNEG
          mbuf += "  reported error negative\n"
        end
        if s & TEST_SF_TOLBAD
          mbuf += "  value not within tolerance of expected value\n"
        end
        return s, mbuf
      end
      
      def test_sf_check_val(mbuf, rval, val, tol)
        s = 0;
        f = test_sf_frac_diff(val, rval)
        if f > TEST_FACTOR * tol
          s |= TEST_SF_TOLBAD
        end
        if s != 0
          buf = sprintf("  expected: %20.16g\n", val)
          mbuf += buf
          buf = sprintf("  obtained: %20.16g\n", rval)
          mbuf += buf
          buf = sprintf("  fracdiff: %20.16g\n", f)
          mbuf += buf
        end
        if s & TEST_SF_TOLBAD
          mbuf +=  "  value not within tolerance of expected value\n"
        end
        return s, mbuf
      end
      
      def test_sf_check_result_relax(mbuf, r, val, tol)
        s = 0;
        f = test_sf_frac_diff(val, r.val)
        
        if f > GSL_MAX_DBL(TEST_SNGL, TEST_FACTOR * tol)
          s |= TEST_SF_INCONS
        end
        if r.err < 0.0
          s |= TEST_SF_ERRNEG
        end
        if f > TEST_FACTOR * tol    
          s |= TEST_SF_TOLBAD
        end
        if s != 0
          buf = sprintf("  expected: %20.16g\n", val)
          mbuf += buf
          buf = sprintf("  obtained: %20.16g   %20.16g  %g\n", r.val, r.err, r.err/(r.val.abs + r.err))
          mbuf += buf
          buf = sprintf("  fracdiff: %20.16g\n", f)
          mbuf += buf
        end
        if s & TEST_SF_INCONS
          mbuf += "  value/expected not consistent MAX(tol,SINGLE_PREC)\n"
        end
        if s & TEST_SF_ERRNEG
          mbuf += "  reported error negative\n"
        end
        if s & TEST_SF_TOLBAD
          mbuf += "  value not within tolerance of expected value\n"
        end
        return s, mbuf
      end
      
      def test_sf_check_return(mbuf, val, expected)
        if val != expected
          buf = sprintf("  unexpected return code: %d\n", val)
          mbuf += buf
          return TEST_SF_RETBAD, mbuf
        else
          return 0, mbuf
        end
      end
      
      def test_sf(r, val, tol, status, expect, desc)
        local_s = 0
        mbuf = ""
        s, mbuf = test_sf_check_result(mbuf, r, val, tol)
        local_s |= s
        s, mbuf = test_sf_check_return(mbuf, status, expect)
        local_s |= s
        GSL::Test::test(local_s, desc)
        if local_s != 0
          print(mbuf)
          printf("  %22.18g  %22.18g\n", r.val, r.err)
        end
        return local_s
      end
      
      def test_sf_val(val, val_in, tol, desc)
        local_s = 0
        mbuf = ""
        s, mbuf = test_sf_check_val(mbuf, val, val_in, tol)
        local_s |= s
        GSL::Test::test(local_s, desc)
        if local_s != 0
          printf("%s", mbuf)
          printf("  %22.18g\n", val)
        end
        return local_s
      end
      
      def test_sf_rlx(r, val_in, tol, status, expect_return, desc)
        local_s = 0
        message_buff = ""
        s, message_buff = test_sf_check_result_relax(message_buff, r, val_in, tol)
        local_s |= s
        s, message_buff = test_sf_check_return(message_buff, status, expect_return)
        local_s |= s
        
        GSL::Test::test(local_s, desc)
        if local_s != 0
          printf("%s", message_buff)
          printf("  %22.18g  %22.18g\n", r.val, r.err)
        end
        return local_s
      end
      
      def test_sf_2(r1, val1, tol1, r2, val2, tol2, status, expect_return, desc)
        int local_s = 0
        message_buff = ""
        s, message_buff = test_sf_check_result(message_buff, r1, val1, tol1)
        local_s |= s
        s, message_buff = test_sf_check_result(message_buff, r2, val2, tol2)
        local_s |= s
        s, message_buff = test_sf_check_return(message_buff, status, expect_return)
        local_s |= s
        
        GSL::Test::test(local_s, desc)
        if local_s != 0
          printf("%s", message_buff)
          printf("  %22.18g  %22.18g\n", r1.val, r1.err)
          printf("  %22.18g  %22.18g\n", r2.val, r2.err)
        end
        return local_s
      end
      
      def test_sf_sgn(r, sgn, val_in, tol, expect_sgn, status, expect_return, desc)
        local_r
        local_s = 0
        message_buff = ""
        local_r.val = sgn
        local_r.err = 0.0
        s, message_buff = test_sf_check_result(message_buff, r, val_in, tol)
        local_s |= s
        s, message_buff = test_sf_check_result(message_buff, local_r, expect_sgn, 0.0)
        local_s |= s
        s, message_buff = test_sf_check_return(message_buff, status, expect_return)
        local_s |= s
        GSL::Test::test(local_s, desc)
        if local_s != 0
          printf("%s", message_buff)
          printf("  %22.18g  %22.18g\n", r.val, r.err)
        end
        return local_s
      end
      
      def TEST_SF(stat, func, args, val_in,  tol, expect_return)
        r, = eval("#{func}#{args}")
        #      p r
        #      p val_in
        status = 0
        stat += test_sf(r, val_in, tol, status, expect_return, "#{func}#{args}")
        return stat
      end
      def TEST_SF_2(stat, func, args, val1, tol1, val2, tol2, expect_return) 
        r, = eval("#{func}#{args}") 
        status = 0
        stat += test_sf_2(r, val1, tol1, r2, val2, tol2, status, 
                          expect_return, "#{func}#{args}")
        return stat
      end
      def TEST_SF_SGN(stat, func, args, val_in, tol, expect_sgn, expect_return)
        r, = eval("#{func}#{args}")
        status = 0
        stat += test_sf_sgn(r, sgn, val_in, tol, expect_sgn, status, 
                            expect_return, "#{func}#{args}")
        return stat
      end
    end
  end
end
