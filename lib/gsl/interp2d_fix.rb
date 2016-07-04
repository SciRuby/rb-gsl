# Reasons for existence of this code:
#
# It so happens that GSL's 2D interpolation methods are somehow swapping X and
# Y co-ordinates and returning a wrong result. This is a bug in GSL, and in order
# to make up for that, we swap the x and y points before passing them to GSL.
# The detect_gsl_interp2d_swapping_bug method tests the inter2d eval function
# to see if the bug exists so that in case a future GSL update fixes the bug,
# this Ruby wrapper will remain unaffected.
#
# For testing whether the bug exists, the eval method is run and tested against
# what is known to be corrent output (ans_expected), and the @@swapped class
# variable is set to true if the result is not consistent with the expected
# output. If @swapped is true, the x and y values are internally swapped by
# the Ruby wrapper before passing them to GSL for processing.
module GSL
  class Interp2d
    module BugDetectHelper
      class << self
        def asymmetric_function(x, y)
          x - 2*y
        end
      end
    end

    def self.detect_gsl_interp2d_swapping_bug
      @@swapped = nil
      x = GSL::Vector.alloc((-4..4).to_a)
      y = GSL::Vector.alloc((-4..4).to_a)

      z = []
      x.each do |xi|
        y.each do |yi|
          z << BugDetectHelper.asymmetric_function(xi, yi)
        end
      end
      z = GSL::Vector.alloc(z)
      i2d = GSL::Interp2d.alloc(GSL::Interp2d::BICUBIC, x, y, z)

      test_x = x[1]
      test_y = y[1]

      ans_normal = i2d.eval(x,y,z,test_x, test_y)
      ans_swapped = i2d.eval(x,y,z,test_y, test_x)
      ans_expected = BugDetectHelper.asymmetric_function(test_x, test_y)
      @@swapped = ans_expected == ans_normal

      @@swapped
    end

    def self.swapped
      @@swapped
    end
  end

  class Spline2d
    @@swapped = GSL::Interp2d.detect_gsl_interp2d_swapping_bug

    def self.swapped
      @@swapped
    end
  end
end

GSL::Interp2d.detect_gsl_interp2d_swapping_bug