require 'test_helper'

class DiffTest < GSL::TestCase

  def setup
    @f = [GSL::Function.alloc { |x| Math.exp(x) }]
    @df = [GSL::Function.alloc { |x| Math.exp(x) }]

    @f << GSL::Function.alloc { |x|
      x >= 0.0 ? x * Math.sqrt(x) : 0.0
    }
    @df << GSL::Function.alloc { |x|
      x >= 0.0 ? 1.5 * Math.sqrt(x) : 0.0
    }

    @f << GSL::Function.alloc { |x|
      x != 0.0 ? Math.sin(1.0 / x) : 0.0
    }
    @df << GSL::Function.alloc { |x|
      x != 0.0 ? -Math.cos(1.0 / x) / (x * x) : 0.0
    }

    @f << GSL::Function.alloc { |x| Math.exp(-x * x) }
    @df << GSL::Function.alloc { |x| -2.0 * x * Math.exp(-x * x) }

    @f << GSL::Function.alloc { |x| x * x }
    @df << GSL::Function.alloc { |x| 2.0 * x }

    @f << GSL::Function.alloc { |x| 1.0 / x }
    @df << GSL::Function.alloc { |x| -1.0 / (x * x) }
  end

  {
    'exp(x)'    => 1.0,
    'x^(3/2)'   => 0.1,
    'sin(1/x)'  => 0.45,
    'exp(-x^2)' => 0.5,
    'x^2'       => 0.0,
    '1/x'       => 10.0
  }.each_with_index { |(f, x), i|
    %w[central forward backward].each { |diff|
      define_method("test_#{diff}_#{i}") {
        expected = @df[i].eval(x)

        result, abserr = @f[i].send("diff_#{diff}", x)
        assert_abs result, expected, abserr, desc = "#{f}, x=#{x}, #{diff} diff"

        refute((result - expected).abs > abserr, '%s, valid error estimate' % desc)
      }
    }
  }

end
