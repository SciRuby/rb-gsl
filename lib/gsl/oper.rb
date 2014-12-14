module GSL::Oper
  def self.included base
    base.class_eval do

      alias :_orig_mul :*
      alias :_orig_div :/

      def *(other)
        return _orig_mul(other) if other.is_a?(Numeric)

        if other.is_a?(GSL::Matrix) || other.is_a?(GSL::Vector) || other.is_a?(GSL::Matrix::Int) ||
           other.is_a?(GSL::Vector::Int) || other.is_a?(GSL::Vector::Complex) || other.is_a?(GSL::Matrix::Complex)
          other.scale(self)
        elsif GSL.have_tensor? && (other.is_a?(GSL::Tensor) || other.is_a?(GSL::Tensor::Int))
          other.scale(self)
        else
          _orig_mul(other)
        end
      end

      def /(other)
        return _orig_div(other) if other.is_a?(Numeric)

        if other.is_a?(GSL::Poly) || other.is_a?(GSL::Poly::Int)
          a = GSL::Poly[1]; a[0] = self
          GSL::Rational.new(a, other)
        elsif other.is_a?(GSL::Vector::Col)
          other.scale(1.0/GSL::pow_2(other.dnrm2))
        elsif other.is_a?(GSL::Vector::Int::Col)
          v = other.to_f
          v.scale(1.0/GSL::pow_2(v.dnrm2))
        else
          _orig_div(other)
        end
      end

    end
  end
end

Fixnum.send(:include, GSL::Oper)
Float.send(:include, GSL::Oper)
