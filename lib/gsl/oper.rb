module GSL::Oper
  def self.included base
    base.class_eval do

      alias :_orig_mul :*
      alias :_orig_div :/

      def *(other)
        return _orig_mul(other) if other.is_a?(Numeric)

        gsl_classes = [GSL::Matrix, GSL::Vector, GSL::Matrix::Int,
          GSL::Vector::Int, GSL::Vector::Complex, GSL::Matrix::Complex]
        if GSL.have_tensor?
          gsl_classes << GSL::Tensor << GSL::Tensor::Int
        end

        case other
        when *gsl_classes
          other.scale(self)
        else
          _orig_mul(other)
        end
      end

      def /(other)
        return _orig_div(other) if other.is_a?(Numeric)

        case other
        when GSL::Poly, GSL::Poly::Int
          a = GSL::Poly[1]; a[0] = self
          GSL::Rational.new(a, other)
        when GSL::Vector::Col
          other.scale(1.0/GSL::pow_2(other.dnrm2))
        when GSL::Vector::Int::Col
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
