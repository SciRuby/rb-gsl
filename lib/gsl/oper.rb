module GSL::Oper

  def self.included(base)
    base.class_eval {
      alias_method :_gsl_oper_original_mul, :*
      alias_method :_gsl_oper_original_div, :/

      def *(other)
        case other
          when Numeric
            _gsl_oper_original_mul(other)
          when GSL::Matrix,          GSL::Vector,
               GSL::Matrix::Int,     GSL::Vector::Int,
               GSL::Vector::Complex, GSL::Matrix::Complex,
               *GSL.have_tensor? ? [GSL::Tensor, GSL::Tensor::Int] : []
            other.scale(self)
          else
            _gsl_oper_original_mul(other)
        end
      end

      def /(other)
        case other
          when Numeric
            _gsl_oper_original_div(other)
          when GSL::Poly, GSL::Poly::Int
            a = GSL::Poly[1]; a[0] = self
            GSL::Rational.new(a, other)
          when GSL::Vector::Col
            other.scale(1.0 / GSL.pow_2(other.dnrm2))
          when GSL::Vector::Int::Col
            v = other.to_f
            v.scale(1.0 / GSL.pow_2(v.dnrm2))
          else
            _gsl_oper_original_div(other)
        end
      end
    }
  end

end

[Fixnum, Float].each { |klass| klass.send(:include, GSL::Oper) }
