class Fixnum
  alias :_orig_mul :*
  alias :_orig_div :/
  def *(other)
    if other.kind_of?(GSL::Matrix) or other.kind_of?(GSL::Vector) or other.kind_of?(GSL::Matrix::Int) or other.kind_of?(GSL::Vector::Int) or other.kind_of?(GSL::Vector::Complex)or other.kind_of?(GSL::Matrix::Complex)
      other.scale(self)
    else
      if GSL.have_tensor?
        if other.kind_of?(GSL::Tensor) or other.kind_of?(GSL::Tensor::Int)
          other.scale(self)
        else
          self._orig_mul(other)
        end
      else
        self._orig_mul(other)
      end
    end
  end

  def /(other)
    if other.kind_of?(GSL::Poly) or other.kind_of?(GSL::Poly::Int)
      a = GSL::Poly[1]; a[0] = self
      GSL::Rational.new(a, other)
    elsif other.kind_of?(GSL::Vector::Col) 
      other.scale(1.0/pow_2(other.dnrm2))
    elsif other.kind_of?(GSL::Vector::Int::Col)
      v = other.to_f
      v.scale(1.0/pow_2(v.dnrm2))
    else
      self._orig_div(other)
    end
  end
end

class Float
  alias :_orig_mul :*
  alias :_orig_div :/

  def *(other)
    if other.kind_of?(GSL::Matrix) or other.kind_of?(GSL::Vector) or other.kind_of?(GSL::Matrix::Int) or other.kind_of?(GSL::Vector::Int) or other.kind_of?(GSL::Vector::Complex)or other.kind_of?(GSL::Matrix::Complex)
      other.scale(self)
    else
      if GSL.have_tensor?
        if other.kind_of?(GSL::Tensor) or other.kind_of?(GSL::Tensor::Int)
          other.scale(self)
        else
          self._orig_mul(other)
        end
      else
        self._orig_mul(other)
      end
    end
  end

  def /(other)
    if other.kind_of?(GSL::Poly) or other.kind_of?(GSL::Poly::Int)
      a = GSL::Poly[1]; a[0] = self
      GSL::Rational.new(a, other)
    elsif other.kind_of?(GSL::Vector::Col) 
      other.scale(1.0/pow_2(other.dnrm2))
    elsif other.kind_of?(GSL::Vector::Int::Col)
      v = other.to_f
      v.scale(1.0/pow_2(v.dnrm2))
    else
      self._orig_div(other)
    end
  end
end
