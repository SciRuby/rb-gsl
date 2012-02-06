require File.join(File.dirname(__FILE__), "nmatrix/nmatrix.so")

class NMatrix
  VERSION = '0.0.1'

  #def inspect
  #
  #end

  # TODO: Make this actually pretty.
  def pretty_print
    raise(NotImplementedError, "can only print rank 2 matrices") unless rank == 2
    (0...shape[0]).each do |i|
      arr = []
      (0...shape[1]).each do |j|
        arr << self[i,j]
      end
      puts arr.join("  ")
    end
    nil
  end


  def inspect
    original_inspect = super
    original_inspect = original_inspect[0...original_inspect.size-1]
    ary = [original_inspect]
    ary << "shape:[#{shape.join(',')}]" << "dtype:#{dtype}" << "stype:#{stype}"

    if stype == :yale
      ary << "capacity:#{capacity}" << "ija:#{__yale_ary__to_s(:ija)}" << "ia:#{__yale_ary__to_s(:ia)}" <<
             "ja:#{__yale_ary__to_s(:ja)}" << "a:#{__yale_ary__to_s(:a)}" << "d:#{__yale_ary__to_s(:d)}" <<
             "lu:#{__yale_ary__to_s(:lu)}" << "yale_size:#{__yale_size__}"
    end

    ary.join(" ") + ">"
  end

  def __yale_ary__to_s(sym)
    ary = self.send("__yale_#{sym.to_s}__".to_sym)
    "[" + ary.collect { |a| a.nil? ? "nil" : a }.join(',') + "]"
  end

  class << self

    def cblas_gemm a, b, c=nil, alpha=1.0, beta=0.0, transpose_a=false, transpose_b=false, m=nil, n=nil, k=nil, lda=nil, ldb=nil, ldc=nil
      raise(ArgumentError, "expected dense NMatrices as first two arguments") unless a.is_a?(NMatrix) && b.is_a?(NMatrix) && a.stype == :dense && b.stype == :dense
      raise(ArgumentError, "expected nil or dense NMatrix as third argument") unless c.nil? || (c.is_a?(NMatrix) && c.stype == :dense)
      raise(ArgumentError, "NMatrix dtype mismatch") unless a.dtype == b.dtype && (c.nil? ? true : a.dtype == c.dtype)

      # First, set m, n, and k, which depend on whether we're taking the transpose of a and b.
      if c.nil?
        if transpose_a # either :transpose or :complex_conjugate
          m ||= a.shape[1]
          k ||= a.shape[0]
        else # no transpose
          m ||= a.shape[0]
          k ||= a.shape[1]
        end
        n ||= transpose_b ? b.shape[0] : b.shape[1]
        c = NMatrix.new([m, n], a.dtype)
      else
        m ||= c.shape[0]
        n ||= c.shape[1]
        k ||= transpose_a ? a.shape[0] : a.shape[1]
      end

      # I think these are independent of whether or not a transpose occurs.
      lda ||= a.shape[1]
      ldb ||= b.shape[1]
      ldc ||= c.shape[1]

      # For argument descriptions, see: http://www.netlib.org/blas/dgemm.f
      NMatrix.__cblas_gemm__(transpose_a, transpose_b, m, n, k, alpha, a, lda, b, ldb, beta, c, ldc)

      return c
    end
  end


end
