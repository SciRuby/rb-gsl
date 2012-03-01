class Array
  # Convert array to an NMatrix
  #
  # == Arguments:
  # <tt>stype</tt> :: Optional storage type (defaults to :dense)
  # <tt>shape</tt> :: Array describing matrix dimensions (or Fixnum for square) -- REQUIRED!
  # <tt>dtype</tt> :: Override data type (e.g., to store a Float as :float32 instead of :float64) -- optional.
  def to_nm *args
    pos   = 0
    stype = args[pos].is_a?(Symbol) ? args[pos].tap { pos += 1} : :dense
    shape = args[pos]; pos += 1
    dtype = begin
      if pos >= args.size
        # TODO: Upcasting.
        if self[0].is_a?(Fixnum)
          :int64
        elsif self[0].is_a?(Float)
          :float64
        elsif self[0].is_a?(Rational)
          :rational128
        elsif self[0].is_a?(Complex)
          :complex128
        end
      else
        args[pos]
      end
    end

    m = NMatrix.new(stype, shape, dtype, self)
  end
end