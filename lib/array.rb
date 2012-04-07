# = NMatrix
#
# A linear algebra library for scientific computation in Ruby.
# NMatrix is part of SciRuby.
#
# NMatrix was originally inspired by and derived from NArray, by
# Masahiro Tanaka: http://narray.rubyforge.org
#
# == Copyright Information
#
# SciRuby is Copyright (c) 2010 - 2012, Ruby Science Foundation
# NMatrix is Copyright (c) 2012, Ruby Science Foundation
#
# Please see LICENSE.txt for additional copyright notices.
#
# == Contributing
#
# By contributing source code to SciRuby, you agree to be bound by
# our Contributor Agreement:
#
# * https://github.com/SciRuby/sciruby/wiki/Contributor-Agreement
#
# == array.rb
#
# Ruby core extensions for NMatrix.

class Array
  # Convert a Ruby Array to an NMatrix.
  #
  # You must provide a shape for the matrix as the first argument.
  #
  # == Arguments:
  # <tt>shape</tt> :: Array describing matrix dimensions (or Fixnum for square) -- REQUIRED!
  # <tt>dtype</tt> :: Override data type (e.g., to store a Float as :float32 instead of :float64) -- optional.
  # <tt>stype</tt> :: Optional storage type (defaults to :dense)
  def to_nm *args
    pos   = 0

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
        end.tap { pos += 1 }
      else
        args[pos].tap { pos += 1 }
      end
    end

    stype = args[pos].is_a?(Symbol) ? args[pos].tap { pos += 1} : :dense


    if stype == :dense
      NMatrix.new(stype, shape, self, dtype)
    else
      NMatrix.new(:dense, shape, self, dtype).cast(:stype, dtype)
    end
  end
end