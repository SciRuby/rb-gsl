require "packable"
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
# == io/matlab/mat_reader.rb
#
# Base class for .mat file reading (Matlab files).
#

require "packable"

module NMatrix::IO::Matlab
  # Class for parsing a .mat file stream.
  #
  # The full format of .mat files is available here:
  # * http://www.mathworks.com/help/pdf_doc/matlab/matfile_format.pdf
  #
  class MatFileReader
    MDTYPE_UNPACK_ARGS = {
      :miINT8 => [Integer, {:signed => true, :bytes => 1}],
      :miUINT8 => [Integer, {:signed => false, :bytes => 1}],
      :miINT16 => [Integer, {:signed => true, :bytes => 2}],
      :miUINT16 => [Integer, {:signed => false, :bytes => 2}],
      :miINT32 => [Integer, {:signed => true, :bytes => 4}],
      :miUINT32 => [Integer, {:signed => false, :bytes => 4}],
      :miSINGLE => [Float, {:precision => :single, :bytes => 4}],
      :miDOUBLE => [Float, {:precision => :double, :bytes => 8}],
      :miINT64 => [Integer, {:signed => true, :bytes => 8}],
      :miUINT64 => [Integer, {:signed => false, :bytes => 8}]
    }

    # Before release v7.1 (release 14) matlab (TM) used the system
    # default character encoding scheme padded out to 16-bits. Release 14
    # and later use Unicode. When saving character data, R14 checks if it
    # can be encoded in 7-bit ascii, and saves in that format if so.
    MDTYPES = [nil, :miINT8, :miUINT8, :miINT16, :miUINT16, :miINT32, :miUINT32, :miSINGLE, nil, :miDOUBLE, nil, nil, :miINT64, :miUINT64, :miMATRIX, :miCOMPRESSED, :miUTF8, :miUTF16, :miUTF32]
    MCLASSES = [nil, :mxCELL, :mxSTRUCT, :mxOBJECT, :mxCHAR, :mxSPARSE, :mxDOUBLE, :mxSINGLE, :mxINT8, :mxUINT8, :mxINT16, :mxUINT16, :mxINT32, :mxUINT32, :mxINT64, :mxUINT64, :mxFUNCTION, :mxOPAQUE, :mxOBJECT_CLASS_FROM_MATRIX_H]

    attr_reader :byte_order

    def initialize stream, options = {}
      raise(ArgumentError, "first arg must be IO") unless stream.is_a?(::IO)

      @stream     = stream
      @byte_order = options[:byte_order] || guess_byte_order
    end

    def guess_byte_order
      :native # assume native, since we don't know what type of file we have.
    end

  protected
    attr_reader :stream
  end
end