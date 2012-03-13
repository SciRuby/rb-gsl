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
# == blas.rb
#
# Generator module for BLAS functions.
#

module Generator
  module Blas
    PREFIXES = {:float32 => 's', :float64 => 'd', :complex64 => 'c', :complex128 => 'z'}

    FUNCTIONS = {:gemm => [:Order, :TransA, :TransB, :M, :N, :K, :alpha, :A, :lda, :B, :ldb, :beta, :C, :ldc],
                 :gemv => [:TransA, :TransB, :M, :N, :alpha, :A, :lda, :X, :incX, :beta, :Y, :incY],
    }

    TYPES = {:M => :size_t, :N => :size_t, :K => :size_t, :lda => :size_t, :incX => :int, :incY => :int}
    POINTERS = [:A, :B, :C, :X, :Y]
  end
end