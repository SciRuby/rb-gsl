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
# == shortcuts.rb
#
# These are shortcuts for NMatrix creation, contributed by Daniel
# Carrera (dcarrera@hush.com) and Carlos Agarie (carlos@onox.com.br).

class NMatrix

  class << self
    # zeros() or zeroes()
    #
    #     Creates a new matrix of zeros with the dimensions supplied as
    #     parameters. Optional parameters include:
    #
    #     * A storage type as the first parameter (default is :dense).
    #     * A dtype as the last parameter (default is :float64).
    #
    # Examples:
    #
    #   zeros(3) # =>  0.0   0.0   0.0
    #
    #   zeros(2,3,:int32) # =>  0  0  0
    #                           0  0  0
    #
    #   zeros(:list, 5, :int32) # =>  0  0  0  0  0
    #
    
    def zeros(*params)
      dtype = params.last.is_a?(Symbol) ? params.pop : :float64
      stype = params.first.is_a?(Symbol) ? params.shift : :dense

      NMatrix.new(stype, params, 0, dtype)
    end

    alias :zeroes :zeros

    # ones()
    #
    #     Creates a :dense matrix of ones with the dimensions supplied
    #     as parameters. Optionaly, one can specify a dtype as the last
    #     parameter (default is :float64).
    #
    # Examples:
    #
    #   ones(3) # =>  1.0   1.0   1.0
    #
    #   ones(2,3,:int32) # =>  1  1  1
    #                          1  1  1
    #
    
    def ones(*params)
        dtype = params.last.is_a?(Symbol) ? params.pop : :float64

        NMatrix.new(params, 1, dtype)
    end

    # identity() or eye()
    #
    #     Creates an identitiy matrix (square matrix rank 2) of the size
    #     supplied as a parameter. Optional parameters include:
    #
    #     * A storage type as the first parameter (default is :dense).
    #     * A dtype as the last parameter (default is :float64).
    #
    # Examples:
    #
    #    eye(3) # =>   1.0   0.0   0.0
    #                  0.0   1.0   0.0
    #                  0.0   0.0   1.0
    #
    #    eye(3, :int32) # =>   1   0   0
    #                          0   1   0
    #                          0   0   1
    #
    #    eye(:yale, 2, :int32) # =>   1   0
    #                                 0   1
    #

    def eye(*params)
      dtype = params.last.is_a?(Symbol) ? params.pop : :float64
      stype = params.first.is_a?(Symbol) ? params.shift : :dense

      n = params[0]
      matrix = zeros(stype, n, n, dtype)
      (0..n-1).each { |i| matrix[i,i] = 1 }
      matrix
    end

    alias :identity :eye
    
    # rand()
    #
    #     Creates a :dense NMatrix with random numbers between 0 and 1 generated
    #     by Random::rand. The parameters to the method are the dimensions of the
    #     matrix.
    #
    # Examples:
    #
    #   rand(3,3) # => 0.4859439730644226   0.1783195585012436
    #                  0.23193766176700592  0.4503345191478729
    #

    def rand(*params)
      rng = Random.new

      product = params.reduce(1) { |prod, n| prod *= n }

      random_values = []
      product.times { |i| random_values << rng.rand }

      NMatrix.new(params, random_values, :float64)
    end

    # seq()
    #
    #     Creates a :dense NMatrix with a sequence of integers starting at
    #     zero until the matrix is filled. The parameters to the method
    #     are the dimensions of the matrix. Optionaly, one can specify a
    #     dtype as the last parameter (default is :float64).
    #
    # Examples:
    #
    #   seq(4) # =>   0   1   2   3
    #
    #   seq(3,3, :float32) # =>  0.0  1.0  2.0
    #                            3.0  4.0  5.0
    #                            6.0  7.0  8.0
    #
    
    def seq(*params)
      dtype = params.last.is_a?(Symbol) ? params.pop : nil
      
      # Already popped the dtype, if one is provided.
      is_vector = (params.size == 1)
      
      product = params.reduce(1) { |prod, n| prod *= n }
      
      # TODO: Is there a cleaner way to handle this?
      if dtype
        if is_vector
          NVector.new(params, (0..product-1).to_a, dtype)
        else
          NMatrix.new(params, (0..product-1).to_a, dtype)
        end
      else
        if is_vector
          NVector.new(params, (0..product-1).to_a)
        else
          NMatrix.new(params, (0..product-1).to_a)
        end
      end
    end

    #########################################
    # FUNCTIONS FOR MATLAB AND IDL REFUGEES #
    #########################################

    #
    # These are functions that replicate existing functionality, but
    # would probably be appreciated by MATLAB or IDL users.
    #

    # indgen() , findgen() , bindgen() , cindgen()
    #
    #      These IDL functions are similar to seq() but less flexible.
    #      They produce one-dimensional vectors:
    #
    #      indgen    -- Integer vector   --   seq(n, :int32)
    #      findgen   -- Float vector     --   seq(n, :float32)
    #      bindgen   -- Byte vector      --   seq(n, :byte)
    #      cindgen   -- Complex vector   --   seq(n, :complex64)
    #

    def indgen(n)
      seq(n, :int32)
    end

    def findgen(n)
      seq(n, :float32)
    end

    def bindgen(n)
      seq(n, :byte)
    end

    def cindgen(n)
      seq(n, :complex64)
    end

    # linspace()
    #
    #      This MATLAB function somewhat resembles seq(), but it differs
    #      enough that is likely to be legitimately useful to non-MATLAB
    #      refugees. This function takes three parameter, the last one
    #      being an integer.
    #
    #      linspace( a, b, n )
    #
    #      This returns a vector with n values from a to b.
    #
    # Ex:  x = linspace(0, pi, 1000)
    #      y = sin(x)
    #

    def linspace(a, b, nsteps)
      #
      # Algorithm:  seq(n) * (b-a)/(n-1) + a
      #
      step = (b-a) * 1.0 / (nsteps - 1)
      seq(nsteps) * NMatrix.new([nsteps], step) + NMatrix.new([nsteps], a)
    end
  end
end


# NMatrix needs to have a succinct way to create a matrix by specifying
# the components directly. This is very usefeul for using NMatrix as an
# advanced calculator, it is useful for learning NMatrix and it is also
# useful for testing language features or developing algorithms.
#
# The N[] function provides a way to create a matrix in a way that is
# very short and very natural, simply by specifying the components in
# the traditional Ruby array syntax.  Optionally, one can specify a
# dtype as the last parameter (default is :float64).
#
# a = N[ 1,2,3,4 ]          =>  1.0  2.0  3.0  4.0
#
# a = N[ 1,2,3,4, :int32 ]  =>  1  2  3  4
#
# a = N[ [1,2,3], [3,4,5] ] =>  1.0  2.0  3.0
#                               3.0  4.0  5.0
#
#
# SYNTAX COMPARISON:
#
#     MATLAB:		a = [ [1 2 3] ; [4 5 6] ]   or  [ 1 2 3 ; 4 5 6 ]
#     IDL:			a = [ [1,2,3] , [4,5,6] ]
#     NumPy:		a = array( [1,2,3], [4,5,6] )
#
#     SciRuby:      a = N[ [1,2,3], [4,5,6] ]
#     Ruby array:   a =  [ [1,2,3], [4,5,6] ]
#

class N
  class << self
    def [](*params)
      dtype = params.last.is_a?(Symbol) ? params.pop : nil

      #
      # First find the dimensions of the array.
      #
      i = 0
      dim = []
      foo = params
      while foo.is_a?(Array)
        dim[i] = foo.length
        foo    = foo[0]
        i += 1
      end

      #
      # Then flatten the array.
      #
      if dtype
        NMatrix.new(dim, params.flatten, dtype)
      else
        NMatrix.new(dim, params.flatten)
      end
    end
  end
end