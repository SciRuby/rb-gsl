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
# == nmatrix.rb
#
# This file loads the C extension for NMatrix, and adds a few
# additional pieces of functionality (e.g., inspect, pretty_print).
# Also provided is NVector, which represents a rank-1 NMatrix in
# vector operations.

# For some reason nmatrix.so ends up in a different place during gem build
if File.exist? "lib/nmatrix/nmatrix.so"
  require File.join(File.dirname(__FILE__), "nmatrix/nmatrix.so") # development
else
  require File.join(File.dirname(__FILE__), "nmatrix.so")         # gem
end
require File.join(File.dirname(__FILE__), "array.rb") # Load Array extensions


class NMatrix
  # Read and write extensions for NMatrix. These are only loaded when needed.
  module IO
    autoload(:Matlab, File.join(File.dirname(__FILE__), 'nmatrix', 'io', 'matlab.rb'))
  end

  # TODO: Make this actually pretty.
  def pretty_print
    raise(NotImplementedError, "can only print rank 2 matrices") unless rank == 2
    (0...shape[0]).each do |i|
      arr = []
      (0...shape[1]).each do |j|
        arr << (self[i,j].nil? ? "nil" : self[i,j])
      end
      puts arr.join("  ")
    end
    nil
  end
  alias :pp :pretty_print


  # Get the complex conjugate of this matrix. See also complex_conjugate! for an in-place operation (provided the dtype
  # is already :complex64 or :complex128).
  #
  # Does not work on list matrices, but you can optionally pass in the type you want to cast to if you're dealing with
  # a list matrix.
  def complex_conjugate(new_stype = nil)
    self.cast(new_stype || self.stype, NMatrix::upcast(dtype, :complex64)).complex_conjugate!
  end

  # Calculate the conjugate transpose of a matrix. If your dtype is already complex, this should
  # only require one copy (for the transpose).
  def conjugate_transpose
    self.transpose.complex_conjugate!
  end

  def hermitian?
    return false if self.rank != 2 || self.shape[0] != self.shape[1]
    if [:complex64, :complex128].include?(self.dtype)
      # TODO: Write much faster Hermitian test in C
      self.eql?(conjugate_transpose)
    else
      symmetric?
    end
  end


  def inspect
    original_inspect = super
    original_inspect = original_inspect[0...original_inspect.size-1]
    original_inspect + inspect_helper.join(" ") + ">"
  end

  def __yale_ary__to_s(sym)
    ary = self.send("__yale_#{sym.to_s}__".to_sym)
    "[" + ary.collect { |a| a.nil? ? "nil" : a }.join(',') + "]"
  end

  class << self

    # Helper function for loading a file in the first sparse format given here:
    #   http://math.nist.gov/MatrixMarket/formats.html
    #
    # Override type specifier (e.g., 'real') using :read_with => :to_f (or any other string-to-numeric conversion
    # function), and with :dtype => :float32 or :dtype => :int8 to force storage in a lesser type.
    def load_matrix_matrix_coordinate_file filename, options = {}
      f = File.new(filename, "r")

      func = options.has_key?(:read_with) ? options[:read_with] : nil
      dtype = options.has_key?(:dtype) ? options[:dtype] : nil

      line = f.gets
      raise(IOError, "incorrect file type specifier") unless line =~ /^%%MatrixMarket\ matrix\ coordinate/
      spec = line.split
      case spec[3]
        when 'real'
          func ||= :to_f
          dtype ||= :float64
        when 'integer'
          func ||= :to_i
          dtype ||= :int64
        when 'complex'
          func ||= :to_complex
          dtype ||= :complex128
        when 'rational'
          func = :to_rational
          dtype ||= :rational128
        else
          raise ArgumentError, "Unrecognized dtype"
      end unless !func.nil? && !dtype.nil?

      line = f.gets
      while line =~ /^%/
        line = f.gets
      end

      rows, cols, entries = line.split.collect { |x| x.to_i }

      matrix = NMatrix.new(:yale, [rows, cols], entries, dtype)

      entries.times do
        i,j,v = line.split
        matrix[i.to_i-1,j.to_i-1] = v.send(func)
      end

      matrix
    end


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

      if a.dtype == :complex64 || a.dtype == :complex128 # NM_COMPLEX64 and NM_COMPLEX128 both require complex alpha and beta
        alpha = Complex.new(1.0, 0.0) if alpha == 1.0
        beta  = Complex.new(0.0, 0.0) if beta == 0.0
      end

      # For argument descriptions, see: http://www.netlib.org/blas/dgemm.f
      NMatrix.__cblas_gemm__(transpose_a, transpose_b, m, n, k, alpha, a, lda, b, ldb, beta, c, ldc)

      return c
    end


    def cblas_gemv a, x, y=nil, alpha=1.0, beta=0.0, transpose_a=false, m=nil, n=nil, lda=nil, incx=nil, incy=nil
      m ||= transpose_a ? a.shape[1] : a.shape[0]
      n ||= transpose_a ? a.shape[0] : a.shape[1]

      lda ||= a.shape[1]
      incx ||= 1
      incy ||= 1

      if a.dtype == :complex64 || a.dtype == :complex128 # NM_COMPLEX64 and NM_COMPLEX128 both require complex alpha and beta
        alpha = Complex.new(1.0, 0.0) if alpha == 1.0
        beta  = Complex.new(0.0, 0.0) if beta == 0.0
      end

      NMatrix.__cblas_gemv__(transpose_a, m, n, alpha, a, lda, x, incx, beta, y, incy)

      return y
    end
  end

protected
  def inspect_helper
    ary = []
    ary << "shape:[#{shape.join(',')}]" << "dtype:#{dtype}" << "stype:#{stype}"

    if stype == :yale
      ary << "capacity:#{capacity}" << "ija:#{__yale_ary__to_s(:ija)}" << "ia:#{__yale_ary__to_s(:ia)}" <<
             "ja:#{__yale_ary__to_s(:ja)}" << "a:#{__yale_ary__to_s(:a)}" << "d:#{__yale_ary__to_s(:d)}" <<
             "lu:#{__yale_ary__to_s(:lu)}" << "yale_size:#{__yale_size__}"
    end

    ary
  end

end

#######################################
#
# SIMPLE IN-LINE CONSTRUCTOR
#
# by Daniel Carrera <dcarrera@hush.com>
#
#######################################
# 
# NMatrix needs to have a succinct way to create a matrix by specifying
# the components directly. This is very usefeul for using NMatrix as an
# advanced calculator, it is useful for learning NMatrix and it is also
# useful for testing language features or developing algorithms.
# 
# The N[] function provides a way to create a matrix in a way that is
# very short and very natural, simply by specifying the components in
# the traditional Ruby array syntax.  Optionaly, one can specify a
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

	def N.[](*params)
		dtype = params.last.is_a?(Symbol) ? params.pop : :float64
		
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
		NMatrix.new( dim, params.flatten, dtype )
	end
end

# This is a specific type of NMatrix in which only one dimension is not 1. Although it is stored as a rank-2, n x 1,
# matrix, it acts as a rank-1 vector of size n. If the @orientation flag is set to :row, it is stored as 1 x n instead
# of n x 1.
class NVector < NMatrix
  def initialize length, *args
    super :dense, [length,1], *args
  end

  # Orientation defaults to column (e.g., [3,1] is a column of length 3). It may also be row, e.g., for [1,5].
  def orientation
    defined?(@orientation) && !@orientation.nil? ? @orientation : :column
  end

  def transpose
    t = super
    t.send :eval, "@orientation = @orientation == :row ? :column : :row"
    t
  end

  def transpose!
    super
    @orientation = @orientation == :row ? :column : :row
    self
  end

  def multiply m
    v = super(m)
    v.send :eval, "@orientation = @orientation == :row ? :column : :row"
    v
  end

  def multiply! m
    super
    @orientation = @orientation == :row ? :column : :row
    self
  end

  def [] i
    @orientation == :row ? super(0,i) : super(i,0)
  end

  def []= i,val
    @orientation == :row ? super(0,i,val) : super(i,0,val)
  end

  def rank; 1; end

  # TODO: Make this actually pretty.
  def pretty_print
    dim = @orientation == :row ? 1 : 0
    arr = []
    (0...shape[dim]).each do |i|
      arr << self[i]
    end
    puts arr.join("  ")
    nil
  end

protected
  def inspect_helper
    ary = super
    ary << "orientation:#{(@orientation || 'column').to_s}"
    ary
  end

end

#######################################
#
# Shortcuts to initialize common matrices.
#
# by Daniel Carrera <dcarrera@hush.com>
#
##########################################

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
	dtype = params.last.is_a?(Symbol) ? params.pop   : :float64
	store = params.first.is_a?(Symbol) ? params.shift : :dense
	
	NMatrix.new(store, params, 0, dtype)
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
    dtype = params.last.is_a?(Symbol) ? params.pop   : :float64
    
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
	dtype = params.last.is_a?(Symbol) ? params.pop   : :float64
	store = params.first.is_a?(Symbol) ? params.shift : :dense
	
	n = params[0]
	m = zeros(store,n,n,dtype)
	(0..n-1).each { |i| m[i,i] = 1 }
	m
end

alias :identity :eye

# seq()
#      
#     Creates a :dense NMatrix with a sequence of integers starting at
#     zero until the matrix is filled. The parameters to the method
#     are the dimensions of the matrix. Optionaly, one can specify a
#     dtype as the last parameter (default is :float64).
#     
# Examples:
# 
#   seq(4) # =>   0.0   1.0   2.0   3.0
#
#   seq(3,3, :int32) # =>  0  1  2
#                          3  4  5
#                          6  7  8
#

def seq(*params)
	dtype = params.last.is_a?(Symbol) ? params.pop   : :float64
	
  product = params.reduce(1) { |prod, n| prod *= n }
  	
	NMatrix.new(params,  (0..prod-1).to_a, dtype )
end

# rand()
#      
#     Creates a :dense NMatrix with random numbers between 0 and 1 generated
#     by Random::rand. The parameters to the method are the dimensions of the
#     matrix.
#     
#     by Carlos Agarie <carlos@onox.com.br>
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
  
  NMatrix.new params, random_values, :float32
end

######################################
#
# FUNCIONS FOR MATLAB AND IDL REFUGEES
#
######################################
#
# These are functions that replicate existing functionality, but
# would probably be appreciated by MATLAB or IDL users.
#

# indgen() , findgen() , bindgen() , cindgen()
#      
#      These IDL functions are similar to seq() but less flexile.
#      They produce on-dimensional matrices:
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