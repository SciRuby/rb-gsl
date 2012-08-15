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
# This file defines the NVector class.

############
# Requires #
############


#######################
# Classes and Modules #
#######################

# This is a specific type of NMatrix in which only one dimension is not 1.
# Although it is stored as a rank-2, n x 1, matrix, it acts as a rank-1 vector
# of size n. If the @orientation flag is set to :row, it is stored as 1 x n
# instead of n x 1.
class NVector < NMatrix
	def initialize(length, *args)
		super(:dense, [length, 1], *args)
    orientation
	end

	# Orientation defaults to column (e.g., [3,1] is a column of length 3). It
	# may also be row, e.g., for [1,5].
	def orientation
		@orientation ||= :column
	end

	def transpose
		super().flip
	end

	def transpose!
		super()
		self.flip
	end

	def multiply(m)
		super(m).flip
	end

	def multiply!(m)
		super().flip
	end

	def [](i)
		case @orientation
		when :row;		super(i, 0)
		when :column;	super(0, i)
		end
  end

	def []=(i, val)
		case @orientation
		when :row;		super(i, 0, val)
		when :column;	super(0, i, val)
		end
	end

	def rank; 1; end

	# TODO: Make this actually pretty.
	def pretty_print
		dim = @orientation == :row ? 1 : 0
		
		puts (0...shape[dim]).inject(Array.new) { |a, i| a << self[i] }.join('  ')
	end

protected
	def inspect_helper
		super() << "orientation:#{self.orientation}"
	end
	
	def flip_orientation
		returning(self) { @orientation = @orientation == :row ? :column : :row }
	end
	alias :flip :flip_orientation
end
