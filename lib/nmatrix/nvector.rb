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

# NMatrix
require 'nmatrix'

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
	end

	# Orientation defaults to column (e.g., [3,1] is a column of length 3). It
	# may also be row, e.g., for [1,5].
	def orientation
		@orientation ||= :column
	end

	def transpose
		
		t = super()
		t.instance_exec { @orientation = @orientation == :row ? :column : :row }
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
		(0...shape[dim]).each { |i| arr << self[i] }
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
