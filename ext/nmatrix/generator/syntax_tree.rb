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
# == syntax_tree.rb
#
# Syntax Tree for mathematical expression parsing. Produces the
# correct forms for each type in NMatrix. Used by the C templates
# in generator.rb.
#

class Identifier < String
  def operate type=nil,exact_type=nil
    [type == :value || type == :object && to_s == "0" ? "RUBY_ZERO" : to_s]
  end

  def depth; 0; end

  # These are just dummy functions which alias to_s

  def operate_complex exact_type=nil
    [to_s]
  end

  def operate_rational exact_type=nil
    [to_s]
  end

  def operate_object exact_type=nil
    [to_s]
  end
end


class SyntaxTree
  # Changing the order of these arrays will change the relative precedence of the operators.
  ASSIGN_OPS = %w{+= -= *=}.map { |o| o.intern } # skipping: >>= <<= /= %= &= ^= |=
  COMP_OPS = %w{<= >= == != < >}.map { |o| o.intern }
  TRANSITIVE_BINARY_OPS = %w{&& || & ~ + *}.map { |o| o.intern }
  BINARY_OPS = %w{>> << && || & ~ - + * / %}.map { |o| o.intern }
  DISALLOWED_OPS = %w{++ -- !}.map { |o| o.intern }
  OPS = ASSIGN_OPS.concat(%w{<= >= == != < > = - + * / %}).map { |o| o.intern } # skipping: >> << && || & | ~
  OPS_REGEX = /[<>=!\+\-\*\/%]/
  EQ = :'='
  FLIP_OPS = {:'<=' => :'>=',
              :'>=' => :'<=',
              :'==' => :'==',
              :'!=' => :'!=',
              :'>' => :'<',
              :'<' => :'>'    }

  # Generally this is only used internally. Instead, call SyntaxTree.parse(operation), which will
  # in turn correctly call SyntaxTree::initialize.
  def initialize op, left, right
    if left == "0" && COMP_OPS.include?(op) # Flip certain operations so that the 0 is on the right.
      @op = FLIP_OPS[op]
      @left = right
      @right = left
    else
      @op = op
      @left = left
      @right = right
    end
  end
  attr_reader :op
  attr_accessor :left, :right


  def is_boolean?
    COMP_OPS.include?(op)
  end


  # Flip an operation, if possible. If not, do nothing.
  def reverse!
    if COMP_OPS.include?(op)
      @op = FLIP_OPS[op]
      tmp = @left
      @left = @right
      @right = tmp
    elsif TRANSITIVE_BINARY_OPS.include?(op)
      tmp = @left
      @left = @right
      @right = tmp
    end
  end

  # Flip an operation
  def reverse
    self.dup.reverse!
  end


  # Expand the syntax tree in place on ASSIGN_OPS. If @op is not one of the ASSIGN_OPS, do nothing.
  def expand!
    if ASSIGN_OPS.include?(@op)
      old_op    = @op
      @op       = EQ
      old_right = @right
      @right = SyntaxTree.new(old_op.to_s[0...old_op.size-1].intern, @left, old_right)
    end
    self
  end


  # Copy and expand the syntax tree on ASSIGN_OPS. e.g., += becomes = and +
  def expand
    self.dup.expand!
  end


  def depth
    (left.depth > right.depth ? left.depth : right.depth) + 1
  end

  # Split the SyntaxTree into a whole bunch of simpler ones appropriate for doing struct-type operations.
  def simplify_split! count=1
    expand!
    ary = [self]
    if op == EQ
      if depth <= 2
        ary
      else
        ary
        if right.left.is_a?(SyntaxTree)
          old_right_left = right.left
          right.left     = Identifier.new("temp#{count}")
          new_right_left = SyntaxTree.new(EQ, right.left, old_right_left)
          count += 1
          ary = ary.concat(new_right_left.simplify_split!(count))
        end

        if right.right.is_a?(SyntaxTree)
          old_right_right = right.right
          right.right     = Identifier.new("temp#{count}")
          new_right_right = SyntaxTree.new(EQ, right.right, old_right_right)
          count += 1
          ary = ary.concat(new_right_right.simplify_split!(count))
        end

        ary.reverse
      end
    else
      ary
    end
  end

  def simplify_split count=1
    dup.simplify_split! count
  end


  # Deep copy.
  def dup
    SyntaxTree.new(op, left.dup, right.dup)
  end


  # Split tree into multiple simpler subtrees and figure out the appropriate operation on each of those.
  #
  # Only needed for rational and complex, which are structs and can't rely on basic operators like * and +.
  def operate_on_subtrees type, exact_type
    i = 0
    best_type = begin
      case type
      when :rational
        :r128
      when :complex
        :c128
      end
    end
    subtrees = simplify_split
    subtrees.map do |subtree|
      i += 1
      i == subtrees.size ? subtree.send("operate_#{type}", exact_type) : subtree.send("operate_#{type}", best_type)
    end
  end


  def operate type, exact_type=nil
    operations =
    case type
    when :int, :float
      ["#{left.operate(type, exact_type).first} #{op} #{right.operate(type, exact_type).first}"]
    when :rational, :complex
      operate_on_subtrees(type, exact_type).flatten
    when :value,:object
      operate_object
    else
      raise ArgumentError, "undefined math expression type '#{type}'"
    end

    # Separate boolean operations by && or ; as appropriate for the type of operation.
    # This is overly simplistic, and won't work for things like x = y && z; but, it will
    # work for the majority of cases we need to deal with in our templates.
    #COMP_OPS.include?(op) ? operations.join(" && ") : operations.join(";\n") + ";"
    operations
  end


protected
  def operate_object
    ruby_right = !right.is_a?(SyntaxTree) && right.to_i.to_s == right ? "INT2FIX(#{right})" : right.operate_object.first
    case op
    when EQ
      ["#{left} #{op} #{ruby_right}"]
    when *COMP_OPS
      ["rb_funcall(#{left.operate_object.first}, rb_intern(\"#{op}\"), 1, #{ruby_right}) == Qtrue"]
    else
      ["rb_funcall(#{left.operate_object.first}, rb_intern(\"#{op}\"), 1, #{ruby_right})"]
    end
  end

  def operate_rational exact_type=:r128
    if right == "0"
      if op == EQ
        ["#{left}.n = 0", " #{left}.d = 1"]
      elsif COMP_OPS.include?(op)
        ["#{left}.n #{op} 0"]
      else
        raise NotImplementedError, "unhandled operation #{op} on rational and 0"
      end
    elsif right == "1"
      if op == EQ
        ["#{left}.n = #{left}.d = 1"]
      elsif COMP_OPS.include?(op)
        ["#{left}.n #{op} #{left}.d"]
      elsif [:'+', :'-'].include?(op)
        ["#{exact_type}_addsub(#{left}.n, #{left}.d, #{left}.d, #{left}.d, '#{op}')"]
      else
        raise NotImplementedError, "unhandled operation #{op} on rational and 1"
      end
    else
      if op == EQ
        ["#{left} = #{right.operate_rational(exact_type).first}"]
      elsif COMP_OPS.include?(op)
        raise NotImplementedError, "unhandled comparison #{op} on rationals"
      elsif [:'+', :'-'].include?(op)
        ["#{exact_type}_addsub(#{left}.n, #{left}.d, #{right}.n, #{right}.d, '#{op}')"]
      elsif [:'*', :'/'].include?(op)
        ["#{exact_type}_muldiv(#{left}.n, #{left}.d, #{right}.n, #{right}.d, '#{op}')"]
      else
        raise NotImplementedError, "unhandled operation #{op} on rationals"
      end
    end
  end

  def operate_complex exact_type=:c128
    if right == "0"
      if op == EQ
        ["#{left}.r = 0", " #{left}.i = 0"]
      elsif op == :'=='
        ["#{left}.r == 0 && #{left}.i == 0"]
      elsif op == :'!='
        ["(#{left}.r != 0 || #{left}.i != 0)"]
      elsif (COMP_OPS - [:'==']).include?(op)
        raise NotImplementedError, "unhandled comparison #{op} on complex and 0"
      else
        raise NotImplementedError, "unhandled operation #{op} on complex and 0"
      end
    else
      if op == EQ
        ["#{left} = #{right.operate_complex(exact_type).first}"]
      elsif op == :'=='
        ["#{left}.r == #{right}.r && #{left}.i == #{right}.i"]
      elsif (COMP_OPS - [:'==']).include?(op)
        raise NotImplementedError, "unhandled comparison #{op} on complex numbers"
      elsif op == :'+'
        ["#{exact_type}_add(#{left}.r, #{left}.i, #{right}.r, #{right}.i)"]
      elsif op == :'-'
        ["#{exact_type}_sub(#{left}.r, #{left}.i, #{right}.r, #{right}.i)"]
      elsif op == :'*'
        ["#{exact_type}_mul(#{left}.r, #{left}.i, #{right}.r, #{right}.i)"]
      elsif op == :'/'
        ["#{exact_type}_div(#{left}.r, #{left}.i, #{right}.r, #{right}.i)"]
      else
        raise NotImplementedError, "unhandled operation #{op} on complex numbers"
      end
    end
  end

public
  class << self

    # Parse an operation template of some kind into a SyntaxTree or, if no operation present,
    # into an Identifier (basically a string)
    def parse str
      SyntaxTree::OPS.each do |op_symbol|
        op = op_symbol.to_s
        pos = str.index(op)

        while !pos.nil? # Need to look for each operator multiple times in a line

          raise(ArgumentError, "disallowed operation '#{op}' in expression") if DISALLOWED_OPS.include?(op_symbol)

          # Is the operator contained within brackets?
          lb = str.rindex('[', pos)
          unless lb.nil?
            rb = str.rindex(']', pos) # Perhaps that bracket is terminated?
            rb = str.index(']', pos+op.size) if rb.nil? || rb < lb # still have to worry; try looking after the operator

            if BINARY_OPS.include?(op_symbol) && !rb.nil? && lb < pos && rb >= pos
              # Operator IS within brackets. Let's start searching again after the brackets.
              #STDERR.puts "operator #{op_symbol} was found, but within brackets (#{lb}, #{pos}, #{rb}) '#{str}'"
              pos = str.index(op, rb+1)
              next
            end
          end

          # If we get this far, the operator was not found within brackets.

          # Split around the operator
          left  = str[0...pos].strip
          right = str[pos+op.size...str.size].strip
          #STDERR.puts "operator #{op_symbol} was found (#{pos}); parsing '#{left}' and '#{right}'"
          return SyntaxTree.new(op_symbol, SyntaxTree.parse(left), SyntaxTree.parse(right))

        end
      end
      STDERR.puts "Making identifier out of #{str}"
      Identifier.new(str)
    end
  end
end
