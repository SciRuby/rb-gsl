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
# == syntax_tree_spec.rb
#
# Tests for the source code generator (syntax_tree.rb in the
# ext/nmatrix/generator dir).
#
require "./ext/nmatrix/generator/syntax_tree.rb"

describe SyntaxTree do
  #it "correctly handles rational == 0" do
  #  SyntaxTree.parse("x == 0").operate(:rational,:r128).should == "x.n == 0"
  #end
  #
  #it "correctly handles rational = 0" do
  #  SyntaxTree.parse("x == 0").operate(:rational,:r128).split("\n").should == ["x.n = 0", "x.d = 1;"]
  #end
  #
  #it "correctly handles complex == 0" do
  #  SyntaxTree.parse("x == 0").operate(:complex,:c128).should == "x.r == 0 && x.i == 0"
  #end
  #
  #it "correctly handles complex = 0" do
  #  SyntaxTree.parse("x = 0").operate(:complex,:c128).split("\n").should == ["x.r = 0", "x.i = 1;"]
  #end
end