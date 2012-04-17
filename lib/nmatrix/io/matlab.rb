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
# == io/matlab.rb
#
# Code for reading and writing Matlab matrix files.
#

require_relative "./matlab/mat5_reader"

class NMatrix
  module IO
    # IO components for Matlab.
    module Matlab
      class << self
        # Attempt to convert a Matlab .mat file's contents to a Ruby object.
        #
        # EXPERIMENTAL. At this time, only supports version 5.
        #
        def load_mat file_path
          NMatrix::IO::Matlab::Mat5Reader.new(File.open(file_path, "rb+")).to_ruby
        end
      end
    end
  end
end