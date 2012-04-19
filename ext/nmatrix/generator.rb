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
# == generator.rb
#
# Module for generating source files.

$RELATIVE_PATH = nil

$IN_MAKEFILE = begin
  dir_pwd_split = Dir.pwd.split('/')
  if dir_pwd_split.size >= 4 && dir_pwd_split[-4] == "tmp" # when running make by hand
    $RELATIVE_PATH = "../../../../"
    true
  elsif dir_pwd_split[-2] == "ext" # when building gem
    $RELATIVE_PATH = File.join(File.dirname(__FILE__), "../..")
    true
  else # when building in development dir
    $RELATIVE_PATH = File.dirname(__FILE__)
    false
  end
end

require File.join($RELATIVE_PATH, "lib/string.rb") # from the Makefile
require File.join($RELATIVE_PATH, "ext/nmatrix/generator/syntax_tree.rb")
require File.join($RELATIVE_PATH, "ext/nmatrix/generator/templater.rb")

class DTypeInfo < Struct.new(:enum, :sizeof, :sym, :id, :type)
  def max_macro
    typename = self.sizeof.to_s
    if typename.include?('_')
      ary = typename.split('_')
      typename = ary[0...ary.size-1].join('')
    end
    typename.upcase + "_MAX"
  end

  def min_macro
    typename = self.sizeof.to_s
    if typename.include?('_')
      ary = typename.split('_')
      typename = ary[0...ary.size-1].join('')
    end
    typename.upcase + "_MIN"
  end

  # What type would this be if we used the maximum number of bytes available?
  def long_dtype
    Generator::DTYPES.select { |x| x.type == self.type }.last
  end
end


class Array
  def max
    found_max   = nil
    self.each_index do |i|
      found_max = self[i] if found_max.nil? || self[i] > found_max
    end
    found_max
  end

  def min
    found_min   = nil
    self.each_index do |i|
      found_min = self[i] if found_min.nil? || self[i] < found_min
    end
    found_min
  end
end


module Generator
  SRC_DIR = File.join("ext", "nmatrix")
  DTYPES = [
      # enum            sizeof        sym           id      type
      [:NM_NONE,        0,            :none,        0,      :none],
      [:NM_BYTE,        :u_int8_t,    :byte,        :b,     :int],
      [:NM_INT8,        :int8_t,      :int8,        :i8,    :int],
      [:NM_INT16,       :int16_t,     :int16,       :i16,   :int],
      [:NM_INT32,       :int32_t,     :int32,       :i32,   :int],
      [:NM_INT64,       :int64_t,     :int64,       :i64,   :int],
      [:NM_FLOAT32,     :float,       :float32,     :f32,   :float],
      [:NM_FLOAT64,     :double,      :float64,     :f64,   :float],
      [:NM_COMPLEX64,   :complex64,   :complex64,   :c64,   :complex],
      [:NM_COMPLEX128,  :complex128,  :complex128,  :c128,  :complex],
      [:NM_RATIONAL32,  :rational32,  :rational32,  :r32,   :rational],
      [:NM_RATIONAL64,  :rational64,  :rational64,  :r64,   :rational],
      [:NM_RATIONAL128, :rational128, :rational128, :r128,  :rational],
      [:NM_ROBJ,        :VALUE,       :object,      :v,     :value],
      [:NM_TYPES,       0,            :dtypes,      0,      :none]
  ].map { |d| DTypeInfo.new(*d) }

  INDEX_DTYPES = DTYPES.select { |dtype| dtype.type == :int && dtype.id != :b }
  INTEGER_DTYPES = DTYPES.select { |dtype| dtype.type == :int }
  RATIONAL_DTYPES = DTYPES.select { |dtype| dtype.type == :rational }
  NONBLAS_DTYPES = DTYPES.select { |dtype| [:int,:rational,:value].include?(dtype.type) }
  COMPLEX_DTYPES = DTYPES.select { |dtype| dtype.type == :complex }
  FLOAT_DTYPES = DTYPES.select { |dtype| dtype.type == :float }
  OBJECT_DTYPES = DTYPES.select { |dtype| dtype.type == :value }
  ACTUAL_DTYPES = DTYPES.select { |dtype| dtype.type != :none }

  YIELD_REGEX = /%%=\ [^%%]*%%/


  DTYPES_ASSIGN = {
      :complex => { # Assign a complex to:
          :complex  => lambda {|l,r| "((#{l}*)p1)->r = ((#{r}*)p2)->r; ((#{l}*)p1)->i = ((#{r}*)p2)->i;" },
          :float    => lambda {|l,r| "*(#{l}*)p1 = ((#{r}*)p2)->r;" },
          :int      => lambda {|l,r| "*(#{l}*)p1 = ((#{r}*)p2)->r;" },
          :rational => lambda {|l,r| "rb_raise(rb_eNotImpError, \"I don't know how to assign a complex to a rational\");"  },
          :value    => lambda {|l,r| "*(VALUE*)p1 = rb_complex_new(rb_float_new(((#{r}*)p2)->r), rb_float_new(((#{r}*)p2)->i));" },
       },
      :float => {
          :complex  => lambda {|l,r| "((#{l}*)p1)->i = 0; ((#{l}*)p1)->r = *(#{r}*)p2;" },
          :float    => lambda {|l,r| "*(#{l}*)p1 = *(#{r}*)p2;" },
          :int      => lambda {|l,r| "*(#{l}*)p1 = *(#{r}*)p2;" },
          :rational => lambda {|l,r| "rb_raise(rb_eNotImpError, \"I don't know how to assign a float to a rational\");" },
          :value    => lambda {|l,r| "*(VALUE*)p1 = rb_float_new(*(#{r}*)p2);" },
      },
      :int => {
          :complex  => lambda {|l,r| "((#{l}*)p1)->i = 0; ((#{l}*)p1)->r = *(#{r}*)p2;" },
          :float    => lambda {|l,r| "*(#{l}*)p1 = *(#{r}*)p2;" },
          :int      => lambda {|l,r| "*(#{l}*)p1 = *(#{r}*)p2;" },
          :rational => lambda {|l,r| "((#{l}*)p1)->d = 1; ((#{l}*)p1)->n = *(#{r}*)p2;" },
          :value    => lambda {|l,r| "*(VALUE*)p1 = INT2NUM(*(#{r}*)p2);" },
      },
      :rational => {
          :complex  => lambda {|l,r| "((#{l}*)p1)->i = 0; ((#{l}*)p1)->r = ((#{r}*)p2)->n / (double)((#{r}*)p2)->d;" },
          :float    => lambda {|l,r| "*(#{l}*)p1 = ((#{r}*)p2)->n / (double)((#{r}*)p2)->d;" },
          :int      => lambda {|l,r| "*(#{l}*)p1 = ((#{r}*)p2)->n / ((#{r}*)p2)->d;" },
          :rational => lambda {|l,r| "((#{l}*)p1)->d = ((#{r}*)p2)->d; ((#{l}*)p1)->n = ((#{r}*)p2)->n;" },
          :value    => lambda {|l,r| "*(VALUE*)p1 = rb_rational_new(INT2FIX(((#{r}*)p2)->n), INT2FIX(((#{r}*)p2)->d));" }
      },
      :value => {
          :complex  => lambda {|l,r| "((#{l}*)p1)->r = REAL2DBL(*(VALUE*)p2); ((#{l}*)p1)->i = IMAG2DBL(*(VALUE*)p2);" },
          :float    => lambda {|l,r| "*(#{l}*)p1 = NUM2DBL(*(VALUE*)p2);"},
          :int      => lambda {|l,r| "*(#{l}*)p1 = NUM2DBL(*(VALUE*)p2);"},
          :rational => lambda {|l,r| "((#{l}*)p1)->n = NUMER2INT(*(VALUE*)p2); ((#{l}*)p1)->d = DENOM2INT(*(VALUE*)p2);" },
          :value    => lambda {|l,r| "*(VALUE*)p1 = *(VALUE*)p2;"}
      }
  }


  class << self

    def decl spec_name, ary
      a = []
      a << "#{spec_name} {"
      ary.each do |v|
        a << "  #{v.to_s},"
      end
      a << "};"
      a.join("\n") + "\n\n"
    end


    def dtypes_err_functions
      str = <<SETFN
static void TypeErr(void) {
  rb_raise(rb_eTypeError, "illegal operation with this type");
}

SETFN
    end

    def dtypes_function_name func, dtype_i, dtype_j = nil
      if dtype_i[:enum] == :NM_NONE || (!dtype_j.nil? && dtype_j[:enum] == :NM_NONE)
        str = "TypeErr"
      else
        str = func.to_s.camelize
        str += "_#{dtype_i[:id]}"
        str += "_#{dtype_j[:id]}" unless dtype_j.nil?
      end
      str
    end

    def dtypes_assign lhs, rhs
      Generator::DTYPES_ASSIGN[ rhs.type ][ lhs.type ].call( lhs.sizeof, rhs.sizeof )
    end



    # Declare a set function for a pair of dtypes
    def dtypes_set_function dtype_i, dtype_j
      str = <<SETFN
static void #{dtypes_function_name(:set, dtype_i, dtype_j)}(size_t n, char* p1, size_t i1, char* p2, size_t i2) {
  for (; n > 0; --n) {
    #{dtypes_assign(dtype_i, dtype_j)}
    p1 += i1; p2 += i2;
  }
}

SETFN
    end

    def dtypes_increment_function dtype_i
      str = <<INCFN
static void #{dtypes_function_name(:increment, dtype_i)}(void* p) { ++(*(#{dtype_i[:sizeof]}*)p); }
INCFN
    end

    def dtypes_upcast
      ary = Array.new(15) { Array.new(15, nil) }
      DTYPES.each_index do |a|
        ad = DTYPES[a]
        (a...DTYPES.size).each do |b|
          bd = DTYPES[b]

          entry = nil

          if ad.type == :none || bd.type == :none
            entry ||= 'NM_NONE'
          elsif bd.type == ad.type
            entry ||= DTYPES[[a,b].max].enum.to_s
          elsif ad.type == :int # to float, complex, rational, or value
            entry ||= DTYPES[[a,b].max].enum.to_s
          elsif ad.enum == :NM_FLOAT32 # to complex or value
            if [:NM_FLOAT64, :NM_COMPLEX64, :NM_COMPLEX128, :NM_ROBJ].include?(bd.enum)
              entry ||= DTYPES[b].enum.to_s
            elsif [:NM_RATIONAL32, :NM_RATIONAL64, :NM_RATIONAL128].include?(bd.enum)
              entry ||= 'NM_FLOAT64'
            else
              entry ||= DTYPES[a].enum.to_s
            end
          elsif ad.enum == :NM_FLOAT64 # to complex or value
            if [:NM_COMPLEX128, :NM_ROBJ].include?(bd.enum)
              entry ||= DTYPES[b].enum.to_s
            elsif bd.enum == :NM_COMPLEX64
              entry ||= 'NM_COMPLEX128'
            else
              entry ||= DTYPES[a].enum.to_s
            end
          elsif ad.type == :rational # to float, complex, or value
            if [:NM_FLOAT64, :NM_COMPLEX128, :NM_ROBJ].include?(bd.enum)
              entry ||= DTYPES[b].enum.to_s
            elsif bd.enum == :NM_FLOAT32
              entry ||= 'NM_FLOAT64'
            elsif bd.enum == :NM_COMPLEX64
              entry ||= 'NM_COMPLEX128'
            else
              entry ||= DTYPES[a].enum.to_s
            end
          elsif ad.type == :complex
            if bd.enum == :NM_ROBJ
              entry ||= DTYPES[b].enum.to_s
            else
              entry ||= DTYPES[a].enum.to_s
            end
          elsif ad.type == :value # always value
            entry ||= DTYPES[a].enum.to_s
          end

          ary[a][b] = ary[b][a] = entry
        end
      end

      res = []
      ary.each_index do |i|
        res << "{ " + ary[i].join(", ") + " }"
      end
      decl("const int8_t Upcast[#{DTYPES.size}][#{DTYPES.size}] =", res) + "\n"
    end

    # binary-style functions, like Set (copy)
    def dtypes_binary_functions_matrix func
      ary = []
      DTYPES.each do |i|
        next if i[:enum] == :NM_TYPES
        bry = []
        DTYPES.each do |j|
          next if j[:enum] == :NM_TYPES
          bry << dtypes_function_name(func, i,j)
        end
        ary << "{ " + bry.join(", ") + " }"
      end
      ary
    end


    def dtypes_increment_functions_array
      ary = []
      DTYPES.each do |i|
        next if i[:enum] == :NM_TYPES
        if [:NM_INT8, :NM_INT16, :NM_INT32, :NM_INT64].include?(i.enum)
          ary << dtypes_function_name(:increment, i)
        else
          ary << dtypes_function_name(:increment, DTYPES[0]) # TypeErr
        end
      end
      ary
    end


    def dtypes_set_functions
      ary = []

      ary << dtypes_err_functions

      DTYPES.each do |dtype_i|
        DTYPES.each do |dtype_j|
          begin
            setfn = dtypes_set_function(dtype_i, dtype_j)
            ary << setfn unless setfn =~ /TypeErr/
          rescue NotImplementedError => e
            STDERR.puts "Warning: #{e.to_s}"
          rescue NoMethodError => e
            # do nothing
          end
        end
      end
      ary << ""
      ary << decl("nm_setfunc_t SetFuncs =", dtypes_binary_functions_matrix(:set))

      ary.join("\n")
    end

    def dtypes_increment_functions
      ary = []

      DTYPES.each do |dtype_i|
        next unless [:NM_INT8, :NM_INT16, :NM_INT32, :NM_INT64].include?(dtype_i.enum)
        incfn = dtypes_increment_function(dtype_i)
        ary << incfn unless incfn =~ /TypeErr/
      end

      ary << ""
      ary << decl("nm_incfunc_t Increment =", dtypes_increment_functions_array)

      ary.join("\n")
    end


    def dtypes_enum
      decl("enum NMatrix_DTypes", DTYPES.map{ |d| d[:enum].to_s })
    end

    def dtypes_sizeof
      decl("const int nm_sizeof[#{DTYPES.size}] =", DTYPES.map { |d| d[:sizeof].is_a?(Fixnum) ? d[:sizeof] : "sizeof(#{d[:sizeof].to_s})"})
    end

    def dtypes_typestring
      decl("const char *nm_dtypestring[] =", DTYPES.map { |d| "\"#{d[:sym].to_s}\"" })
    end


    def make_file filename, &block
      STDERR.puts "generated #{filename}"
      f = File.new(filename, "w")
      file_symbol = filename.split('.').join('_').upcase

      f.puts "/* Automatically created using generator.rb - do not modify! */"

      f.puts "#ifndef #{file_symbol}\n# define #{file_symbol}\n\n"
      yield f
      f.puts "\n#endif\n\n"
      f.close
    end


    def make_dtypes_c
      make_file "dtypes.c" do |f|
        f.puts dtypes_sizeof
        f.puts dtypes_typestring
        f.puts dtypes_upcast
      end
    end


    def make_dtypes_h
      make_file "dtypes.h" do |f|
        f.puts dtypes_enum
      end
    end


    def make_dfuncs_c
      make_file "dfuncs.c" do |f|
        f.puts '#include <ruby.h>'
        f.puts '#include "nmatrix.h"' + "\n\n"
        f.puts dtypes_set_functions
        f.puts dtypes_increment_functions
      end
    end


    # Read templates given by +names+ from <tt>SRC_DIR/relative_path</tt>, and output them to a filename described by
    # +output_name+.
    #
    # == Example
    #
    #    make_templated_c './smmp', 'header', %w{numbmm transp bstoy ytobs}, "smmp1.c", {:TYPE => RATIONAL_DTYPES, :INT => INDEX_DTYPES}
    #
    # TODO: index dtype is unsigned!
    # That means instead of int8_t, we should be doing uint8_t. But can't always do that because the Fortran code
    # occasionally starts at index -1. Stupid Fortran! Someone needs to go through and fix the code by hand.
    #
    # TODO: Write tests to confirm that the signedness isn't screwing up stuff.
    #
    def make_templated_c relative_path, header_name, names, output_name, subs = {:TYPE => INDEX_DTYPES}

      # First print the header once
      `cat #{$RELATIVE_PATH}/#{SRC_DIR}/#{relative_path}/#{header_name}.template.c > ./#{output_name}` unless header_name.nil?

      subs[:TYPE].each do |type|
        if subs.has_key?(:INT)
          subs[:INT].each do |int|
            names.each do |name|
              template "#{$RELATIVE_PATH}/#{SRC_DIR}/#{relative_path}/#{name}.template.c", output_name, :TYPE => type, :INT => int
            end
          end
        else
          names.each do |name|
            template "#{$RELATIVE_PATH}/#{SRC_DIR}/#{relative_path}/#{name}.template.c", output_name, :TYPE => type
          end
        end
      end
    end


    # Evaluate one-line Ruby statements embedded in a template.
    def gsub_yield line, t, dtype, line_number=nil, filename=nil
      match      = line.match YIELD_REGEX
      while !match.nil?

        statement = match[0][4...-2]
        result = self.send :eval, statement, binding, filename, line_number
        line["%%= #{statement}%%"] = result.to_s

        match      = line.match YIELD_REGEX
      end
      line
    end


    def gsub_expression_re re, line, t, dtype, line_number=nil, filename=nil
      match      = line.match re
      while !match.nil?
        expression = match[0][t.size+3...-2]
        operation  = SyntaxTree.parse(expression)

        begin
          operation_output = operation.operate(dtype.type, dtype.id)

          # Correctly join together the lines of output operations and insert them into the template line
          if operation.is_boolean?
            line["%%#{t} #{expression}%%"] = operation_output[0]
          else
            line["%%#{t} #{expression}%%"] = operation_output.join(";\n") + ";"
          end

        rescue NotImplementedError => e
          STDERR.puts "Error: #{e.inspect}"
          raise(SyntaxError, "possible NotImplementedError (#{dtype.type}) in template #{filename}: #{line_number}: \"#{expression}\"")
        rescue IndexError
          raise(StandardError, "string not matched: '%%#{t} #{expression}%%'")
        end

        match      = line.match re
      end
      line
    end


    # Replace a pseudo-mathematical expression with an actual one with dtypes taken into account.
    def gsub_expression line, t, dtype, line_number=nil, filename=nil
      gsub_expression_re /%%#{t}\ .*?%%/, line, t, dtype, line_number, filename
    end

    def gsub_expression_long line, t, dtype, line_number=nil, filename=nil
      gsub_expression_re /%%#{t}_LONG\ .*?%%/, line, "#{t}_LONG", dtype.long_dtype, line_number, filename
    end

    # Takes a list of declarations and cleans it for insertion in a header file.
    #
    # * Removes inline keyword
    # * Removes static functions
    # * Removes variable names
    def process_declarations declarations
      declarations.map do |d|
        process_declaration d
      end.compact
    end

    # Helper for process_declarations that works on a single function prototype.
    #
    # * Removes variable names
    # * Removes inline keyword
    # * Returns nil if declaration is static, otherwise returns corrected prototype
    def process_declaration declaration
      tokens = declaration.split(' ')

      return nil if tokens.include?('static')
      declaration = tokens.delete_if { |t| t == 'inline'}.join(' ')

      tokens = declaration.split('(')
      arg_list = tokens.last.split(')').first

      # Remove variable names
      args = arg_list.split(',')
      args = args.map do |arg|
        arg_tokens = arg.strip.split(' ')
        arg_tokens[0...arg_tokens.size-1].join(' ')
      end

      tokens[tokens.size-1] = args.join(',') + ')'

      tokens.join('(') + ";"
    end


    # Replaces sub_int_real and sub_int.
    #
    # Allows more flexible substitutions. Pass a hash of templates, e.g., {:INT => INDEX_DTYPES[0], :REAL => RATIONAL_DTYPES[1]}, and
    # it'll produce all possible combinations thereof.
    #
    # At some point we should probably just switch to erb. This just started growing and pretty soon I realized
    # erb would likely have been a better option. Oh well.
    def template template_filepath, output_filepath, types = {}
      raise(ArgumentError, "expected substitution templates") if types.size == 0

      # Keep track of all declarations in this template
      declarations = []

      # Process the current declaration
      block_level = 0
      declaration = ""
      decl_probably_finished = false

      in_comment = false

      #STDERR.puts "output_filepath = #{output_filepath}; Dir.pwd = #{Dir.pwd}"

      output   = File.new output_filepath, "a" # append
      template = File.new template_filepath, "r"

      line_count = 1

      while line = template.gets
        line.chomp!

        types.each_pair do |t_sym,dtype|
          t = t_sym.to_s

          if in_comment && line.include?("*/")
            m = line.split("*/", 1)
            m.shift
            line = m.first
            in_comment = false
          end

          # Are we in a multi-line C-style comment?
          unless in_comment
            # Ignore C-style single-line comments
            while m = line.match(/\/\*[^\*\/]*\*\//)
              line = m.pre_match + m.post_match
            end

            if line.include?("/*")
              line = line.split("/*").first
              in_comment = true
            end
          end

          # Remove C++-style comments
          line = line.split("//")[0] if line.include?("//")

          #STDERR.puts "Processing #{template_filepath}: #{line}"
          if line.include?("%%#{t}")
            line.gsub! "%%#{t}%%", dtype.sizeof.to_s
            line.gsub! "%%#{t}_ABBREV%%", dtype.id.to_s
            line.gsub! "%%#{t}_MAX%%", dtype.max_macro
            line.gsub! "%%#{t}_LONG%%", dtype.long_dtype.sizeof.to_s #e.g., int64 instead of int8 for temporary variables

            # Get any mathematical expressions that need to be translated
            line = gsub_expression(line, t, dtype, line_count, template_filepath)

            # Do the same for temp variables (which are often going to be more bytes)
            line = gsub_expression_long(line, t, dtype, line_count, template_filepath)
          end

          # Deal with any Ruby statements in the template.
          if line.include?("%%=")
            line = gsub_yield(line, t, dtype, line_count, template_filepath)
          end
        end

        unless in_comment
          # If we're not in a block, we should look for a function prototype.
          if block_level == 0
            maybe_prototype = line.split('{')[0] || ""
            if maybe_prototype !~ /;/
              declaration += maybe_prototype
            end

            paren_level = declaration.scan(/\(/).size
            if paren_level > 0 && paren_level == declaration.scan(/\)/).size
              decl_probably_finished = true
            else
              decl_probably_finished = false
            end
          end

          # Keep track of the block level to make sure we can identify function prototypes.
          block_level += line.scan(/{/).size

          # Found {, so prototype is probably finished. Add it to the declarations list
          if block_level > 0 && decl_probably_finished
            declarations << declaration
            declaration = ""
            decl_probably_finished = false
          end
          block_level -= line.scan(/}/).size
        end

        line_count += 1
        output.puts line
      end

      output.close

      declarations
    end

  end
end

if $IN_MAKEFILE
  Generator.make_dtypes_h
  Generator.make_dtypes_c
  Generator.make_dfuncs_c

  #
  # Order matters for these templates! Many functions are static.
  #

  Generator::Templater.new('smmp1.c', :in => 'yale', :boilerplate => 'smmp1_header') do |c|
    # 1-type interface functions for SMMP
    c.template 'smmp1', :TYPE => Generator::INDEX_DTYPES
    # 2-type interface functions for SMMP
    c.template 'smmp2', :TYPE => Generator::ACTUAL_DTYPES, :INT => Generator::INDEX_DTYPES

    c.update_header 'nmatrix'
  end

  Generator::Templater.new('smmp2.c', :in => 'yale', :boilerplate => 'smmp2_header') do |c|
    # 1-type SMMP functions from Fortran
    c.template 'symbmm', :TYPE => Generator::INDEX_DTYPES

    c.template 'complexmath', :in => 'shared', :TYPE => Generator::COMPLEX_DTYPES

    # Elementwise operations
    c.template 'elementwise_op', :TYPE => Generator::ACTUAL_DTYPES

    # 2-type SMMP functions from Fortran and selection sort
    c.template %w{numbmm transp sort_columns elementwise}, :TYPE => Generator::ACTUAL_DTYPES, :INT => Generator::INDEX_DTYPES

    c.update_header 'nmatrix'
  end

  Generator::Templater.new('blas.c', :in => 'dense', :boilerplate => 'blas_header') do |c|
    c.template 'rationalmath', :in => 'shared', :TYPE => Generator::RATIONAL_DTYPES
    c.template 'complexmath', :in => 'shared', :TYPE => Generator::COMPLEX_DTYPES

    # Functions derived from BLAS but adapted for rationals, integers, and Ruby objects
    c.template %w{gemm gemv}, :TYPE => Generator::NONBLAS_DTYPES

    # Elementwise operations
    c.template 'elementwise', :TYPE => Generator::ACTUAL_DTYPES

    c.update_header 'nmatrix'
  end

end