class DTypeInfo < Struct.new(:enum, :sizeof, :sym, :id, :type,  :gemm); end

module Generator
  SRC_DIR = File.join("ext", "nmatrix")
  DTYPES = [
      # dtype enum      sizeof        label/symbol/string   num-class
      [:NM_NONE,        0,            :none,        0,      :none,        :igemm],
      [:NM_BYTE,        :u_int8_t,    :byte,        :b,     :int,         :igemm],
      [:NM_INT8,        :int8_t,      :int8,        :i8,    :int,         :igemm],
      [:NM_INT16,       :int16_t,     :int16,       :i16,   :int,         :igemm],
      [:NM_INT32,       :int32_t,     :int32,       :i32,   :int,         :igemm],
      [:NM_INT64,       :int64_t,     :int64,       :i64,   :int,         :igemm],
      [:NM_FLOAT32,     :float,       :float32,     :f32,   :float,       :cblas_sgemm],
      [:NM_FLOAT64,     :double,      :float64,     :f64,   :float,       :cblas_dgemm],
      [:NM_COMPLEX64,   :complex64,   :complex64,   :c64,   :complex,     :cblas_cgemm],
      [:NM_COMPLEX128,  :complex128,  :complex128,  :c128,  :complex,     :cblas_zgemm],
      [:NM_RATIONAL32,  :rational32,  :rational32,  :r32,   :rational,    :rgemm],
      [:NM_RATIONAL64,  :rational64,  :rational64,  :r64,   :rational,    :rgemm],
      [:NM_RATIONAL128, :rational128, :rational128, :r128,  :rational,    :rgemm],
      [:NM_ROBJ,        :VALUE,       :object,      :v,     :value,       :vgemm],
      [:NM_TYPES,       0,            :dtypes,      0,      :none,         nil]
  ].map { |d| DTypeInfo.new(*d) }

  DTYPES_ASSIGN = {
      :complex => { # Assign a complex to:
          :complex  => lambda {|l,r| "((#{l}*)p1)->r = ((#{r}*)p2)->r; ((#{l}*)p1)->i = ((#{r}*)p2)->i;" },
          :float    => lambda {|l,r| "*(#{l}*)p1 = ((#{r}*)p2)->r;" },
          :int      => lambda {|l,r| "*(#{l}*)p1 = ((#{r}*)p2)->r;" },
          :rational => lambda {|l,r| "rb_raise(rb_eNotImpError, \"I don't know how to assign a complex to a rational\");"  },
          :value    => lambda {|l,r| "*(VALUE*)p1 = rb_complex_new(((#{r}*)p2)->r, ((#{r}*)p2)->i);" },
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
          :value    => lambda {|l,r| "*(VALUE*)p1 = rb_rational_new(((#{r}*)p2)->n, ((#{r}*)p2)->d);" }
      },
      :value => {
          :complex  => lambda {|l,r| "((#{l}*)p1)->r = NUM2REAL(*(VALUE*)p2); ((#{l}*)p1)->i = NUM2IMAG(*(VALUE*)p2);" },
          :float    => lambda {|l,r| "*(#{l}*)p1 = NUM2DBL(*(VALUE*)p2);"},
          :int      => lambda {|l,r| "*(#{l}*)p1 = NUM2DBL(*(VALUE*)p2);"},
          :rational => lambda {|l,r| "((#{l}*)p1)->n = NUM2NUMER(*(VALUE*)p2); ((#{l}*)p1)->d = NUM2DENOM(*(VALUE*)p2);" },
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


    def dtypes_set_function_ident dtype_i, dtype_j
      dtype_i[:enum] == :NM_NONE || dtype_j[:enum] == :NM_NONE ? "TypeErr" : "Set_#{dtype_i[:id]}_#{dtype_j[:id]}"
    end

    def dtypes_assign lhs, rhs
      Generator::DTYPES_ASSIGN[ rhs.type ][ lhs.type ].call( lhs.sizeof, rhs.sizeof )
    end


    # Declare a set function for a pair of dtypes
    def dtypes_set_function dtype_i, dtype_j
      str = <<SETFN
static void #{dtypes_set_function_ident(dtype_i, dtype_j)}(size_t n, char* p1, size_t i1, char* p2, size_t i2) {
  for (; n > 0; --n) {
    #{dtypes_assign(dtype_i, dtype_j)}
    p1 += i1; p2 += i2;
  }
}

SETFN
    end


    def dtypes_set_functions_matrix
      ary = []
      DTYPES.each do |i|
        next if i[:enum] == :NM_TYPES
        bry = []
        DTYPES.each do |j|
          next if j[:enum] == :NM_TYPES
          bry << dtypes_set_function_ident(i,j)
        end
        ary << "{ " + bry.join(", ") + " }"
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
      ary << decl("nm_setfunc_t SetFuncs =", dtypes_set_functions_matrix)

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
      end
    end


    # Read templates given by +names+ from <tt>SRC_DIR/relative_path</tt>, and output them to a filename described by
    # +output_name+.
    #
    # == Example
    #
    #    make_templated_c './smmp', 'header', %w{numbmm transp bstoy ytobs}, "smmp1.c"
    #
    # TODO: index dtype is unsigned!
    # That means instead of int8_t, we should be doing uint8_t. But can't always do that because the Fortran code
    # occasionally starts at index -1. Stupid Fortran! Someone needs to go through and fix the code by hand.
    #
    # TODO: Write tests to confirm that the signedness isn't screwing up stuff.
    #
    # TODO: Make templates work with complex and rational types too.
    #
    def make_templated_c relative_path, header_name, names, output_name, subs=1

      # First print the header once
      `cat #{Dir.pwd}/../../../../#{SRC_DIR}/#{relative_path}/#{header_name}.template.c > ./#{output_name}` unless header_name.nil?

      DTYPES.each do |index_dtype|
        next unless [:NM_INT8, :NM_INT16, :NM_INT32, :NM_INT64].include?(index_dtype.enum)

        if subs == 1
          names.each do |name|
            sub_int relative_path, name, output_name, index_dtype
          end
        else
          DTYPES.each do |dtype|
            next unless [:NM_FLOAT32, :NM_FLOAT64].include?(dtype.enum)
            names.each do |name|
              sub_int_real relative_path, name, output_name, index_dtype, dtype
            end
          end
        end

      end
    end

    def sub_int_real relative_path, name, output_name, int_dtype, real_dtype
      cmd = ["#{Dir.pwd}/../../../../#{SRC_DIR}/#{relative_path}/#{name}.template.c",
             "sed s/%%INT_ABBREV%%/#{int_dtype.id}/g",
             "sed s/%%INT%%/#{int_dtype.sizeof}/g",
             "sed s/%%REAL_ABBREV%%/#{real_dtype.id}/g",
             "sed s/%%REAL%%/#{real_dtype.sizeof}/g" ]

      raise(ArgumentError, "no such template '#{name}'; looked at #{cmd[0]}") unless File.exist?(cmd[0])

      `cat #{cmd.join(' | ')} >> ./#{output_name}`
    end


    def sub_int relative_path, name, output_name, dtype
      cmd = ["#{Dir.pwd}/../../../../#{SRC_DIR}/#{relative_path}/#{name}.template.c",
             "sed s/%%INT_ABBREV%%/#{dtype.id}/g",
             "sed s/%%INT%%/#{dtype.sizeof}/g"]

      raise(ArgumentError, "no such template '#{name}'; looked at #{cmd[0]}") unless File.exist?(cmd[0])

      `cat #{cmd.join(' | ')} >> ./#{output_name}`
    end


  end
end

Generator.make_dtypes_h
Generator.make_dtypes_c
Generator.make_dfuncs_c
Generator.make_templated_c './smmp', 'blas_header', ['blas1'], 'blas1.c', 1 # 1-type interface functions for SMMP
Generator.make_templated_c './smmp', 'blas_header', ['blas2'], 'blas2.c', 2 # 2-type interface functions for SMMP
Generator.make_templated_c './smmp', 'smmp_header', ['symbmm'], 'smmp1.c', 1 # 1-type SMMP functions from Fortran
Generator.make_templated_c './smmp', 'smmp_header', ['numbmm', 'transp'], 'smmp2.c', 2 # 2-type SMMP functions from Fortran