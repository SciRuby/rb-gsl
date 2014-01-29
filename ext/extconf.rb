require 'mkmf'

def gsl_config_arg(arg)
  yield arg_config("--with-gsl-#{arg}") {
    sh = 'sh ' if RUBY_PLATFORM =~ /mingw/
    IO.popen("#{sh}gsl-config --#{arg}") { |f| f.gets.chomp }
  }, lambda { |val| puts "checking gsl #{arg}... #{val}"; val }
rescue => err
  abort "*** ERROR: missing required library to compile this module: #{err}"
end

def gsl_file_path(path = '.')
  File.expand_path(path, File.dirname(__FILE__))
end

def gsl_def(const, value = nil)
  $defs << "-D#{const}#{"=#{value}" if value}"
end

def gsl_have_header(library, header)
  have_library(library) if have_header(header)
end

def gsl_have_library(func)
  have_func(func) if have_library('gsl', func)
end

$CFLAGS += " -Wall -I#{gsl_file_path('../include')}"

gsl_config_arg(:version) { |version, check|
  gsl_def(:GSL_VERSION, check[version])

  ver = version.split('.').map { |x| x.to_i }

  later = lambda { |other|
    ary = other.split('.').map { |x| x.to_i }

    gte = ver[0] > ary[0] ? true  :
          ver[0] < ary[0] ? false :
          ver[1] > ary[1] ? true  :
          ver[1] < ary[1] ? false :
          ver.size < ary.size ? false :
          ver.size == 3 && ary.size == 3 ? ver[2] >= ary[2] : true

    ary.pop && ary[-1] += 1 if ary.last == 90

    gte && gsl_def("GSL_#{ary.join('_')}_LATER")
  }

  raise 'Ruby/GSL requires gsl-0.9.4 or later.' unless later['0.9.4']

  gsl_def(:GSL_1_4_9_LATER) if later['1.4.90']

  %w[
    1.0 1.1 1.1.1 1.2 1.3 1.4 1.5.90 1.7.90
    1.8.90 1.9.90 1.11 1.12.90 1.14 1.15
  ].each { |v| later[v] }
}

gsl_config_arg(:cflags) { |cflags, check|
  $CFLAGS += ' ' + check[cflags]
}

gsl_config_arg(:libs) { |libs, check|
  libs.tr!(File::PATH_SEPARATOR, ' ')

  dir_config('cblas')
  dir_config('atlas')

  if have_library('cblas') && have_library('atlas')
    libs.gsub!('-lgslcblas', '-lcblas -latlas')
  end

  $LOCAL_LIBS += ' ' + check[libs]
}

have_header('ruby/io.h')
have_func('round')

%w[alf qrngextra rngextra tensor].each { |library|
  gsl_have_header(library, "#{library}/#{library}.h")
}

gsl_have_header('bundle_method', 'gsl/gsl_multimin_fsdf.h')
gsl_have_header('cqp',           'gsl/gsl_cqp.h')
gsl_have_header('jacobi',        'jacobi.h')
gsl_have_header('ndlinear',      'ndlinear/gsl_multifit_ndlinear.h')
gsl_have_header('ool',           'ool/ool_version.h')

gsl_have_library('gsl_eigen_francis')
gsl_have_library('gsl_poly_solve_quartic')

gsl_def(:HAVE_GNU_GRAPH) if find_executable('graph')

####

#narray_config = dir_config("narray")
narray_config = dir_config('narray',$sitearchdir,$sitearchdir)
# Try to find narray with RubyGems
begin
  require 'rubygems'
  na_gemspec=Gem::Specification.find_by_path('narray.h')
  if na_gemspec
    narray_config = na_gemspec.full_gem_path
    $CPPFLAGS = " -I#{File.join(narray_config, na_gemspec.require_path)} "+$CPPFLAGS
    $LOCAL_LIBS = " -L#{File.join(narray_config, 'src')}" + $LOCAL_LIBS
  end
rescue LoadError
end
have_narray_h = have_header("narray.h")
if narray_config
  if RUBY_PLATFORM =~ /cygwin|mingw/
#    have_library("narray") || raise("ERROR: narray import library is not found") 
  have_library("narray")
  end
end

unless arg_config('--disable-tamu-anova')
tamu_anova_config = dir_config('tamu_anova',$sitearchdir,$sitearchdir)
have_header("tamu_anova/tamu_anova.h")
if tamu_anova_config
  have_library("tamuanova")
#  if RUBY_PLATFORM =~ /cygwin|mingw/
#    have_library("tamuanova") || raise("ERROR: tamu_anova import library is not found")
#  end
end
end

File.open(gsl_file_path("../lib/gsl.rb"), "w") do |file|
  file.print("require('gsl/version')\n")
  if have_narray_h
    file.print("require('narray')\n")
  end
#  file.print("require('rb_gsl')\ninclude GSL\n")
  file.print("require('rb_gsl')\n")  
  file.print("require('gsl/oper.rb')\n")
end

File.open(gsl_file_path("../lib/rbgsl.rb"), "w") do |file|
  if have_narray_h
    file.print("require('narray')\n")
  end
  file.print("require('rb_gsl')\n")
  file.print("require('gsl/oper.rb')\n")
end

####

Dir.chdir(gsl_file_path) {
  $objs = Dir['*.c'].map { |f| f.sub('.c', '.o') }.sort - %w[
    block_source.o matrix_source.o poly_source.o tensor_source.o vector_source.o
  ]
}

create_makefile('rb_gsl')
