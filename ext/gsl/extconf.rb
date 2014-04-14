require 'mkmf'

def gsl_config_arg(arg)
  yield arg_config("--with-gsl-#{arg}") {
    sh = 'sh ' if RUBY_PLATFORM =~ /mingw/
    IO.popen("#{sh}gsl-config --#{arg}") { |f| f.gets.chomp }
  }, lambda { |val| puts "checking gsl #{arg}... #{val}"; val }
rescue => err
  abort "*** ERROR: missing required library to compile this module: #{err}"
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

def gsl_dir_config(target)
  dir_config(target, $sitearchdir, $sitearchdir)
end

$CFLAGS += ' -Wall -Iinclude'

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

  if enable_config('atlas')
    dir_config('cblas')
    dir_config('atlas')

    if have_library('cblas') && have_library('atlas')
      libs.gsub!('-lgslcblas', '-lcblas -latlas')
    end
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

narray = gsl_dir_config('narray')

begin
  require 'rubygems'

  if spec = Gem::Specification.find_by_path('narray.h')
    $LOCAL_LIBS = "-L#{File.join(narray = spec.full_gem_path, 'src')} " + $LOCAL_LIBS
    $CPPFLAGS   = "-I#{File.join(narray, spec.require_path)} "          + $CPPFLAGS
  end
rescue LoadError
end

have_header('narray.h')
have_library('narray') if narray && RUBY_PLATFORM =~ /cygwin|mingw/

if !arg_config('--disable-tamu-anova') && gsl_dir_config('tamu_anova')
  gsl_have_header('tamuanova', 'tamu_anova/tamu_anova.h')
end

$objs = Dir["#{File.dirname(__FILE__)}/*.c"].map { |f| File.basename(f, '.c') << '.o' }.
  sort - %w[block matrix poly tensor vector].map { |f| "#{f}_source.o" }

create_makefile('gsl/gsl_native')
