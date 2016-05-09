require 'mkmf'

def gsl_config_arg(arg)
  yield arg_config("--with-gsl-#{arg}") {
    sh = 'sh ' if RUBY_PLATFORM =~ /mingw/
    IO.popen("#{sh}gsl-config --#{arg}") { |f| f.gets.chomp }
  }, lambda { |val| puts "checking gsl #{arg}... #{val}"; val }
rescue => err
  abort "*** ERROR: missing required library to compile this module: #{err}"
end

# Function derived from NArray's extconf.rb.
def create_conf_h(file) #:nodoc:
  print "creating #{file}\n"
  File.open(file, 'w') do |hfile|
    header_guard = file.upcase.sub(/\s|\./, '_')

    hfile.puts "#ifndef #{header_guard}"
    hfile.puts "#define #{header_guard}"
    hfile.puts

    # FIXME: Find a better way to do this:
    hfile.puts "#define RUBY_2 1" if RUBY_VERSION >= '2.0'

    for line in $defs
      line =~ /^-D(.*)/
      match_data = $1.dup

      if match_data.match(/GSL_VERSION*/)
        hfile.printf "#define #{match_data.to_s.split('=').join(' ')}\n"
      else
        hfile.printf "#define %s 1\n", match_data
      end
    end

    hfile.puts
    hfile.puts "#endif"
  end
end

def gsl_def(const, value = nil)
  value = "=#{value}" if value
  $defs << "-D#{const}#{value}"
end

def gsl_have_header(library, header)
  have_library(library) if have_header(header)
end

def gsl_have_library(func)
  have_func(func) if have_library('gsl', func)
end

def gsl_dir_config(target, idir = nil, ldir = idir)
  dir_config(target, idir || $sitearchdir, ldir || $sitearchdir)
end

def gsl_gem_config(target, dir = 'ext')
  path = begin
    require 'rubygems'

    # For some weird reason finding narray and nmatrix headers works with specific
    # functions only. Might need to fix nmatrix/narray to place headers in correct
    # locations?
    if target == 'narray'
      gem 'narray'
      spec = Gem::Specification.find_by_path("#{target}")
      File.join(spec.full_gem_path, dir) if spec
    else
      gem 'nmatrix'
      spec = Gem::Specification.find_all_by_name("#{target}").compact
      File.join(spec[0].require_path)
    end

  rescue LoadError
  end

  gsl_dir_config(target, path)

  $LOCAL_LIBS += " -l:#{target}.so" if arg_config("--force-link-#{target}")   ||
                                       $CFLAGS.include?('-Wl,--no-undefined') ||
                                       $LDFLAGS.include?('-Wl,--no-undefined')
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

  raise 'Ruby/GSL requires gsl-1.15 or later.' unless later['1.15']

  %w[1.15 1.16 2.0 2.1].each { |v| later[v] }
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

external_libs = []
external_libs << 'narray' if ENV['NARRAY']
external_libs << 'nmatrix' if ENV['NMATRIX']

external_libs.each do |library|
  gsl_gem_config(library)
  have_header("#{library}.h")
  have_header("nmatrix_config.h") if library == 'nmatrix'
  have_library(library) if RUBY_PLATFORM =~ /cygwin|mingw/  
end

unless arg_config('--disable-tamu-anova')
  gsl_dir_config('tamu_anova')
  gsl_have_header('tamuanova', 'tamu_anova/tamu_anova.h')
end

have_struct_member('gsl_multifit_fdfsolver', 'J', 'gsl/gsl_multifit_nlin.h')
have_func('gsl_sf_mathieu_a_e',  'gsl/gsl_sf_mathieu.h');
have_func('gsl_sf_mathieu_b_e',  'gsl/gsl_sf_mathieu.h');
have_func('gsl_sf_mathieu_ce_e', 'gsl/gsl_sf_mathieu.h');
have_func('gsl_sf_mathieu_se_e', 'gsl/gsl_sf_mathieu.h');
have_func('gsl_sf_mathieu_Mc_e', 'gsl/gsl_sf_mathieu.h');
have_func('gsl_sf_mathieu_Ms_e', 'gsl/gsl_sf_mathieu.h');

create_conf_h('gsl_config.h')
create_makefile('gsl_native')
