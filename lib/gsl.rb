begin
  require 'narray' if ENV['NARRAY']
rescue LoadError
end

begin
  require 'nmatrix/nmatrix' if ENV['NMATRIX']
rescue LoadError
end

require 'gsl_native'
require 'gsl/version'
require 'gsl/oper'
require 'gsl/interp2d_fix' if GSL::GSL_VERSION >= '2.0'