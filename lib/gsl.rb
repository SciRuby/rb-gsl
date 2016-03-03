begin
  require 'narray'
rescue LoadError
end

begin
  require 'nmatrix/nmatrix'
rescue LoadError
end

require 'gsl_native'
require 'gsl/version'
require 'gsl/oper'
