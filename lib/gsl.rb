begin
  require 'narray'
rescue LoadError
end

begin
  require "gsl/#{RUBY_VERSION[/\d+.\d+/]}/gsl_native"
rescue LoadError
  require 'gsl/gsl_native'
end

require 'gsl/version'
require 'gsl/oper'
