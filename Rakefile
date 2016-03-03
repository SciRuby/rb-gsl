require File.expand_path(%q{../lib/gsl/version}, __FILE__)

require 'bundler/setup'
require 'rubygems/package_task'
require 'rake/extensiontask'
require 'rake/testtask'

Bundler::GemHelper.install_tasks

Rake::TestTask.new do |t|
  t.libs << 'test'
  t.libs << 'test/gsl'
  file_list = [ 'test/*.rb', 'test/gsl/*.rb']  
  if ENV['NMATRIX']
    t.libs    << 'test/gsl/nmatrix_tests' 
    file_list << 'test/gsl/nmatrix_tests/*_test.rb'
  end
 
  t.test_files = FileList[*file_list]
end

spec = eval(IO.read('gsl.gemspec'))
Gem::PackageTask.new(spec).define
Rake::ExtensionTask.new(:gsl_native, spec)

task :default => [:compile, :test]
