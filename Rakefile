require File.expand_path(%q{../lib/gsl/version}, __FILE__)

require 'rubygems'
require 'rubygems/package_task'
require 'bundler'

Bundler::GemHelper.install_tasks

begin
  Bundler.setup(:default, :development)
rescue Bundler::BundlerError => e
  $stderr.puts e.message
  $stderr.puts "Run `bundle install` to install missing gems"
  exit e.status_code
end

require 'rake'
require 'rake/testtask'

Rake::TestTask.new do |t|
  t.libs << "test"
  t.libs << "test/gsl"
  t.test_files = FileList['test/*.rb', 'test/gsl/*.rb']
end

gemspec = eval(IO.read("rb-gsl.gemspec"))
Gem::PackageTask.new(gemspec).define

task default: :test
