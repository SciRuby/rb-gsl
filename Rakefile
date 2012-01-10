# -*- ruby -*-

require 'rubygems'
require 'hoe'
require 'rspec/core/rake_task'

Hoe.plugin :compiler
Hoe.plugin :bundler
Hoe.plugin :git
# Hoe.plugin :compiler
# Hoe.plugin :gem_prelude_sucks
# Hoe.plugin :inline
# Hoe.plugin :racc
# Hoe.plugin :rubyforge

Hoe.spec 'nmatrix' do
  developer('John Woods', 'john.o.woods@gmail.com')

  # self.rubyforge_name = 'nmatrixx' # if different than 'nmatrix'
end

RSpec::Core::RakeTask.new(:spec)

task :default => :spec

# vim: syntax=ruby
