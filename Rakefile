# -*- ruby -*-

require 'rubygems'
require 'hoe'
require 'pathname'
require 'rspec/core/rake_task'

Hoe.plugin :compiler
Hoe.plugin :bundler
Hoe.plugin :git
# Hoe.plugin :compiler
# Hoe.plugin :gem_prelude_sucks
# Hoe.plugin :inline
# Hoe.plugin :racc
# Hoe.plugin :rubyforge

BASEDIR = Pathname( __FILE__ ).dirname.relative_path_from( Pathname.pwd )
SPECDIR = BASEDIR + 'spec'

VALGRIND_OPTIONS = [
        "--num-callers=50",
        "--error-limit=no",
        "--partial-loads-ok=yes",
        "--undef-value-errors=no",
]
VALGRIND_MEMORYFILL_OPTIONS = [
        "--freelist-vol=100000000",
        "--malloc-fill=6D",
        "--free-fill=66 ",
]

GDB_OPTIONS = []


h = Hoe.spec 'nmatrix' do
  developer('John Woods', 'john.o.woods@gmail.com')
  # self.rubyforge_name = 'nmatrixx' # if different than 'nmatrix'
end

RSpec::Core::RakeTask.new(:spec)

task :console do |task|
  cmd = [ 'irb', "-r './lib/nmatrix.rb'" ]
  run *cmd
end

task :default => :spec

def run *cmd
  sh(cmd.join(" "))
end

namespace :spec do
  # partial-loads-ok and undef-value-errors necessary to ignore
  # spurious (and eminently ignorable) warnings from the ruby
  # interpreter

  RSPEC_CMD = [ 'ruby', '-S', 'rspec', '-Ilib:ext', SPECDIR ]

  desc "Run the specs under GDB."
  task :gdb => [ :compile ] do |task|
          cmd = [ 'gdb' ] + GDB_OPTIONS
          cmd += [ '--args' ]
          cmd += RSPEC_CMD
          run( *cmd )
  end

  desc "Run the specs under Valgrind."
  task :valgrind do
          cmd = [ 'valgrind' ] + VALGRIND_OPTIONS
          cmd += RSPEC_CMD
          run( *cmd )
  end
end

# vim: syntax=ruby
