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
  self.require_ruby_version ">=1.9"
  self.developer('John Woods', 'john.o.woods@gmail.com')
  self.post_install_message = <<-EOF
***********************************************************
Welcome to SciRuby: Tools for Scientific Computing in Ruby!

                     *** WARNING ***
Please be aware that NMatrix is in ALPHA status. If you're
thinking of using NMatrix to write mission critical code,
such as for driving a car or flying a space shuttle, you
may wish to choose other software (for now).

NMatrix requires a C compiler, and has been tested only
with GCC 4.6.1. We are happy to accept contributions
which improve the portability of this project.

Also required is ATLAS. Most Linux distributions and Mac
versions include ATLAS, but you may wish to compile it
yourself.

More explicit instructions for NMatrix and SciRuby should
be available on the SciRuby website, sciruby.com, or
through our mailing list (which can be found on our web-
site).

Thanks for trying out NMatrix! Happy coding!

***********************************************************
EOF
  #self.need_rdoc = false
  self.readme_file = 'README.rdoc'
  # self.rubyforge_name = 'nmatrixx' # if different than 'nmatrix'
end

RSpec::Core::RakeTask.new(:spec)

task :console do |task|
  cmd = [ 'irb', "-r './lib/nmatrix.rb'" ]
  run *cmd
end

#namespace :console do
#  CONSOLE_CMD = ['irb', "-r './lib/nmatrix.rb'"]
#  desc "Run console under GDB."
#  task :gdb => [ :compile ] do |task|
#          cmd = [ 'gdb' ] + GDB_OPTIONS
#          cmd += [ '--args' ]
#          cmd += CONSOLE_CMD
#          run( *cmd )
#  end
#
#  desc "Run console under Valgrind."
#  task :valgrind => [ :compile ] do |task|
#          cmd = [ 'valgrind' ] + VALGRIND_OPTIONS
#          cmd += CONSOLE_CMD
#          run( *cmd )
#  end
#end

task :default => :spec

def run *cmd
  sh(cmd.join(" "))
end

namespace :spec do
  # partial-loads-ok and undef-value-errors necessary to ignore
  # spurious (and eminently ignorable) warnings from the ruby
  # interpreter

  RSPEC_CMD = [ 'ruby', '-S', 'rspec', '-Ilib:ext', SPECDIR ]

  #desc "Run the spec for generator.rb"
  #task :generator do |task|
  #  run 'rspec spec/generator_spec.rb'
  #end

  desc "Run specs under GDB."
  task :gdb => [ :compile ] do |task|
          cmd = [ 'gdb' ] + GDB_OPTIONS
          cmd += [ '--args' ]
          cmd += RSPEC_CMD
          run( *cmd )
  end

  desc "Run specs under Valgrind."
  task :valgrind => [ :compile ] do |task|
          cmd = [ 'valgrind' ] + VALGRIND_OPTIONS
          cmd += RSPEC_CMD
          run( *cmd )
  end
end


namespace :clean do
  task :so do |task|
    tmp_path = "tmp/#{RUBY_PLATFORM}/nmatrix/#{RUBY_VERSION}"
    chdir tmp_path do
      if RUBY_PLATFORM =~ /mswin/
        `nmake soclean`
      else
        mkcmd = ENV['MAKE'] || %w[gmake make].find { |c| system("#{c} -v >> /dev/null 2>&1") }
        `#{mkcmd} soclean`
      end
    end
  end
end

# vim: syntax=ruby
