# -*- encoding: utf-8 -*-
require File.dirname(__FILE__) + '/lib/gsl/version'
require 'date'

Gem::Specification.new do |s|
  s.name = 'gsl'
  s.version = GSL::RUBY_GSL_VERSION
  s.date = Date.today.to_s

  s.require_paths = %w(lib)
  s.authors = ['Yoshiki Tsunesada', 'David MacMahon', 'Jens Wille', 'Daniel Mendler']
  s.summary = 'Ruby interface to the GNU Scientific Library'
  s.description = 'Ruby/GSL is a Ruby interface to the GNU Scientific Library, for numerical computing with Ruby'
  s.email = 'mail@daniel-mendler.de'
  s.extensions  = Dir['ext/**/extconf.rb']
  s.extra_rdoc_files = Dir['**/*.rdoc']
  s.files = `git ls-files`.split("\n")
  s.homepage = 'http://github.com/SciRuby/rb-gsl'
  s.licenses = ['GPL-2.0']
  s.rdoc_options = ['--title', "Ruby/GSL (#{GSL::RUBY_GSL_VERSION})", '--charset', 'UTF-8', '--line-numbers', '--all', '--main', 'index.rdoc', '--root', 'rdoc']
  s.required_ruby_version = '>= 1.9.3'
  s.requirements = ['GSL (http://www.gnu.org/software/gsl/)']

  s.post_install_message = %{
    #{s.name} can be installed with or without narray support. Please install
    narray before and reinstall #{s.name} if it is missing.

    #{s.name} is also now compatible with NMatrix. Please install nmatrix before
    installing #{s.name}.
  }

  s.add_development_dependency 'rake-compiler', '>= 0'
  s.add_development_dependency 'rake', '>= 0'
  s.add_development_dependency 'test-unit', '>= 0'
  s.add_development_dependency 'bundler', '~> 1.11'
end
