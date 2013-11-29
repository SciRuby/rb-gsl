# -*- encoding: utf-8 -*-
require 'date'

Gem::Specification.new do |s|
  s.name = 'rb-gsl'
  s.version = '1.16.0.6'
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
  s.required_ruby_version = '>= 1.9.3'
  s.requirements = ['GSL (http://www.gnu.org/software/gsl/)']

  s.post_install_message = 'rb-gsl has been replaced by gsl'
  s.add_runtime_dependency 'gsl'
end
