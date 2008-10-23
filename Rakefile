require 'rubygems'
require 'rake/gempackagetask'

spec = Gem::Specification.new do |s|
  # Basics
  s.name = 'gsl'
  s.version = File.readlines('VERSION')[0].chomp
  s.summary = 'Ruby interface to GSL'
  s.description = 'RubyGSL is a Ruby interface to the GNU Scientific Library, for numerical computing with Ruby'
  #s.platform = Gem::Platform::Ruby
  s.required_ruby_version = '>= 1.8.1'
  s.requirements << 'GSL (http://www.gnu.org/software/gsl/)'
  # plotlib?
  s.add_dependency('narray', '>= 0.5.9')

  # About
  s.authors = ['Yoshiki Tsunesada', 'David MacMahon']
  s.email = 'y-tsunesada@mm.em-net.ne.jp'
  s.homepage = 'http://rb-gsl.rubyforge.org/'
  s.rubyforge_project = 'rb-gsl' 

  # Files, Libraries, and Extensions
  s.files = FileList[
    'README',
    'VERSION',
    'Rakefile',
    'ext/*',
    'lib/**/*',
    'include/*'
  ].to_a
  s.require_paths = ['lib', 'lib/gsl', 'lib/ool', 'ext']
  #s.autorequire = nil
  #s.bindir = 'bin'
  #s.executables = []
  #s.default_executable = nil

  # C compilation
  s.extensions = %w[ ext/extconf.rb ]

  # Documentation TODO
  #s.rdoc_options = []
  #s.has_rdoc = false
  #s.extra_rdoc_files = []

  # Testing TODO
  #s.test_files = []
end

Rake::GemPackageTask.new(spec) do |pkg|
  pkg.need_zip = true
  pkg.need_tar = true
end

task :default => :gem
