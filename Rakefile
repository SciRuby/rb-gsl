require 'rubygems'
require 'rubygems/package_task'
require 'rdoc/task'

RB_GSL_VERSION = File.readlines('VERSION')[0].chomp

spec = Gem::Specification.new do |s|
  # Basics
  s.name = 'rb-gsl'
  s.version = RB_GSL_VERSION
  s.summary = 'Ruby interface to GNU Scientific Library'
  s.description = 'Ruby/GSL is a Ruby interface to the GNU Scientific Library, for numerical computing with Ruby'
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
    'AUTHORS',
    'COPYING',
    'ChangeLog',
    'README',
    'Rakefile',
    'setup.rb',
    'THANKS',
    'VERSION',
    'examples/**/*',
    'ext/extconf.rb',
    'ext/*.c',
    'lib/**/*',
    'include/*',
    'rdoc/*',
    'tests/**/*'
  ].to_a
  s.require_paths = ['lib', 'lib/gsl', 'lib/ool', 'ext']
  #s.autorequire = nil
  #s.bindir = 'bin'
  #s.executables = []
  #s.default_executable = nil

  # C compilation
  s.extensions = %w[ ext/extconf.rb ]

  # Documentation
  s.has_rdoc = true
  s.rdoc_options = [
    '--title', 'Ruby/GSL',
    '--main', 'rdoc/index.rdoc',
    '--exclude', 'ext/',
    '--exclude', 'include/',
    '--exclude', 'lib/',
  ]
  s.extra_rdoc_files = FileList['rdoc/*'].to_a

  # Testing TODO
  #s.test_files = []
end

Rake::PackageTask.new('rb-gsl', RB_GSL_VERSION) do |pkg|
  pkg.need_zip = true
  pkg.need_tar = true
  pkg.package_files = spec.files
end

Gem::PackageTask.new(spec) do |pkg|
  pkg.need_zip = false
  pkg.need_tar = false
end

task :default => [:package, :gem]

# --------------------------------------------------------------------
# Create a task to build the RDOC documentation tree.

desc "Create the RDoc html files"
Rake::RDocTask.new("rdoc") { |rdoc|
  rdoc.rdoc_dir = 'html'
  rdoc.title    = 'Ruby/GSL'
  rdoc.main     = 'rdoc/index.rdoc'
  rdoc.options << '--exclude' << 'ext/'
  rdoc.options << '--exclude' << 'include/'
  rdoc.options << '--exclude' << 'lib/'
  rdoc.rdoc_files.include('rdoc/*.rdoc')
}

desc "Publish the RDoc files on RubyForge"
task :pub_rdoc => ["html/index.html"] do
  mkdir_p "emptydir"
  sh "scp -rq html/* www.rubyforge.org:/var/www/gforge-projects/rb-gsl/."
  rm_r "emptydir"
end
