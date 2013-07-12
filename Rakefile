require 'rubygems'
require 'rubygems/package_task'
require 'rdoc/task'

RB_GSL_VERSION = File.readlines('VERSION')[0].chomp

spec = gemspec = eval(IO.read("gsl.gemspec"))

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
