require("rake/testtask")
require("rake/clean")

PROJECT = "Ruby/GSL"
PROJECT_VERSION = "1.10.0"
MY_NAME = "Yoshiki Tsunesada"
MY_EMAIL = "y-tsunesada@mm.em-net.ne.jp"
PROJECT_SUMMARY = "Ruby extension for GSL, GNU Scientific Library"
UNIX_NAME = "rb-gsl"
RUBYFORGE_USER = ENV["RUBYFORGE_USER"] || "ytsunesada"

EXT_DIR = "ext"
HAVE_EXT = File.directory?(EXT_DIR)
EXTCONF_FILES = EXT_DIR + "/extconf.rb"
EXT_SOURCES = FileList["#{EXT_DIR}/*.c"] - ["#{EXT_DIR}/block_source", "#{EXT_DIR}/vector_source", "#{EXT_DIR}/matrix_source", "#{EXT_DIR}/tentor_source","#{EXT_DIR}/poly_source.c"]

CLOBBER.include("#{EXT_DIR}/*.{so,bundle,dll,o}", "#{EXT_DIR}/Makefile")
CLOBBER.include(".config")

CONFIG_OPTS = ENV["CONFIG"]
if HAVE_EXT
  file_create ".config" do
    ruby "setup.rb config #{CONFIG_OPTS}"
  end

  desc "Configure and make extension. " +
    "The CONFIG variable is passed to 'setup.rb config'"
  task "make-ext" => ".config" do
    ruby "setup.rb -q setup"
  end
end

task "default" => ["make-ext"]

