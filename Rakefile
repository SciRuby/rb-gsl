require File.expand_path(%q{../lib/gsl/version}, __FILE__)

begin
  require 'hen'

  note = '[Ruby 2.x and GSL 1.16 compatible fork]'

  Hen.lay! {{
    :gem => {
      :name          => %q{rb-gsl},
      :version       => GSL::RB_GSL_VERSION,
      :summary       => %Q{Ruby interface to the GNU Scientific Library #{note}},
      :description   => %Q{Ruby/GSL is a Ruby interface to the GNU Scientific Library, for numerical computing with Ruby #{note}},
      :authors       => ['Yoshiki Tsunesada', 'David MacMahon', 'Jens Wille'],
      :email         => %q{jens.wille@gmail.com},
      :license       => %q{GPL-2.0},
      :homepage      => :blackwinter,
      :dependencies  => [['narray', '>= 0.5.9']],
      :requirements  => ['GSL (http://www.gnu.org/software/gsl/)'],

      :extra_files => FileList['examples/**/*', 'rdoc/*'].to_a,

      :extension => { :cross_compile => false },

      :required_ruby_version => '>= 1.8.7'
    },
    :rdoc => {
      :title      => 'Ruby/GSL{version: (v%s)}',
      :rdoc_files => FileList['rdoc/*'].to_a,
      :main       => 'index.rdoc',
      :root       => 'rdoc'
    },
    :test => {
      :libs => %w[lib test]
    }
  }}
rescue LoadError => err
  warn "Please install the `hen' gem. (#{err})"
end
