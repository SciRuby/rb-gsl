# Ruby/GSL, a Ruby interface to GSL (GNU Scientific library)

Permission is granted to copy, distribute and/or modify this document under
the terms of the GNU Free Documentation License.

## Description

Ruby/GSL is a Ruby interface to the [GNU Scientific Library](https://gnu.org/software/gsl/)
(GSL), for numerical computing with [Ruby](http://www.ruby-lang.org/).


## Installation

Ruby/GSL may be installed as a Ruby Gem by simply running

  gem install rb-gsl

Note that the GSL libraries must already be installed before Ruby/GSL
can be installed:

Debian/Ubuntu: +libgsl0-dev+
Fedora/SuSE:   +gsl-devel+
Gentoo:        +sci-libs/gsl+
OS X:          <tt>brew install gsl</tt>

It is recommended to install the [GNU plotutils](https://gnu.org/software/plotutils/plotutils.html)
package. Some of the example scripts in the +examples/+ directory use the
+graph+ utility included in the package to plot the results. Windows cygwin
binaries of <tt>GNU plotutils</tt> and related packages are available
[here](http://gnuwin32.sourceforge.net/packages/plotutils.htm).


## Reference

The [Ruby/GSL reference manual](link:rdoc/ref_rdoc.html) follows and borrows
large parts of the GSL reference manual.


## Examples

See scripts in +examples/+ and +test/+ directories.


## Related Projects

* [ruby-gsl](http://ruby-gsl.sourceforge.net/): Another Ruby binding,
  developed by Arno Erpenbeck.


## Licence

Ruby/GSL is free software: you can redistribute it and/or modify it
under the terms of the GNU General Public License.
This library is distributed in the hope that it will be useful, but
WITHOUT ANY WARRANTY.


## Bug Reports

Any bug reports are welcome. If you encounter bugs in Ruby/GSL, please
report them on GitHub(https://github.com/SciRuby/rb-gsl/issues).


## Links

Documentation: https://sciruby.github.com/rb-gsl
Source code:   https://github.com/SciRuby/rb-gsl
RubyGem:       https://rubygems.org/gems/rb-gsl
Bug tracker:   https://github.com/SciRuby/rb-gsl/issues
Travis CI:     https://travis-ci.org/SciRuby/rb-gsl


## Authors

* Yoshiki Tsunesada <y-tsunesada at mm dot em-net dot ne dot jp> (July, 2004)
* David MacMahon <davidm@astro.berkeley.edu> (November, 2010)
* Jens Wille <mailto:jens.wille@gmail.com> (November, 2013)
