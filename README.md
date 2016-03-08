# Ruby/GSL, a Ruby interface to GSL (GNU Scientific library)

Permission is granted to copy, distribute and/or modify this document under
the terms of the GNU Free Documentation License.

## Description

Ruby/GSL is a Ruby interface to the [GNU Scientific Library](https://gnu.org/software/gsl/)
(GSL), for numerical computing with [Ruby](http://www.ruby-lang.org/).

Ruby/GSL is compatible with GSL versions upto 2.1.

## Usage with GSL 2.1

As of this release, GSL 2.1 has not made it's way into the Debian stable repositories. Hence, after compiling GSL 2.1 from source, you will need to set the installation location in your `LD_LIBRARY_PATH` variable. After following standard GSL 2.1 installation procedures, you should do:
  export LD_LIBRARY_PATH=/usr/local/lib

The need to do this should not arise if GSL has been installed from `apt-get`.

## Installation

Ruby/GSL may be installed as a Ruby Gem by simply running

  gem install gsl

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

## NMatrix and NArray usage

Ruby/GSL works with [NMatrix](https://github.com/SciRuby/nmatrix) and [NArray](https://github.com/masa16/narray) for a variety of methods. See the docs for a detailed list.

### Basic Installation

In order to use rb-gsl with NMatrix you must first set the `NMATRIX` environment variable and then install rb-gsl:
  gem install nmatrix
  export NMATRIX=1
  gem install rb-gsl

This will compile rb-gsl with NMatrix specific functions.

For using rb-gsl with NArray:
  gem install narray
  export NARRAY=1
  gem install rb-gsl

Note that setting both `NMATRIX` and `NARRAY` variables will lead to undefined behaviour. Only one can be used at a time.

### NMatrix basic usage

Convert an NMatrix 1D vector to GSL::Vector:
``` ruby
require 'gsl'

nm = NMatrix.new([5], [1,2,3,4,5], dtype: :float64)
#=> [1.0, 2.0, 3.0, 4.0, 5.0]
nm.to_gslv
# => GSL::Vector
# [ 1.000e+00 2.000e+00 3.000e+00 4.000e+00 5.000e+00 ]
```

Convert an integer 2D NMatrix to GSL::Matrix::Int:
``` ruby

require 'gsl'
nm = NMatrix.new([3,3], [2]*9, dtype: :int32)
#=> 
#[
#  [2, 2, 2]   [2, 2, 2]   [2, 2, 2] ]
nm.to_gslm
#=> GSL::Matrix::Int
#[ 2 2 2 
#  2 2 2 
#  2 2 2 ]
```

Convert GSL::Vector to 1D NMatrix:
``` ruby
g = GSL::Vector.alloc(1,2,3,4)
# => GSL::Vector
# [ 1.000e+00 2.000e+00 3.000e+00 4.000e+00 ]
g.to_nm
# => [1.0, 2.0, 3.0, 4.0]
```

`to_nm` can be used on all sorts of `GSL::Vector` and `GSL::Matrix` objects to convert them to NMatrix.

For a detailed list of methods that are compatible with NMatrix, see 'nmatrix' in the [docs](https://sciruby.github.com/rb-gsl).

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

## Testing and Contributing

If you wish to make contributions, run the following commands to clone and test the gem on your local machine:
  git clone https://github.com/SciRuby/rb-gsl.git
  cd rb-gsl
  bash test.sh

This will run tests with and without NMatrix/NArray.

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
