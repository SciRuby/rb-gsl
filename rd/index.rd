=begin
= Ruby/GSL

== Description
((<RubyGSL|URL:http://rubyforge.org/projects/rb-gsl/>)) is a ruby interface to the ((<GSL|URL:http://www.gnu.org/software/gsl/>)) (GNU Scientific Library), for numerical computing with ((<Ruby|URL:http://www.ruby-lang.org/en/>)). 

#Ruby/GSL is developed with Ruby 1.8.3 and GSL 1.7 on 
#((<MacOS X|URL:http://www.apple.com/macosx/overview/>)) 10.3.7 (Darwin 7.7.0).
#This also runs under Linux, and ((<Cygwin|URL:http://cygwin.com/>)).

== Download
  * ((<RubyForge|URL:http://rubyforge.org/frs/?group_id=285>))

== Installation
   (1) Get and install ((<GSL|URL:http://www.gnu.org/software/gsl/#downloading>)). Make sure the command "gsl-config" is in command search path.
   (2) ((<Download|URL:http://rubyforge.org/frs/?group_id=285>)) Ruby/GSL, ungzip and untar the archive (({rb-gsl-xxx.tar.gz})). 
   (3) ((% % cd rb-gsl-xxx/%))
   (4) ((% % ruby setup.rb config%))
   (5) ((% % ruby setup.rb setup%))
   (6) ((% % ruby setup.rb install%)) (as root)

  * It is recommended to install the
    ((<(({GNU plotutils}))|URL:http://www.gnu.org/software/plotutils/plotutils.html>))
    package.
    Some of the example scripts in the 'examples/' directory use the (({graph})) 
    utility included in the package to plot the results. 
    Windows-cygwin binaries of (({GNU plotutils})) and 
    related packages are available from 
    ((<here|URL:http://rustam.uwp.edu/support>)).

== ((<Screenshot|URL:screenshot.html>))

== ((<Reference|URL:ref.html>))

== Examples
See scripts in (({examples/})) and (({tests/})) directories.

== Related Projects
* ((<ruby-gsl:|URL:http://ruby-gsl.sourceforge.net/>))
  Another Ruby biding, developed by Arno Erpenbeck.

== Licence
Ruby/GSL is free software: you can redistribute it and/or modify it
under the terms of the GNU General Public License.
This library is distributed in the hope that it will be useful, but
WITHOUT ANY WARRANTY.

== Bug Report
If you encounter bugs in Ruby/GSL, please e-mail to me, or
submit reports from ((<RubyForge page|URL:http://rubyforge.org/projects/rb-gsl/>)).

== Author
Yoshiki Tsunesada

Jul/2004

This document is generated with ((<RDTool|RAA:RDTool>)).

=end
