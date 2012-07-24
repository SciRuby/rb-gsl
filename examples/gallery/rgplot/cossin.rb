#!/usr/bin/env ruby
#begin
#  require 'rubygems'
#  require_gem "gnuplot"         # Try using rubygem
#ensure
  require 'gnuplot'             # No gem, use traditional require
#end

require "gsl"

# Add the to_gplot method to Vector since its not already built in.  This
# might be worthwhile adding to the core GSL stuff.

x = GSL::Vector.linspace(0, 2*M_PI, 100)
s = Sf::sin(x)
c = Sf::cos(x)

# Now generate the actual plot
Gnuplot::open do |gp|
  Gnuplot::Plot.new( gp ) do |plot|

    plot.title "GSL plotting example"
    plot.data = [
      Gnuplot::DataSet.new( [x, c] ) { |ds|
        ds.title = "cos(x)"
        ds.with = "lines"
      },
      Gnuplot::DataSet.new( [x, s] ) { |ds|
        ds.title = "sin(x)"
        ds.with = "lines"
      }
    ]
  end
end

