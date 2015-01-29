#!/usr/bin/env ruby
# Turn on warnings
$-w = true

require 'gnuplot'
require 'gsl'
require 'gsl/gnuplot';

# Plot using gnuplot
Gnuplot.open do |gp|
  Gnuplot::Plot.new( gp ) do |plot|
  
    plot.xrange "[0:10]"
    plot.yrange "[-1.5:1.5]"
    plot.title  "Sin Wave Example"
    plot.xlabel "x"
    plot.ylabel "sin(x)"
    plot.pointsize 3
    plot.grid 

    x = GSL::Vector[0..10]
    y = GSL::Sf::sin(x)

    plot.data = [
      Gnuplot::DataSet.new( "sin(x)" ) { |ds|
        ds.with = "lines"
        ds.title = "String function"
        ds.linewidth = 4
      },
      
      Gnuplot::DataSet.new( [x, y] ) { |ds|
        ds.with = "linespoints"
        ds.title = "Array data"
      }
    ]

  end
end
