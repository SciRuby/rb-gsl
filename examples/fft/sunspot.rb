#!/usr/bin/env ruby
# Analysis of the solar activity of 11-years cycle
# from the number of sunspots.
# This example is taken from the MATLAB user's manual Chap 13.

require("gsl")

year, sunspot = GSL::Vector.filescan("sunspot.dat")
N = year.size

ffted = sunspot.fft

power = GSL.sqrt(ffted[1..(N-2)].to_complex2.abs2)*2/N
freq = GSL::Vector.linspace(1, N/2, power.size)/N
period = 1.0/freq
GSL::graph(period, power, "-C -g 3 -x 0 40 -X 'Period [year]'")
