#!/usr/bin/env ruby
require("gsl")

N = 2048
SAMPLING = 1000   # 1 kHz
TMAX = 1.0/SAMPLING*N
FREQ1 = 50
FREQ2 = 120
t = GSL::Vector.linspace(0, TMAX, N)
x = GSL::Sf::sin(2*Math::PI*FREQ1*t) + GSL::Sf::sin(2*Math::PI*FREQ2*t)
y = x.fft

y2 = y.subvector(1, N-2).to_complex2
mag = y2.abs
phase = y2.arg
f = GSL::Vector.linspace(0, SAMPLING/2, mag.size)
GSL::graph(f, mag, "-C -g 3 -x 0 200 -X 'Frequency [Hz]'")
