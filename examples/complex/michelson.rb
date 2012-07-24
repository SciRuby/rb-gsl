#!/usr/bin/env ruby
# Frequency response of a Michelson interferometer
#
require("gsl")

class Michelson
  CLIGHT = GSL::CONST::MKSA::SPEED_OF_LIGHT

  def initialize(len, wavelength)
    @len = len.to_f
    @laserfreq = 2.0*Math::PI*CLIGHT/wavelength.to_f
  end

  def response(f)
    if f.class != GSL::Vector; f = GSL::Vector[f]; end
    r = GSL::Vector::Complex[f.size]
    f.each_index do |i|
      freq = f[i]
      omega = 2.0*Math::PI*freq
      a = @len*omega/CLIGHT
      r[i] = 2.0*@laserfreq/omega*Math::sin(a)*GSL::Complex::exp([0, -a])
    end
    r
  end
  attr_reader :len, :laserfreq
end

# 3km
mi1 = Michelson.new(3000, 1e-6)
# 75km
mi2 = Michelson.new(75000, 1e-6)

f = GSL::Vector.linspace(10, 100000, 1000)
r1 = mi1.response(f).amp
r2 = mi2.response(f).amp
GSL::graph(f, r2/1e11, r1/1e11, "-C -g 3 -l x -l y -y 1e-2 1e2 -L 'Transfer function of Michelson interferometer'")
