#!/usr/bin/env ruby
# Frequency response of a Fabry-Perot Michelson interferometer
#
require("gsl")
CLIGHT = GSL::CONST::MKSA::SPEED_OF_LIGHT

class Michelson
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
      #      r[i] = 2.0*@laserfreq/omega*Math::sin(a)*GSL::Complex[0, -a].exp
      r[i] = 2.0*@laserfreq/omega*Math::sin(a)*GSL::Complex::exp(0, -a)
    end
    r
  end
  attr_reader :len, :laserfreq
end

class FabryPerot
  # rF: Reflectance of the front mirror
  # rE: Reflectance of the end mirror
  # tF: Transmittance of the front mirror
  def initialize(len, wavelength, rF, rE)
    @len = len.to_f
    @laserfreq = 2.0*Math::PI*CLIGHT/wavelength.to_f
    @tF = Math::sqrt(1.0 - rF*rF)
    @rF = rF
    @rE = rE
  end
  def response(f)
    if f.class != GSL::Vector; f = GSL::Vector[f]; end
    r = GSL::Vector::Complex[f.size]
    f.each_index do |i|
      freq = f[i]
      omega = 2.0*Math::PI*freq
      a = 2.0*@len*omega/CLIGHT
      r[i] = @tF*@tF*@rE/(1.0 - @rF*@rE)/(1.0 - @rF*@rE*GSL::Complex::exp(0, -a))
    end
    r
  end
  attr_reader :len, :laserfreq, :tF, :rF, :rE
end

class FPMI
  def initialize(len, wavelength, rF, rE)
    @mi = Michelson.new(len, wavelength)
    @fp = FabryPerot.new(len, wavelength, rF, rE)
  end
  def response(f)
    @mi.response(f)*@fp.response(f)
  end
  attr_reader :mi, :fp
end

# Fabry-Perot Michelson IFO of the same length
fpmi = FPMI.new(3000, 1e-6, 0.85, 1.0)

f = GSL::Vector.linspace(10, 100000, 1000)
rmi = fpmi.mi.response(f).amp/1e11
rfpmi = fpmi.response(f).amp/1e11
GSL::graph(f, rmi, rfpmi, "-C -g 3 -l x -l y -y 1e-3 1e2 -L 'Red: MI 3km, Green: FPMI 3km'")
