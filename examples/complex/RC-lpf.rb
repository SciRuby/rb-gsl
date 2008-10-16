#!/usr/bin/env ruby
# Transfer function of a RC low-pass filter
require("gsl")
include GSL::CONST::NUM

class RC_LPF
  def initialize(r, c)
    @r = r
    @c = c
    @fp = 1.0/(2*Math::PI*r*c)
    @omegap = 1.0/(r*c)
    puts("Pole at #{fp} Hz")
  end

  def Vout(f)
    if f.class == GSL::Vector
      out = GSL::Vector::Complex[f.size]
      i = 0
      f.each do |freq|
        a = 1.0/GSL::Complex[1.0, freq/@fp]
        out[i] = a
        i += 1
      end
      return out
    else
      1.0/GSL::Complex(1.0, f/@fp)
    end
  end

  attr_reader :r, :c, :fp
end

# Create RC filter
R = 1*KILO    # 1 [kOhm]
C = 1*MICRO   # 1 [muF]
lpf = RC_LPF.new(R, C)

# Frequency, 1Hz - 10kHz, 100 divisions
f = GSL::Vector.logspace2(1, 10*KILO, 100)

# Transfer function
tf = lpf.Vout(f)

GSL::graph([f, tf.amp.dB], "-C -l x -g 3 -y -100 10 -X 'Frequency [Hz]' -Y 'dB' --toggle-rotate-y-label")
GSL::graph([f, tf.phase/Math::PI], "-C -l x -g 3 -X 'Frequency [Hz]' -Y 'Phase [Pi]' --toggle-rotate-y-label")


