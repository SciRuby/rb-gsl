
require("gsl")
include GSL
include GSL::CONST
include GSL::CONST::CGSM
include Math

module Cosmology

  # Inverse of H(z), H(z): Hubble parameter determined by Friedmann equation
  InvHz = Function.alloc { |z, unv|
    matter = unv.matter
    radiation = unv.radiation
    lambda = unv.lambda         # Cosmological constant
    k = unv.k                   # Curvature
    zplus1 = z + 1.0
    z12 = zplus1*zplus1
    1.0/sqrt(k*z12 + lambda + matter*z12*zplus1 + radiation*z12*z12)
  }
  WSpace = GSL::Integration::Workspace.alloc(1000)

  class Universe
    H2H0 = 3.24e-18
    G = GRAVITATIONAL_CONSTANT
    RAD_DENSITY = 8.0*pow_5(PI)*pow_4(BOLTZMANN)/15.0/pow_3(PLANCKS_CONSTANT_H*SPEED_OF_LIGHT)
    #  p SPEED_OF_LIGHT/(0.71*100*1e5*1e3)*13

    def initialize(matter = 0.3, lambda = 0.7, cmbT = 2.7, h = 0.7)
      # Hubble parameter at present
      @h = h
      @H0 = H2H0*h

      # Density parameters
      @matter = matter
      @lambda = lambda
      rad = RAD_DENSITY*pow_4(cmbT)
      @radiation = rad*8*PI*G/3/@H0/@H0/SPEED_OF_LIGHT/SPEED_OF_LIGHT
      @Omega = @matter + @lambda + @radiation
      @k = 1.0 - @Omega
      @q = @matter/2.0 + @radiation - @lambda  # Decceleration parameter
    end

    def comoving_distance(z)
      InvHz.set_params(self)
      InvHz.integration_qags([0.0, z], WSpace)[0]
    end

    # r: comoving distance
    def conformal_factor(r)
      tmp = sqrt(@k)
      if @Omega > 1.0
        sin(tmp*r)/tmp
      elsif @Omega < 1.0
        sinh(tmp*r)/tmp
      else
        r
      end
    end

    def luminosity_distance(z)
      r = self.comoving_distance(z)
      s = self.conformal_factor(r)
      [r, (1.0 + z)*s]  # comoving distance, luminosity distance
    end

    attr_reader :h
    attr_reader :H0
    attr_reader :matter
    attr_reader :lambda
    attr_reader :radiation
    attr_reader :Omega
    attr_reader :k
    attr_reader :q
  end
end
