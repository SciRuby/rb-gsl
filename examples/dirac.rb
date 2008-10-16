#!/usr/bin/env ruby
require("gsl")
include GSL::Dirac

I = GSL::Complex[0, 1]

p Pauli1
p Pauli2
p Pauli3

p (Pauli1*Pauli2 - Pauli2*Pauli1)/I/2

p Dirac.anticommute(Pauli1, Pauli1)
p Dirac.anticommute(Pauli1, Pauli2)

#p Eye2
#p Eye4

p Beta

p Alpha1
p Alpha2
p Alpha3

p Alpha1*Alpha1
p Alpha2*Alpha2
p Alpha3*Alpha3

p GSL::Dirac.anticommute(Alpha1, Alpha1)
p GSL::Dirac.anticommute(Alpha1, Alpha2)

p GSL::Dirac.anticommute(Gamma0, Gamma0)
p GSL::Dirac.anticommute(Gamma1, Gamma1)
p GSL::Dirac.anticommute(Gamma2, Gamma2)
p GSL::Dirac.anticommute(Gamma3, Gamma3)
p GSL::Dirac.anticommute(Gamma3, Gamma0)

p IEye4*Gamma0*Gamma1*Gamma2*Gamma3

p Pauli2
p Pauli2.conjugate
p Pauli2
p Pauli3
p Pauli3.dagger

p Lambda1
p Lambda2
p Lambda3
p Lambda4
p Lambda5
p Lambda6
p Lambda7
p Lambda8



