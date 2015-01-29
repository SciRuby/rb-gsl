#!/usr/bin/env ruby
# Calculate the redshift-distance relation
# for Friedmann universe
require("./cosmology")
include Cosmology

# Matter = 1, Lambda = 0, Radiation = 0
u10 = Universe.new(1.0, 0.0, 0.0)

# Matter = 0.27, Lambda = 0.73, Radiation = 0
u37 = Universe.new(0.27, 0.73, 0.0)

# Matter = 0.0, Lambda = 1.0, Radiation = 0
u01 = Universe.new(0.0, 1.0, 0.0)

begin
  file = File.open("friedmann.dat", "w")
  z = 0.01
  while z <= 3
    # comoving distance, luminosity distance
    chi10, ldist10 = u10.luminosity_distance(z)
    chi37, ldist37 = u37.luminosity_distance(z)
    chi01, ldist01 = u01.luminosity_distance(z)
    file.printf("%e %e %e %e %e %e %e\n",
                z, chi10, ldist10, chi37, ldist37, chi01, ldist01)
    z += 0.1
  end
ensure
  file.close
end

#puts("Data file \"friedmann.dat\" created.")
#puts("Try \"gnuplot -persist friedmann.gp\"")
system("gnuplot -persist friedmann.gp")
File.delete("friedmann.dat")
