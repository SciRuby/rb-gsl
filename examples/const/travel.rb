#!/usr/bin/env ruby
require("gsl")

include GSL::CONST::MKSA
puts("In MKSA unit")

c  = SPEED_OF_LIGHT;
au = ASTRONOMICAL_UNIT;
minutes = MINUTE;

# distance stored in meters
r_earth = 1.00 * au;
r_mars  = 1.52 * au;

t_min = (r_mars - r_earth) / c;
t_max = (r_mars + r_earth) / c;

printf("light travel time from Earth to Mars:\n");
printf("c = %e [m/s]\n", c)
printf("AU = %e [m]\n", au)
printf("minutes = %e [s]\n", minutes)
printf("minimum = %.1f minutes\n", t_min / minutes);
printf("maximum = %.1f minutes\n\n", t_max / minutes);


include GSL::CONST::CGSM
puts("In CGSM unit")

c  = SPEED_OF_LIGHT;
au = ASTRONOMICAL_UNIT;
minutes = MINUTE;

# distance stored in meters
r_earth = 1.00 * au;
r_mars  = 1.52 * au;

t_min = (r_mars - r_earth) / c;
t_max = (r_mars + r_earth) / c;

printf("light travel time from Earth to Mars:\n");
printf("c = %e [cm/s]\n", c)
printf("AU = %e [cm]\n", au)
printf("minutes = %e [s]\n", minutes)
printf("minimum = %.1f minutes\n", t_min / minutes);
printf("maximum = %.1f minutes\n", t_max / minutes);
