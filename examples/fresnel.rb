#!/usr/bin/env ruby
require("gsl")

# Abramovitz & Stegun
fc01 = 0.0999975
fs01 = 0.0005236

printf("fresnel_c(0.1):\n")
printf("Expect: %2.7f \t Calculated: %2.7f \n ", fc01, GSL::fresnel_c(0.1))
printf("fresnel_s(0.1):\n")
printf("Expect: %2.7f \t Calculated: %2.7f \n ", fs01, GSL::fresnel_s(0.1))
