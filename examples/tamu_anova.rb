#!/usr/bin/env ruby
require("gsl")

data = GSL::Vector[88.60,73.20,91.40,68.00,75.20,63.00,53.90,69.20,
50.10,71.50,44.90,59.50,40.20,56.30,38.70,31.00,
39.60,45.30,25.20,22.70]
factor = GSL::Vector::Int[1,1,1,1,1,2,2,2,2,2,3,3,3,3,3,4,4,4,4,4]

table = GSL::TAMU_ANOVA::Table.oneway(data, factor, 4)
table.print

data = GSL::Vector[45.50,45.30,45.40,44.40,44.60,43.90,44.60,44.00,44.20,
43.90,44.70,44.20,44.00,43.80,44.60,43.10,46.00,45.90,
44.80,46.20,45.10,45.50]
factor = GSL::Vector::Int[1,1,1,1,1,1,1,1,2,2,2,2,2,2,2,2,3,3,3,3,3,3]

table = GSL::TAMU_ANOVA::Table.oneway(data, factor, 3)
table.print
