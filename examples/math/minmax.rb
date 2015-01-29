#!/usr/bin/env ruby
require("gsl")

p GSL::MAX(1, 2)
p GSL::MAX(1, 2.0)
p GSL::MAX(3, 1)
p GSL::MAX(3.0, 1)

p GSL::MIN(1, 2)
p GSL::MIN(1.0, 2.0)
p GSL::MIN(3, 1.0)
p GSL::MIN(3.0, 1)

p GSL::MAX_INT(5, 6)
p GSL::MAX_INT(5, 6.0)
p GSL::MAX_INT(6, 5)
p GSL::MAX_INT(6.0, 5)

p GSL::MIN_INT(2, 2)
p GSL::MIN_INT(2.0, 2.0)
p GSL::MIN_INT(7, 2.0)
p GSL::MIN_INT(7.0, 2)
