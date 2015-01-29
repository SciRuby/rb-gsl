#!/usr/bin/env ruby
require("gsl")
include GSL

puts("GSL_POSINF")
p GSL::POSINF
p GSL::POSINF.class
p isinf(GSL::POSINF)
p isinf?(GSL::POSINF)
p isnan(GSL::POSINF)
p isnan?(GSL::POSINF)
p finite(GSL::POSINF)
p finite?(GSL::POSINF)

puts("GSL_NEGINF")
p GSL::NEGINF
p GSL::NEGINF.class
p isinf(GSL::NEGINF)
p isinf?(GSL::NEGINF)
p isnan(GSL::NEGINF)
p isnan?(GSL::NEGINF)
p finite(GSL::NEGINF)
p finite?(GSL::NEGINF)

puts("GSL_NAN")
p GSL::NAN
p GSL::NAN.class
p isinf(GSL::NAN)
p isinf?(GSL::NAN)
p isnan(GSL::NAN)
p isnan?(GSL::NAN)
p finite(GSL::NAN)
p finite?(GSL::NAN
)
