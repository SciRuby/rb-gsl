#!/usr/bin/env ruby
require("gsl")
include Math

f = GSL::Function::alloc{ |x|
  exp(-x)/sqrt(x)
}

xmin = 0.0
xmax = 1.0
limit = 1000

=begin
puts("QNG")
p f.integration_qng(xmin, xmax, 0.0, 1.0e-7)
p f.integration_qng(xmin, xmax)
p f.integration_qng([xmin, xmax])
p f.integration_qng([xmin, xmax], [0.0, 1.0e-7])
p f.integration_qng([xmin, xmax], 0.0, 1.0e-7)
p f.integration_qng(xmin, xmax, [0.0, 1.0e-7])
=end

puts("QAG")
p f.integration_qag(xmin, xmax, 0.0, 1.0e-7, limit, GSL::Integration::GAUSS15)
p f.integration_qag([xmin, xmax], 0.0, 1.0e-7, limit, GSL::Integration::GAUSS15)
p f.integration_qag([xmin, xmax], [0.0, 1.0e-7], limit, GSL::Integration::GAUSS15)
p f.integration_qag(xmin, xmax, [0.0, 1.0e-7], limit, GSL::Integration::GAUSS15)
p f.integration_qag(xmin, xmax, GSL::Integration::GAUSS15)
p f.integration_qag([xmin, xmax], GSL::Integration::GAUSS15)
w = GSL::Integration::Workspace.alloc(2000)
p f.integration_qag(xmin, xmax, [0.0, 1.0e-7], limit, GSL::Integration::GAUSS15, w)
p f.integration_qag(xmin, xmax, 0.0, 1.0e-7, GSL::Integration::GAUSS15, w)
p f.integration_qag([xmin, xmax], 0.0, 1.0e-7, GSL::Integration::GAUSS15, w)
p f.integration_qag([xmin, xmax], [0.0, 1.0e-7], GSL::Integration::GAUSS15, w)
p f.integration_qag([xmin, xmax], [0.0, 1.0e-7], limit, GSL::Integration::GAUSS15, w)

puts("QAGS")
p f.integration_qags(xmin, xmax)
p f.integration_qags([xmin, xmax])
p f.integration_qags(xmin, xmax, 0.0, 1e-7)
p f.integration_qags(xmin, xmax, [0.0, 1e-7])
p f.integration_qags([xmin, xmax], [0.0, 1e-7])
p f.integration_qags([xmin, xmax], 0.0, 1e-7)
p f.integration_qags([xmin, xmax], 0.0, 1e-7, limit)
p f.integration_qags([xmin, xmax], [0.0, 1e-7], limit)
p f.integration_qags(xmin, xmax, [0.0, 1e-7], limit)
p f.integration_qags([xmin, xmax], 0.0, 1e-7, limit, w)
p f.integration_qags(xmin, xmax, 0.0, 1e-7, limit, w)
p f.integration_qags(xmin, xmax, [0.0, 1e-7], w)
p f.integration_qags([xmin, xmax], [0.0, 1e-7], w)

p f.integration_qags([xmin, xmax], limit)
p f.integration_qags(xmin, xmax, limit)
p f.integration_qags([xmin, xmax], w)
p f.integration_qags(xmin, xmax, w)
p f.integration_qags(xmin, xmax, limit, w)
p f.integration_qags([xmin, xmax], limit, w)

puts("QAGP")
p f.integration_qagp([xmin, xmax])
p f.integration_qagp([xmin, xmax], [0.0, 1e-7])
p f.integration_qagp([xmin, xmax], 0.0, 1e-7)
p f.integration_qagp([xmin, xmax], 0.0, 1e-7, limit)
p f.integration_qagp([xmin, xmax], [0.0, 1e-7], limit)
p f.integration_qagp([xmin, xmax], 0.0, 1e-7, limit, w)
p f.integration_qagp([xmin, xmax], [0.0, 1e-7], w)

p f.integration_qagp([xmin, xmax], limit)
p f.integration_qagp([xmin, xmax], w)
p f.integration_qagp([xmin, xmax], limit, w)

