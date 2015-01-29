#!/usr/bin/env ruby
require("gsl")

NDATA = 1000
NBINS = 100
MAX = 8

rng = GSL::Rng.alloc()
data = GSL::Ran.gaussian(rng, 1.0, NDATA) + 2
h = GSL::Histogram.alloc(NBINS, [-MAX, MAX])
h.fill(data)

hpdf = GSL::Histogram::Pdf.alloc(h)

rng2 = GSL::Rng.alloc()
h2 = GSL::Histogram.alloc(NBINS, [-MAX, MAX])
NDATA2 = 10000
for i in 0...NDATA2 do
  val = hpdf.sample(rng2.uniform())
  h2.fill(val)
end

GSL::graph(h, h2)

sum = h.sum()      # NDATA
sum2 = h2.sum()    # NDATA2
GSL::graph(h, h2.scale(sum/sum2))
