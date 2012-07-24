#!/usr/bin/env ruby
require("gsl")
require("nimage")

def mandel(w, h)
  zoom = 3.5
  z = (GSL::Matrix::Complex[w, 1].indgen!/w - 0.65)*zoom + (GSL::Matrix::Complex[1, h].indgen!/h - 0.5)*zoom*GSL::Complex[0, 1]
  c = z.clone
  a = GSL::Vector::Int[w, h]
  idx = GSL::Vector::Int[h, w].indgen!

  for i in 1..30
    z = z*z + c
    idx_t,idx_f = (z.abs>2).where2
    a[idx[idx_t]] = i
    break if idx_f.size==0
    idx = idx[idx_f]
    z = z[idx_f]
    c = c[idx_f]
  end
  a
end

NImage.show mandel(400,400).to_na_ref

print "Hit return key..."
STDIN.getc
