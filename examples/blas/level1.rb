#!/usr/bin/env ruby
require("gsl")
include GSL
include GSL::Blas

x = GSL::Vector[1, 2, 3]
y = GSL::Vector[4, 5, 6]

p x.blas_ddot(y)
p GSL::Blas.ddot(x, y)

xz = GSL::Vector::Complex.alloc([1, 0], [2, 0], [3, 0])
yz = GSL::Vector::Complex.alloc([4, 0], [5, 0], [6, 0])
p zdotu(xz, yz)

p x.nrm2
p GSL::Blas.dnrm2(x)
p xz.nrm2
p GSL::Blas.dznrm2(xz)

p x.asum
p GSL::Blas.dasum(x)
p xz.asum
p GSL::Blas.dzasum(xz)

x.swap(y)
p x
p y

GSL::Blas.dswap(x, y)
p x
p y

xz.swap(yz)
p xz
p yz

GSL::Blas.zswap(xz, yz)
p xz
p yz

p x.axpy(2, y)
p y
p GSL::Blas.daxpy(2, x, y)
p y

p x.axpy!(2, y)
p y
p GSL::Blas.daxpy!(2, x, y)
p y

y = GSL::Vector[4, 5, 6]

az = GSL::Complex.alloc(2, 0)
p xz.zaxpy(az, yz)
p yz
p GSL::Blas.zaxpy(az, xz, yz)
p yz

az = GSL::Complex[2, 0]
p xz.axpy!(az, yz)
p yz
p GSL::Blas.zaxpy!(az, xz, yz)
p yz

yz = GSL::Vector::Complex.alloc([4, 0], [5, 0], [6, 0])

p x.scal(2)
p x

p x.scal!(2)
p x

x = GSL::Vector[1, 2, 3]

p xz.zscal(az)
p xz.zdscal(2)

p y
p GSL::Blas.drot(x, y, 2, 3)
p drot(x, y, 2, 3)
