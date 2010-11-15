=begin

= Ruby/GSL Reference
(See also ((<Gnu Scientific Library -- Reference Manual|URL:http://www.gnu.org/software/gsl/manual/html_node/>)))

== Front Matter
This document describes the modules, classes and the methods of Ruby/GSL.  This
includes cut-and-paste from the
((<GNU Scientific Library -- Reference Manual|URL:http://www.gnu.org/software/gsl/manual/html_node/>)),
and documents of the extention libraries.

=== Copyright of the GSL Reference

Copyright 1996, 1997, 1998, 1999, 2000, 2001, 2002, 2003, 2004, 2005, 2006,
2007, 2008 The GSL Team.

Permission is granted to copy, distribute and/or modify this document under the
terms of the GNU Free Documentation License, Version 1.3 or any later version
published by the Free Software Foundation; with the Invariant Sections being
"GNU General Public License" and "Free Software Needs Free Documentation", the
Front-Cover text being "A GNU Manual", and with the Back-Cover Text being (a)
(see below). A copy of the license is included in the section entitled "GNU
Free Documentation License".  (a) The Back-Cover Text is: "You have the freedom
to copy and modify this GNU Manual."

=== Copyright of this reference

2009,2010 Yoshiki Tsunesada, David MacMahon

Permission is granted to copy, distribute and/or modify this document under
the terms of the GNU Free Documentation License.

== Ruby/GSL Reference

  (1) ((<Introduction|URL:intro.html>))
  (2) ((<Using RubyGSL|URL:use.html>))
  (3) ((<Error Handling|URL:ehandling.html>))
  (4) ((<Mathematical Functions|URL:math.html>))
  (5) ((<Complex Numbers|URL:complex.html>))
  (6) ((<Polynomials|URL:poly.html>))
  (7) ((<Special Functions|URL:sf.html>))
  (8) ((<Vectors|URL:vector.html>)) and ((<Matrices|URL:matrix.html>))
  (9) ((<Permutations|URL:perm.html>))
  (10) ((<Combinations|URL:combi.html>))
  (11) ((<Multiset|URL:sort.html>)) (GSL-1.14)
  (12) ((<Sorting|URL:sort.html>))
  (13) ((<BLAS Support|URL:blas.html>))
  (14) ((<Linear Algebra|URL:linalg.html>))
  (15) ((<Eigen Systems|URL:eigen.html>))
  (16) ((<Fast Fourier Transform|URL:fft.html>))
  (17) ((<Numerical Integration|URL:integration.html>))
  (18) ((<Random Number Generation|URL:rng.html>))
  (19) ((<Quasi-Random Sequences|URL:qrng.html>))
  (20) ((<Random Number Distributions|URL:randist.html>))
  (21) ((<Statistics|URL:stats.html>))
  (22) ((<1d-Histograms|URL:hist.html>)), ((<2d-Histograms|URL:hist2d.html>)) and ((<3d-Histograms|URL:hist3d.html>))
  (23) ((<N-tuples|URL:ntuple.html>))
  (24) ((<Monte-Carlo Integration|URL:monte.html>))
  (25) ((<Simulated Annealing|URL:siman.html>))
  (26) ((<Ordinary Differential Equations|URL:odeiv.html>))
  (27) ((<Interpolation|URL:interp.html>))
  (28) ((<Numerical Differentiation|URL:diff.html>))
  (29) ((<Chebyshev Approximations|URL:cheb.html>))
  (30) ((<Series Acceleration|URL:sum.html>))
  (31) ((<Wavelet Transforms|URL:wavelet.html>)) (GSL-1.6 feature)
  (32) ((<Discrete Hankel Transforms|URL:dht.html>))
  (33) ((<One dimensional Root-Finding|URL:roots.html>))
  (34) ((<One dimensional Minimization|URL:min.html>))
  (35) ((<Multidimensional Root-Finding|URL:multiroot.html>))
  (36) ((<Multidimensional Minimization|URL:multimin.html>))
  (37) ((<Least-Squares Fitting|URL:fit.html>))
  (38) ((<Nonlinear Least-Squares Fitting|URL:nonlinearfit.html>))
  (39) ((<Basis Splines|URL:bspline.html>))  
  (40) ((<Physical Constants|URL:const.html>))
  (41) ((<Graphics|URL:graph.html>))
  (42) Supported GSL add-on packages
       (1) ((<rngextra|URL:rngextra.html>))
       (2) ((<Tensor manipulations|URL:tensor.html>))
       (3) OOL: Open Optimization library (see examples/ool/*.rb)
       (4) CQP and Bundle (see examples/multimin/cqp.rb, bundle.rb)
       (5) quartic
       (6) jacobi (see examples/jacobi/*.rb)
       (7) ((<NDLINEAR: multi-linear, multi-parameter least squares fitting|URL:ndlinear.html>))
       (8) ((<ALF: associated Legendre polynomials|URL:alf.html>))
  (43) ((<NArray compatibilities|URL:narray.html>))
  (44) ((<Changes since RubyGSL 1.10.3|URL:changes.html>))
  (44) ((<GNU Free Documentation Licence|URL:http://www.gnu.org/software/gsl/manual/html_node/GNU-Free-Documentation-License.html>))

((<next|URL:intro.html>))

((<top|URL:index.html>))
=end
