=begin
= NArray compatibilities
=== Contents:
(1) ((<Data type conversions|URL:narray.html#1>))
(2) ((<Methods which accepts NArray|URL:narray.html#2>))

== Data type conversions
=== GSL to NArray

--- GSL::Vector#to_na
--- GSL::Vector#to_nvector
--- GSL::Matrix#to_na
--- GSL::Matrix#to_nmatrix
    Convert GSL objects to NArray. The data contained by the GSL objects
    are copied to a newly allocated memory block of the NArray objects created.


--- GSL::Vector#to_na_ref
--- GSL::Vector#to_nvector_ref
--- GSL::Matrix#to_na_ref
--- GSL::Matrix#to_nmatrix_ref
    Create NArray-ref objects from GSL data. The memory block of the GSL
    objects are shared with the NArray-ref objects.

    Example:
      irb(main):009:0> v = Vector::Int[0..5]
      => GSL::Vector::Int
      [ 0 1 2 3 4 5 ]
      irb(main):010:0> na = v.to_nvector_ref
      => NVector(ref).int(6): 
      [ 0, 1, 2, 3, 4, 5 ]
      irb(main):011:0> na[3] = 999
      => 999
      irb(main):012:0> v
      => GSL::Vector::Int
      [ 0 1 2 999 4 5 ]

=== NArray to GSL
--- NArray#to_gv
--- NArray#to_gm
    Convert NArray objects to (({GSL::Vector})) or (({GSL::Matrix})). 
    The data contained by the NArray objects
    are copied to a newly allocated memory block of the GSL objects created.

--- NArray#to_gv_view
--- NArray#to_gm_view
    Create (({GSL::Vector::View})) or (({GSL::Matrix::View})) objects from NArray. 
    The memory block of the NArray objects are shared with the View objects.

    Example:
      irb(main):024:0> na = NArray[0, 1, 2, 3, 4, 5]
      => NArray.int(6): 
      [ 0, 1, 2, 3, 4, 5 ]
      irb(main):025:0> b = na.to_gv_int_view
      => GSL::Vector::Int::View
      [ 0 1 2 3 4 5 ]
      irb(main):026:0> b[2] = -99
      => -99
      irb(main):027:0> na
      => NArray.int(6): 
      [ 0, 1, -99, 3, 4, 5 ]

== Methods which accepts NArray
=== (({GSL})) module
--- GSL::graph()
--- GSL::log1p(x)
--- GSL::expm1(x)
--- GSL::hypot(x, y)
--- GSL::acosh(x)
--- GSL::asinh(x)
--- GSL::atanh(x)
--- GSL::pow(x, a)
--- GSL::pow_int(x, n)
--- GSL::pow_2(x), ..., GSL::pow_9(x)

=== (({GSL::Sf})) module
((<Any|URL:sf.html>))

=== (({GSL::Linalg})) module
--- GSL::Linalg::LU.decomp(na)
--- GSL::Linalg::LU.solve(lu, b)
--- GSL::Linalg::LU.svx(lu, bx)
--- GSL::Linalg::LU.det(lu, sign)
--- GSL::Linalg::LU.lndet(lu)
--- GSL::Linalg::LU.invert(lu, perm)
--- GSL::Linalg::QR.decomp(m)
--- GSL::Linalg::QR.solve(qr, tau, b)	
--- GSL::Linalg::QR.svx(qr, tau, bx)
--- GSL::Linalg::SV.decomp(m)
--- GSL::Linalg::SV.solve(u, v, s, b)
--- GSL::Linalg::SV.svx(u, v, s, bx)
--- GSL::Linalg::Cholesky.decomp(m)
--- GSL::Linalg::Cholesky.solve(u, v, s, b)
--- GSL::Linalg::Cholesky.svx(u, v, s, bx)
--- GSL::Linalg::HH.solve(m, b)
--- GSL::Linalg::HH.svx(m, bx)

=== (({GSL::Eigen})) module
--- GSL::Eigen::symm(na)
--- GSL::Eigen::symmv(na)

=== (({GSL::FFT})) module
((<Many|URL:FFT.html>))

=== (({GSL::Function})) class
--- GSL::Function#eval
--- GSL::Function#deriv_central(x, h)
--- GSL::Function#deriv_forward(x, h)
--- GSL::Function#deriv_backward(x, h)
--- GSL::Function#diff_central(x, h)
--- GSL::Function#diff_forward(x, h)
--- GSL::Function#diff_backward(x, h)

=== (({GSL::Ran})) and (({GSL::Cdf})) module
((<Many|URL:randist.html>))

=== (({GSL::Stats})) module
((<Any|URL:stats.html>))

=== (({GSL::Interp})) and (({GSL::Spline})) class
--- GSL::Interp#init
--- GSL::Interp#eval
--- GSL::Spline#init
--- GSL::Spline#eval

=== (({GSL::Deriv})) and (({GSL::Diff})) module
--- GSL::Deriv.central(f, x, h)
--- GSL::Deriv.forward(f, x, h)
--- GSL::Deriv.backward(f, x, h)
--- GSL::Diff.central(f, x, h)
--- GSL::Diff.forward(f, x, h)
--- GSL::Diff.backward(f, x, h)

=== (({GSL::Cheb})) class
--- GSL::Cheb#eval(x)
--- GSL::Cheb#eval_n(n, x)

=== (({GSL::Wavelet})) class
((<Many|URL:wavelet.html>))

((<prev|URL:tensor.html>))

((<Reference index|URL:ref.html>))
((<top|URL:index.html>))

=end
