=begin
= Linear Algebra

Contents:
(1) ((<LU Decomposition|URL:linalg.html#1>))
(2) ((<QR Decomposition|URL:linalg.html#2>))
(3) ((<QR Decomposition with Column Pivoting|URL:linalg.html#4>))
(4) ((<Singular Value Decomposition|URL:linalg.html#5>))
(5) ((<Cholesky Decomposition|URL:linalg.html#6>))
(6) ((<Tridiagonal Decomposition of Real Symmetric Matrices|URL:linalg.html#6>))
(7) ((<Tridiagonal Decomposition of Hermitian Matrices|URL:linalg.html#7>))
(8) ((<Hessenberg Decomposition of Real Matrices|URL:linalg.html#8>))
(9) ((<Hessenberg-Triangular Decomposition of Real Matrices|URL:linalg.html#9>))
(10) ((<Bidiagonalization|URL:linalg.html#10>))
(11) ((<Householder Transformations|URL:linalg.html#11>))
(12) ((<Householder solver for linear systems|URL:linalg.html#12>))
(13) ((<Tridiagonal Systems|URL:linalg.html#13>))
(14) ((<Balancing|URL:linalg.html#14>))
(15) ((<NArray|URL:linalg.html#15>))

== LU Decomposition
--- GSL::Linalg::LU.decomp(A)
--- GSL::Matrix#LU_decomp
    These method calculate the LU decomposition of the matrix. The returned
    value is an array of (({[LU, perm, sign]})).

    Examples:

    (1) Singleton method of the (({GSL::Linalg::LU})) module

          irb(main):012:0> m = Matrix[1..9, 3, 3]
          => GSL::Matrix: 
          [ 1.000e+00 2.000e+00 3.000e+00 
            4.000e+00 5.000e+00 6.000e+00 
            7.000e+00 8.000e+00 9.000e+00 ]
          irb(main):013:0> lu, perm, sign = Linalg::LU.decomp(m)

    (2) Instance method of (({GSL::Matrix})) class

          irb(main):012:0> lu, perm, sign = m.LU_decomp

--- GSL::Linalg::LU.solve(A, b)
--- GSL::Linalg::LU.solve(lu, perm, b)
--- GSL::Matrix#LU_solve(b)
--- GSL::Linalg::LUMatrix#solve(perm, b)
    The following is an example to solve a linear system 

       A x = b, b = [1, 2, 3, 4]

    using LU decomposition.

    (1) Singleton method of the (({GSL::Linalg::LU})) module

          A = Matrix[[0.18, 0.60, 0.57, 0.96], [0.41, 0.24, 0.99, 0.58],
                     [0.14, 0.30, 0.97, 0.66], [0.51, 0.13, 0.19, 0.85]]
          lu, perm, sign = A.LU_decomp         
          b = Vector[1, 2, 3, 4]
          x = Linalg::LU.solve(lu, perm, b)    

    (2) Instance method of (({GSL::Linalg::LUMatrix})) class

          lu, perm, sign = A.LU_decomp  # lu is an instance of Linalg::LUMatrix class
          b = Vector[1, 2, 3, 4]
          x = lu.solve(perm, b)    

    (3) Solve directly

          x = Linalg::LU.solve(A, b)  # LU decomposition is calculated internally (A is not modified)
 
--- GSL::Linalg::LU.svx(A, b)
--- GSL::Linalg::LU.svx(lu, perm, b)
--- GSL::Matrix#svx(b)
--- GSL::Linalg::LUMatrix#svx(perm, b)
    These solve the system Ax = b. The input vector ((|b|)) is overwitten by
    the solution ((|x|)).

--- GSL::Linalg::LU.refine(A, lu, perm, b, x)
    This method applys an iterative improvement to ((|x|)), 
    the solution of ((|A x = b|)), using the LU decomposition of ((|A|)).

--- GSL::Linalg::LU.invert(A)
--- GSL::Linalg::LU.invert(lu, perm)
--- GSL::Matrix#invert
--- GSL::Linalg::LUMatrix#invert(perm)
    These computes and returns the inverse of the matrix.

--- GSL::Linalg::LU.det(A)
--- GSL::Linalg::LU.det(lu, signum)
--- GSL::Matrix#det
--- GSL::Linalg::LUMatrix#det(signum)
    These methods return the determinant of the matrix.

=== ((<Complex LU decomposition|URL:linalg_complex.html>))

== QR decomposition

--- GSL::Linalg::QR_decomp(A)
--- GSL::Linalg::QR.decomp(A)
--- GSL::Matrix#QR_decomp
    These compute QR decomposition of the matrix and return an array [QR, tau].

    (1) Singleton method of the module (({GSL::Linalg}))
          qr, tau = Linalg::QR_decomp(m)
          p qr.class                 # GSL::Linalg::QRMatrix, subclass of GSL::Matrix
          p tau.class                # GSL::Linalg::TauVector, subclass of GSL::Vector
    (2) Singleton method of the module (({GSL::Linalg:QR}))
          qr, tau = Linalg::QR.decomp(m)
    (3) Instance method of (({GSL::Matrix}))
          qr, tau = m.QR_decomp

--- GSL::Linalg::QR.solve(A, b)
--- GSL::Linalg::QR.solve(QR, tau, b)
--- GSL::Matrix#QR_solve(b)
--- GSL::Linalg::QRMatrix#solve(tau, b)
    Solve the system A x = b using the QR decomposition.

    * Ex1:
        m = Matrix.alloc(...)
        b = Vector.alloc(...)
        x = Linalg::QR.solve(m, b)
    * Ex2:
        x = m.QR_solve(b)
    * Ex3:
        qr, tau = Linalg::QR.decomp(m)   # or m.QR_decomp
        x = Linalg::QR.solve(qr, tau, b)
    * Ex4:
        qr, tau = m.QR_decomp
        x = qr.solve(tau, b)
 
--- GSL::Linalg::QR.svx(A, x)
--- GSL::Linalg::QR.svx(QR, tau, x)
--- GSL::Matrix#QR_svx(x)
--- GSL::Linalg::QRMatrix#svx(tau, x)
    Solve the system A x = b. The input vector ((|x|)) is first give by
    the right-hand side vector ((|b|)), and is overwritten by the solution.

--- GSL::Linalg::QR.unpack(QR, tau)
--- GSL::Linalg::QRMatrix#unpack(tau)
    Unpack the encoded QR decomposition ((|QR,tau|)) and return an array
    (({[Q, R]})).

    Ex:
      irb(main):010:0> m = Matrix[1..9, 3, 3]
      => GSL::Matrix: 
      [ 1.000e+00 2.000e+00 3.000e+00 
        4.000e+00 5.000e+00 6.000e+00 
        7.000e+00 8.000e+00 9.000e+00 ]
      irb(main):011:0> qr, tau = m.QR_decomp
      irb(main):012:0> q, r = qr.unpack(tau)  
      irb(main):013:0> q*r               # Reconstruct the metrix m
      => GSL::Matrix: 
      [ 1.000e+00 2.000e+00 3.000e+00 
        4.000e+00 5.000e+00 6.000e+00 
        7.000e+00 8.000e+00 9.000e+00 ]

--- GSL::Linalg::QR.QRsolve(Q, R, tau)
    This method solves the system (({R x = Q^T b})) for (({x})). 
    It can be used when the QR decomposition of a matrix is available 
    in unpacked form as ((|Q,R|)).

== QR Decomposition with Column Pivoting
--- GSL::Linalg::QRPT.decomp(A)
--- GSL::Matrix#QRPT_decomp
    These methods factorize the M-by-N matrix ((|A|)) into the QRP^T decomposition A = Q R P^T, and return an array (({[QR, tau, perm, signum]})).

    * Ex1:
        require("gsl")
        include GSL::Linalg
        m = Matrix.alloc(...)
        qr, tau, perm = QRPT.decomp(m)
        p qr.class                 # GSL::Linalg::QRPTMatrix, subclass of GSL::Matrix

    * Ex2:
        qr, tau, perm = m.QROT_decomp

--- GSL::Linalg::QRPT.decomp2(A)
--- GSL::Matrix#QRPT_decomp2
    These return an array (({[Q, R, tau, perm, signum]})).

    * Ex
        q, r, tau, perm = QRPT.decomp2(m)
        p q.class                  <----- GSL::Linalg::QMatrix
        p r.class                  <----- GSL::Linalg::RMatrix

--- GSL::Linalg::QRPT.solve(m, b)
--- GSL::Linalg::QRPT.solve(qr, tau, perm, b)
--- GSL::Matrix#QRPT_solve(A, b)
--- GSL::Linalg::QRPQMatrix#solve(qr, tau, perm, b)
    These methods solve the system (({A x = b})) using the QRP^T decomposition of 
    ((|A|)) into ((|QR, tau, perm|)). The solution ((|x|)) is returned as a Vector.

    * Ex1:
        m = Matrix.alloc(...)
        qr, tau, perm = m.QRPT_decomp
        b = Vector.alloc([1, 2, 3, 4])
        x = Linalg::QRPT.solve(qr, tau, perm, b)
    * Ex2:
        x = Linalg::QRPT.solve(m, b)
    * Ex3:
        x = qr.solve(tau, p, b)
    * Ex4:
        x = m.QRPT_solve(b)

--- GSL::Linalg::QRPT.svx(m, b)
--- GSL::Linalg::QRPT.svx(qr, tau, perm, b)
--- GSL::Matrix#QRPT_svx(A, b)
    These methods solve the system (({A x = b})) using the QRP^T decomposition of 
    ((|A|)) into ((|QR, tau, perm|)). The input ((|b|)) is overwritten by the solution
    ((|x|)).

--- GSL::Linalg::QRPT.QRsolve(q, r, tau, perm, b)
    This method solves the system (({R P^T x = Q^T b})) for x. 
    It can be used when the QR decomposition of a matrix is available in 
    unpacked form as ((|q, r|)) obtained by the method (({decomp2})).

    * Ex:
       q, r, tau, perm = QRPT_decomp2
       x = Linalg::QRPT.QRsolve(q, r, perm, b)

--- GSL::Linalg::QRPT.update(q, r, perm, u, v)
--- GSL::Linalg::QRPT.Rsolve(qr, perm, b)
--- GSL::Linalg::QRPTMatrix#Rsolve(perm, b)
--- GSL::Linalg::QRPT.Rsvx(qr, perm, b)
--- GSL::Linalg::QRPTMatrix#Rsvx(perm, b)

== Singular Value Decomposition
--- GSL::Linalg::SV.decomp(A[, work])
--- GSL::Matrix#SV_decomp([work])
    These methods factorize the M-by-N matrix ((|A|)) into the singular value 
    decomposition (({A = U S V^T})) using  the Golub-Reinsch SVD algorithm,
    and return an array (({[U, V, S]})).

    Ex: 
       irb(main):020:0* m = Matrix[[3, 5, 2], [5, 1, 4], [7, 6, 3]]
       => GSL::Matrix: 
       [ 3.000e+00 5.000e+00 2.000e+00 
         5.000e+00 1.000e+00 4.000e+00 
         7.000e+00 6.000e+00 3.000e+00 ]
       irb(main):021:0> u, v, s = m.SV_decomp  # u, v: Matrix, s: Vector (singular values)
       irb(main):022:0> u*u.trans              # u is orthnormal
       => GSL::Matrix: 
       [  1.000e+00  2.452e-17 -4.083e-16 
          2.452e-17  1.000e+00 -3.245e-16 
         -4.083e-16 -3.245e-16  1.000e+00 ]
       irb(main):023:0> v*v.trans              # v is also orthnormal
       => GSL::Matrix: 
       [  1.000e+00  3.555e-17 -1.867e-16 
          3.555e-17  1.000e+00 -1.403e-16 
         -1.867e-16 -1.403e-16  1.000e+00 ]
       irb(main):024:0> u*Matrix.diagonal(s)*v.trans # Reconstruct the matrix
       => GSL::Matrix: 
       [ 3.000e+00 5.000e+00 2.000e+00 
         5.000e+00 1.000e+00 4.000e+00 
         7.000e+00 6.000e+00 3.000e+00 ]

--- GSL::Linalg::SV.decomp_mod(A)
--- GSL::Matrix#SV_decomp_mod
    These compute the SVD using the modified Golub-Reinsch algorithm, 
    which is faster for M>>N.

--- GSL::Linalg::SV.decomp_jacobi(A)
--- GSL::Matrix#SV_decomp_jacobi
    These compute the SVD using one-sided Jacobi orthogonalization. 
    The Jacobi method can compute singular values to higher relative accuracy 
    than Golub-Reinsch algorithms.

--- GSL::Linalg::SV.solve(A, b)
--- GSL::Linalg::SV.solve(U, V, S, b)
--- GSL::Matrix#SV_solve(b)
    These methods solve the system (({A x = b})) using the singular value 
    decomposition ((|U, S, V|)) of ((|A|)).

    * Ex1:
        m = Matrix.alloc(...)
        b = Vector.alloc(...)
        u, v, s = GSL::Linalg::SV.decomp(m)
        x = GSL::Linalg::SV.solve(u, v, s, b)
    * Ex2:
        x = GSL::Linalg::SV.solve(m, b)
    * Ex3:
        x = m.SV_solve(b)

== Cholesky Decomposition
A symmetric, positive definite square matrix ((|A|)) has a Cholesky decomposition 
into a product of a lower triangular matrix L and its transpose L^T, 
as ((|A = L L^T|)). This is sometimes referred to as taking the square-root of a 
matrix. The Cholesky decomposition can only be carried out when all the eigenvalues 
of the matrix are positive. This decomposition can be used to convert the linear 
system ((|A x = b|)) into a pair of triangular systems (((|L y = b, L^T x = y|))), 
which can be solved by forward and back-substitution.

--- GSL::Linalg::Cholesky.decomp(A)
    This method factorizes the positive-definite square matrix ((|A|)) 
    into the Cholesky decomposition ((|A = L L^T|)). 
    The upper triangular part of the matrix returned contains L^T, the diagonal terms 
    being identical for both L and L^T. If the matrix is not positive-definite 
    then the decomposition will fail.

    Ex: 
      irb(main):006:0> m = Matrix.pascal(3)
      => GSL::Matrix
      [  1.000e+00  1.000e+00  1.000e+00 
         1.000e+00  2.000e+00  3.000e+00 
         1.000e+00  3.000e+00  6.000e+00 ]
      irb(main):007:0> c = Linalg::Cholesky.decomp(m)
      => GSL::Linalg::Cholesky::CholeskyMatrix
      [  1.000e+00  1.000e+00  1.000e+00 
         1.000e+00  1.000e+00  2.000e+00 
         1.000e+00  2.000e+00  1.000e+00 ]
      irb(main):008:0> l = c.lower
      => GSL::Matrix
      [  1.000e+00  0.000e+00  0.000e+00 
         1.000e+00  1.000e+00  0.000e+00 
         1.000e+00  2.000e+00  1.000e+00 ]
      irb(main):009:0> (l*l.trans) == m
      => true

--- GSL::Linalg::Cholesky.solve(cholesky, b)
--- GSL::Linalg::Cholesky.svx(cholesky, x)
    These methods solve the system ((|A x = b|)) using the Cholesky decomposition 
    of ((|A|)) into the matrix ((|cholesky|)) given by (({GSL::Linalg::Cholesky.decomp})).

== Tridiagonal Decomposition of Real Symmetric Matrices 
--- GSL::Linalg::Symmtd::decomp(A)
    Factorizes the symmetric square matrix ((|A|)) into the symmetric 
    tridiagonal decomposition Q T Q^T, and returns the results
    ((|(A', tau)|)). On output the diagonal and subdiagonal part of the 
    matrix ((|A'|))  contain the tridiagonal matrix ((|T|)). 
    The remaining lower triangular part of the matrix ((|A'|)) contains 
    the Householder vectors which, together with the Householder 
    coefficients ((|tau|)), encode the orthogonal matrix ((|Q|)). 
    This storage scheme is the same as used by LAPACK. 
    The upper triangular part of ((|A|)) is not referenced. 
--- GSL::Linalg::Symmtd::unpack(A', tau)
    Unpacks the encoded symmetric tridiagonal decomposition ((|(A', tau)|)) 
    obtained from (({GSL::Linalg::Symmtd::decomp})) into the orthogonal 
    matrix ((|Q|)), the vector of diagonal elements ((|diag|)) 
    and the vector of subdiagonal elements ((|subdiag|)). 
--- GSL::Linalg::Symmtd::unpack_T(A', tau)
    Unpacks the diagonal and subdiagonal of the encoded symmetric tridiagonal 
    decomposition ((|(A', tau)|)) obtained from 
    (({GSL::Linalg::Symmtd::decomp})) into the vectors 
    ((|diag|)) and ((|subdiag|)). 

== Tridiagonal Decomposition of Hermitian Matrices 
--- GSL::Linalg::Hermtd::decomp(A)
    Factorizes the hermitian matrix ((|A|)) into the symmetric tridiagonal 
    decomposition U T U^T, and returns the result as ((|(A', tau)|)). 
    On output the real parts of the diagonal and subdiagonal part of the
    matrix ((|A'|)) contain the tridiagonal matrix ((|T|)). 
    The remaining lower triangular part of the matrix ((|A'|)) contains 
    the Householder vectors which, together with the Householder 
    coefficients ((|tau|)), encode the orthogonal matrix ((|Q|)). 
    This storage scheme is the same as used by LAPACK. 
    The upper triangular part of ((|A|)) and imaginary parts of the diagonal 
    are not referenced. 

--- GSL::Linalg::Hermtd::unpack(A', tau)
    Unpacks the encoded tridiagonal decomposition ((|(A', tau)|)) 
    obtained from (({GSL::Linalg::Hermtd::decomp})) into the unitary matrix 
    ((|U|)), the real vector of diagonal elements ((|diag|)) and 
    the real vector of subdiagonal elements ((|subdiag|)). 

--- GSL::Linalg::Hermtd::unpack_T(A', tau)
    Unpacks the diagonal and subdiagonal of the encoded tridiagonal 
    decomposition ((|(A, tau)|)) obtained from the 
    (({GSL::Linalg::Hermtd::decomp})) 
    into the real vectors ((|diag|)) and ((|subdiag|)). 

== Hessenberg Decomposition of Real Matrices 
--- GSL::Linalg::Hessenberg::decomp(A)
--- GSL::Linalg::hessenberg_decomp(A)
    Computes the Hessenberg decomposition of the matrix ((|A|)) 
    by applying the similarity transformation ((|H = U^T A U|)), and returns
    the result as ((|(A', tau|)). On output, ((|H|)) is stored in the upper 
    portion of ((|A'|)). The information required to construct the matrix 
    ((|U|)) is stored in the lower triangular portion of ((|A'|)). 
    ((|U|)) is a product of N - 2 Householder matrices. 
    The Householder vectors are stored in the lower portion of ((|A'|)) 
    (below the subdiagonal) and the Householder coefficients are stored 
    in the vector ((|tau|)).

--- GSL::Linalg::Hessenberg::unpack(A', tau)
--- GSL::Linalg::hessenberg_unpack(A', tau)
    Constructs the orthogonal matrix ((|U|)) and returns it 
    from the information stored  in the Hessenberg matrix ((|A'|)) 
    along with the vector ((|tau|)). ((|A'|)) and ((|tau|)) 
    are outputs from (({GSL::Linalg::Hessenberg::decomp})).

--- GSL::Linalg::Hessenberg::unpack_accum(A', tau, V = I)
--- GSL::Linalg::hessenberg_unpack_accum(A', tau, V = I)
    This method is similar to (({GSL::Linalg::Hessenberg::unpack})), 
    except it accumulates the matrix ((|U|)) into ((|V|)), so that 
    ((|V' = VU|)), and returns ((|V|)). Setting V to the identity matrix 
    provides the same result (({GSL::Linalg::Hessenberg::unpack})). 

--- GSL::Linalg::Hessenberg::set_zero(A')
--- GSL::Linalg::hessenberg_set_zero(A')
    Sets the lower triangular portion of ((|A'|)), below the subdiagonal, 
    to zero. 
    It is useful for clearing out the Householder vectors after calling 
    (({GSL::Linalg::Hessenberg::decomp})).

== Hessenberg-Triangular Decomposition of Real Matrices 
-- GSL::Linalg::hesstri_decomp(A, B)
-- GSL::Linalg::hesstri_decomp(A, B, work)
-- GSL::Linalg::hesstri_decomp(A, B, U, V)
-- GSL::Linalg::hesstri_decomp(A, B, U, V, work)
   Compute the Hessenberg-Triangular decomposition of the matrix pair 
   ((|(A, B)|)), and return ((|(H, R|)).
   If U and V are provided (they may be null), the similarity 
   transformations are stored in them. ((|work|)) is an additional workspace
   of length ((|N|)).

-- GSL::Linalg::hesstri_decomp!(A, B)
-- GSL::Linalg::hesstri_decomp!(A, B, work)
-- GSL::Linalg::hesstri_decomp!(A, B, U, V)
-- GSL::Linalg::hesstri_decomp!(A, B, U, V, work)
   Compute the Hessenberg-Triangular decomposition of the matrix pair 
   ((|(A, B)|)). On output, ((|H|)) is stored in ((|A|)), 
   and ((|R|)) is stored in ((|B|)).
   If U and V are provided (they may be null), the similarity 
   transformations are stored in them. ((|work|)) is an additional workspace
   of length ((|N|)).

== Bidiagonalization
--- GSL::Linalg::Bidiag::decomp!(A)
--- GSL::Linalg::bidiag_decomp!(A)
--- GSL::Linalg::Bidiag::decomp(A)
--- GSL::Linalg::bidiag_decomp(A)

--- GSL::Linalg::Bidiag::unpack
--- GSL::Linalg::bidiag_unpack
--- GSL::Linalg::Bidiag::unpack2
--- GSL::Linalg::bidiag_unpack2
--- GSL::Linalg::Bidiag::unpack_B
--- GSL::Linalg::bidiag_unpack_B

== Householder Transformations
--- GSL::Linalg::Householder::transform(v)
--- GSL::Linalg::HH::transform(v)
--- GSL::Vector#householder_transform
    These methods prepare a Householder transformation P = I - tau v v^T 
    which can be used to zero all the elements of the input vector except the first. 
    On output the transformation is stored in the vector ((|v|)) 
    and the scalar tau is returned.

--- GSL::Linalg::Householder::hm(tau, v, A)
--- GSL::Linalg::HH::hm(tau, v, A)
    These methods apply the Householder matrix P defined by the scalar 
    ((|tau|)) and the vector ((|v|)) to the left-hand side of the matrix ((|A|)). 
    On output the result P A is stored in ((|A|)).

--- GSL::Linalg::Householder::mh(tau, v, A)
--- GSL::Linalg::HH::mh(tau, v, A)
    These methods apply the Householder matrix P defined by the scalar ((|tau|)) 
    and the vector ((|v|)) to the right-hand side of the matrix ((|A|)). 
    On output the result A P is stored in ((|A|)).

--- GSL::Linalg::Householder::hv(tau, v, w)
--- GSL::Linalg::HH::hv(tau, v, w)
    These methods apply the Householder transformation P defined by the scalar 
    ((|tau|)) and the vector ((|v|)) to the vector ((|w|)). 
    On output the result P w is stored in ((|w|)).

== Householder solver for linear systems
--- GSL::Linalg::HH.solve(A, b)
--- GSL::Matrix#HH_solve(b)
    These methods solve the system (({A x = b})) directly using Householder 
    transformations. The matrix ((|A|)) is not modified.
--- GSL::Linalg::HH.solve!(A, b)
--- GSL::Matrix#HH_solve!(b)
    These methods solve the system (({A x = b})) directly using Householder 
    transformations. The matrix ((|A|)) is destroyed by the 
    Householder transformations.

--- GSL::Linalg::HH.svx(A, b)
--- GSL::Matrix#HH_svx(b)
    These methods solve the system (({A x = b})) in-place directly using Householder 
    transformations. The input vector ((|b|)) is replaced by the solution.

== Tridiagonal Systems
--- GSL::Linglg::solve_tridiag(diag, e, f, b)
--- GSL::Linglg::Tridiag::solve(diag, e, f, b)
    These methods solve the general N-by-N system A x = b where ((|A|)) 
    is tridiagonal ( N >= 2). The super-diagonal and sub-diagonal vectors ((|e|)) 
    and ((|f|)) must be one element shorter than the diagonal vector ((|diag|)). 
    The form of ((|A|)) for the 4-by-4 case is shown below,
         A = ( d_0 e_0  0   0  )
             ( f_0 d_1 e_1  0  )
             (  0  f_1 d_2 e_2 )
             (  0   0  f_2 d_3 )

--- GSL::Linglg::solve_symm_tridiag(diag, e, b)
--- GSL::Linglg::Tridiag::solve_symm(diag, e, b)
    These methods solve the general N-by-N system A x = b where ((|A|)) is 
    symmetric tridiagonal ( N >= 2). The off-diagonal vector ((|e|)) must 
    be one element shorter than the diagonal vector ((|diag|)). 
    The form of ((|A|)) for the 4-by-4 case is shown below,
         A = ( d_0 e_0  0   0  )
             ( e_0 d_1 e_1  0  )
             (  0  e_1 d_2 e_2 )
             (  0   0  e_2 d_3 )

--- GSL::Linglg::solve_cyc_tridiag(diag, e, f, b)
--- GSL::Linglg::Tridiag::solve_cyc(diag, e, f, b)
    These methods solve the general N-by-N system A x = b where ((|A|)) is 
    cyclic tridiagonal ( N >= 3). The cyclic super-diagonal and sub-diagonal 
    vectors ((|e|)) and ((|f|)) must have the same number of elements as the 
    diagonal vector ((|diag|)). The form of ((|A|)) for the 4-by-4 case is shown below,
         A = ( d_0 e_0  0  f_3 )
             ( f_0 d_1 e_1  0  )
             (  0  f_1 d_2 e_2 )
             ( e_3  0  f_2 d_3 )

--- GSL::Linglg::solve_symm_cyc_tridiag(diag, e, b)
--- GSL::Linglg::Tridiag::solve_symm_cyc(diag, e, b)
    These methods solve the general N-by-N system A x = b where ((|A|)) 
    is symmetric cyclic tridiagonal ( N >= 3). The cyclic off-diagonal vector ((|e|)) 
    must have the same number of elements as the diagonal vector ((|diag|)). 
    The form of ((|A|)) for the 4-by-4 case is shown below,
         A = ( d_0 e_0  0  e_3 )
             ( e_0 d_1 e_1  0  )
             (  0  e_1 d_2 e_2 )
             ( e_3  0  e_2 d_3 )

== Balancing
The process of balancing a matrix applies similarity transformations to 
make the rows and columns have comparable norms. This is useful, 
for example, to reduce roundoff errors in the solution of eigenvalue problems. 
Balancing a matrix ((|A|)) consists of replacing ((|A|)) with a similar matrix 
where ((|D|)) is a diagonal matrix whose entries are powers of the floating 
point radix. 

--- GSL::Linalg::balance(A)
--- GSL::Linalg::balance(A, D)
    Calculates the balanced counterpart of ((|A|)) and the diagonal elements 
    of the similarity transformation. The result is returned as ((|(A', D)|)).

--- GSL::Linalg::balance!(A)
--- GSL::Linalg::balance!(A, D)
    Replaces the matrix ((|A|)) with its balanced counterpart and 
    stores the diagonal elements of the similarity transformation into 
    the vector ((|D|)). 


== NArray
The following Linalg methods can handle NArray objects:
* GSL::Linalg::
  * LU::
    * LU.decomp(m)
    * LU.solve(lu, b)
    * LU.svx(lu, bx)
    * LU.det(lu, sign)
    * LU.lndet(lu)
    * LU.invert(lu, perm)
  * QR::
    * QR.decomp(m)
    * QR.solve(qr, tau, b)	
    * QR.svx(qr, tau, bx)
  * SV::
    * SV.decomp(m)
    * SV.solve(u, v, s, b)
    * SV.svx(u, v, s, bx)
  * Cholesky::
    * Cholesky.decomp(m)
    * Cholesky.solve(u, v, s, b)
    * Cholesky.svx(u, v, s, bx)
  * HH::
    * HH.solve(m, b)
    * HH.svx(m, bx)

((<prev|URL:blas.html>))
((<next|URL:eigen.html>))

((<Reference index|URL:ref.html>))
((<top|URL:index.html>))

=end
