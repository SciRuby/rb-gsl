=begin

== Cholesky decomposition (>= GSL-1.10)
A symmetric, positive definite square matrix ((|A|)) has 
a Cholesky decomposition into a product of a lower triangular matrix 
((|L|)) and its transpose ((|L^T|)).
This is sometimes referred to as taking the square-root of a matrix. 
The Cholesky decomposition can only be carried out when all the eigenvalues 
of the matrix are positive. This decomposition can be used to convert the 
linear system ((|A x = b|)) into a pair of triangular systems 
((|L y = b, L^T x = y|)),
which can be solved by forward and back-substitution. 

--- GSL::Linalg::Complex::Cholesky::decomp(A)
--- GSL::Linalg::Complex::cholesky_decomp(A)
    Factorize the positive-definite square matrix ((|A|)) into the 
    Cholesky decomposition ((|A = L L^H|)).
    On input only the diagonal and lower-triangular part of the matrix ((|A|)) 
    are needed. The diagonal and lower triangular part of the returned matrix
    contain the matrix ((|L|)). The upper triangular part of the 
    returned matrix contains L^T, and
    the diagonal terms being identical for both L and L^T. 
    If the input matrix is not positive-definite then the decomposition 
    will fail, returning the error code (({GSL::EDOM})). 

--- GSL::Linalg::Complex::Cholesky::solve(chol, b, x)
--- GSL::Linalg::Complex::cholesky_solve(chol, b, x)
    Solve the system ((|A x = b|)) using the Cholesky decomposition 
    of ((|A|)) into the matrix ((|chol|)) given by 
    (({GSL::Linalg::Complex::Cholesky::decomp})).

--- GSL::Linalg::Complex::Cholesky::svx(chol, x)
--- GSL::Linalg::Complex::cholesky_svx(chol, x)
    Solve the system ((|A x = b|)) in-place using the Cholesky decomposition 
    of ((|A|)) into the matrix ((|chol|)) given by 
    (({GSL::Linalg::Complex::Cholesky::decomp})). On input ((|x|)) 
    should contain the right-hand side ((|b|)), 
    which is replaced by the solution on output. 

((<back|URL:linalg.html>))
=end
