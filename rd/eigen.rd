=begin
= Eigensystems
=== Contentes
(1) ((<Modules and classes|URL:eigen.html#1>))
(2) ((<Real Symmetric Matrices|URL:eigen.html#2>))
(3) ((<Complex Hermitian Matrices|URL:eigen.html#3>))
(4) ((<Real Nonsymmetric Matrices|URL:eigen.html#4>)) (>= GSL-1.9)
(5) ((<Real Generalized Symmetric-Definite Eigensystems|URL:eigen.html#5>)) (>= GSL-1.10)
(6) ((<Complex Generalized Hermitian-Definite Eigensystems|URL:eigen.html#6>)) (>= GSL-1.10)
(7) ((<Real Generalized Nonsymmetric Eigensystems|URL:eigen.html#7>)) (>= GSL-1.10)
(8) ((<Sorting Eigenvalues and Eigenvectors |URL:eigen.html#8>))

== Modules and classes

* GSL
  * Eigen
    * EigenValues < Vector
    * EigenVectors < Matrix
    * Symm (Module)
      * Workspace (Class)
    * Symmv (Module)
      * Workspace (Class)
    * Nonsymm (Module, >= GSL-1.9)
      * Workspace (Class)
    * Nonsymmv (Module, >= GSL-1.9)
      * Workspace (Class)
    * Gensymm (Module, >= GSL-1.10)
      * Workspace (Class)
    * Gensymmv (Module, >= GSL-1.10)
      * Workspace (Class)
    * Herm (Module)
      * Workspace (Class)
    * Hermv (Module)
      * Workspace (Class)
      * Vectors < Matrix::Complex
    * Genherm (Module, >= GSL-1.10)
      * Workspace (Class)
    * Genhermv (Module, >= GSL-1.10)
      * Workspace (Class)
    * Gen (Module, >= GSL-1.10)
      * Workspace (Class)
    * Genv (Module, >= GSL-1.10)
      * Workspace (Class)

== Real Symmetric Matrices, GSL::Eigen::Symm module
=== Workspace classes
--- GSL::Eigen::Symm::Workspace.alloc(n)
--- GSL::Eigen::Symmv::Workspace.alloc(n)
--- GSL::Eigen::Herm::Workspace.alloc(n)
--- GSL::Eigen::Hermv::Workspace.alloc(n)

=== Methods to solve eigensystems
--- GSL::Eigen::symm(A)
--- GSL::Eigen::symm(A, workspace)
--- GSL::Matrix#eigen_symm
--- GSL::Matrix#eigen_symm(workspace)
    These methods compute the eigenvalues of the real symmetric matrix. 
    The workspace object ((|workspace|)) can be omitted.

--- GSL::Eigen::symmv(A)
--- GSL::Matrix#eigen_symmv
    These methods compute the eigenvalues and eigenvectors of the real symmetric 
    matrix, and return an array of two elements:
    The first is a (({GSL::Vector})) object which stores all the eigenvalues. 
    The second is a (({GSL::Matrix object})), whose columns contain 
    eigenvectors.

    (1) Singleton method of the (({GSL::Eigen})) module, (({GSL::Eigen::symm}))

          m = GSL::Matrix.alloc([1.0, 1/2.0, 1/3.0, 1/4.0], [1/2.0, 1/3.0, 1/4.0, 1/5.0],
                             [1/3.0, 1/4.0, 1/5.0, 1/6.0], [1/4.0, 1/5.0, 1/6.0, 1/7.0])
          eigval, eigvec = Eigen::symmv(m)

    (2) Instance method of (({GSL::Matrix})) class

          eigval, eigvec = m.eigen_symmv

== Complex Hermitian Matrices
--- GSL::Eigen::herm(A)
--- GSL::Eigen::herm(A, workspace)
--- GSL::Matrix::Complex#eigen_herm
--- GSL::Matrix::Complex#eigen_herm(workspace)
    These methods compute the eigenvalues of the complex hermitian matrix. 

--- GSL::Eigen::hermv(A)
--- GSL::Eigen::hermv(A, workspace)
--- GSL::Matrix::Complex#eigen_hermv
--- GSL::Matrix::Complex#eigen_hermv(workspace

== Real Nonsymmetric Matrices (>= GSL-1.9)

--- GSL::Eigen::Nonsymm.alloc(n)
    This allocates a workspace for computing eigenvalues of n-by-n real 
    nonsymmetric matrices. The size of the workspace is O(2n).

--- GSL::Eigen::Nonsymm::params(compute_t, balance, wspace)
--- GSL::Eigen::Nonsymm::Workspace#params(compute_t, balance)
    This method sets some parameters which determine how the eigenvalue 
    problem is solved in subsequent calls to (({GSL::Eigen::nonsymm})).
    If ((|compute_t|)) is set to 1, the full Schur form (({T})) will be 
    computed by gsl_eigen_nonsymm. If it is set to 0, (({T})) will not be 
    computed (this is the default setting). 
    Computing the full Schur form (({T})) requires approximately 1.5-2 times 
    the number of flops.

    If ((|balance|)) is set to 1, a balancing transformation is applied to 
    the matrix prior to computing eigenvalues. This transformation is designed 
    to make the rows and columns of the matrix have comparable norms, and can 
    result in more accurate eigenvalues for matrices whose entries vary widely 
    in magnitude. See section Balancing for more information. Note that the
    balancing transformation does not preserve the orthogonality of the Schur 
    vectors, so if you wish to compute the Schur vectors with 
    (({GSL::Eigen::nonsymm_Z})) you will obtain the Schur vectors of the 
    balanced matrix instead of the original matrix. The relationship will be 
    where Q is the matrix of Schur vectors for the balanced matrix, and (({D})) 
    is the balancing transformation. Then (({GSL::Eigen::nonsymm_Z})) will 
    compute a matrix (({Z})) which satisfies with (({Z = D Q})). 
    Note that (({Z})) will not be orthogonal. For this reason, balancing is
    not performed by default.

--- GSL::Eigen::nonsymm(m, eval, wspace)
--- GSL::Eigen::nonsymm(m)
--- GSL::Matrix#eigen_nonsymm()
--- GSL::Matrix#eigen_nonsymm(wspace)
--- GSL::Matrix#eigen_nonsymm(eval, wspace)
    These methods compute the eigenvalues of the real nonsymmetric matrix (({m})) 
    and return them, or store in the vector ((|eval|)) if it given. 
    If (({T})) is desired, it is stored in (({m})) on output, however the lower 
    triangular portion will not be zeroed out. Otherwise, on output, the diagonal 
    of (({m})) will contain the 1-by-1 real eigenvalues and 2-by-2 complex 
    conjugate eigenvalue systems, and the rest of (({m})) is destroyed. 

--- GSL::Eigen::nonsymm_Z(m, eval, Z, wspace)
--- GSL::Eigen::nonsymm_Z(m)
--- GSL::Matrix#eigen_nonsymm_Z()
--- GSL::Matrix#eigen_nonsymm(eval, Z, wspace)
    These methods are identical to (({GSL::Eigen::nonsymm})) except they also 
    compute the Schur vectors and return them (or store into (({Z}))).

--- GSL::Eigen::Nonsymmv.alloc(n)
    Allocates a workspace for computing eigenvalues and eigenvectors 
    of n-by-n real nonsymmetric matrices. The size of the workspace is O(5n).
--- GSL::Eigen::nonsymm(m)
--- GSL::Eigen::nonsymm(m, wspace)
--- GSL::Eigen::nonsymm(m, eval, evec)
--- GSL::Eigen::nonsymm(m, eval, evec, wspace)
--- GSL::Matrix#eigen_nonsymmv()
--- GSL::Matrix#eigen_nonsymmv(wspace)
--- GSL::Matrix#eigen_nonsymmv(eval, evec)
--- GSL::Matrix#eigen_nonsymmv(eval, evec, wspace)
    Compute eigenvalues and right eigenvectors of the n-by-n real nonsymmetric 
    matrix. The computed eigenvectors are normalized to have Euclidean norm 1.
    On output, the upper portion of ((|m|)) contains the Schur form ((|T|)). 

== Real Generalized Symmetric-Definite Eigensystems (GSL-1.10)
The real generalized symmetric-definite eigenvalue problem is to 
find eigenvalues ((|lambda|)) and eigenvectors ((|x|)) such that 
where ((|A|)) and ((|B|)) are symmetric matrices, and ((|B|)) 
is positive-definite. This problem reduces to the standard symmetric eigenvalue
problem by applying the Cholesky decomposition to ((|B|)): 
Therefore, the problem becomes ((|C y = lambda y|)) 
where ((|C = L^{-1} A L^{-t}|)) is symmetric, and ((|y = L^t x|)). 
The standard symmetric eigensolver can be applied to the matrix ((|C|)). 
The resulting eigenvectors are backtransformed to find the vectors of the 
original problem. The eigenvalues and eigenvectors of the generalized 
symmetric-definite eigenproblem are always real. 

--- GSL::Eigen::Gensymm.alloc(n)
--- GSL::Eigen::Gensymm::Workspace.alloc(n)
    Allocates a workspace for computing eigenvalues of n-by-n real 
    generalized symmetric-definite eigensystems. 
    The size of the workspace is O(2n). 
--- GSL::Eigen::gensymm(A, B, w)
    Computes the eigenvalues of the real generalized symmetric-definite matrix 
    pair ((|A, B|)), and returns them as a (({GSL::Vector})), 
    using the method outlined above. On output, B contains its Cholesky
    decomposition.
--- GSL::Eigen::gensymmv(A, B, w)
    Computes the eigenvalues and eigenvectors of the real generalized 
    symmetric-definite matrix pair ((|A, B|)), and returns 
    them as a (({GSL::Vector})) and a (({GSL::Matrix})). 
    The computed eigenvectors are normalized to have unit magnitude. 
    On output, ((|B|)) contains its Cholesky decomposition.

== Complex Generalized Hermitian-Definite Eigensystems (>= GSL-1.10)
The complex generalized hermitian-definite eigenvalue problem is to 
find eigenvalues ((|lambda|)) and eigenvectors ((|x|)) such that 
where ((|A|)) and ((|B|)) are hermitian matrices, and ((|B|)) 
is positive-definite. Similarly to the real case, this can be reduced to 
((|C y = lambda y|)) where ((|C = L^{-1} A L^{-H}|)) is hermitian, 
and ((|y = L^H x|)). The standard hermitian eigensolver can be applied to 
the matrix ((|C|)). The resulting eigenvectors are backtransformed 
to find the vectors of the original problem. 
The eigenvalues of the generalized hermitian-definite eigenproblem are always 
real. 

--- GSL::Eigen::Genherm.alloc(n)
    Allocates a workspace for computing eigenvalues of n-by-n complex 
    generalized hermitian-definite eigensystems. 
    The size of the workspace is O(3n). 
--- GSL::Eigen::genherm(A, B, w)
    Computes the eigenvalues of the complex generalized hermitian-definite 
    matrix pair ((|A, B|)), and returns them as a (({GSL::Vector})), 
    using the method outlined above. 
    On output, ((|B|)) contains its Cholesky decomposition.
--- GSL::Eigen::genherm(A, B, w)
    Computes the eigenvalues and eigenvectors of the complex generalized 
    hermitian-definite matrix pair ((|A, B|)), 
    and returns them as a (({GSL::Vector})) and a (({GSL::Matrix::Complex})). 
    The computed eigenvectors are normalized to have unit magnitude. 
    On output, ((|B|)) contains its Cholesky decomposition.

== Real Generalized Nonsymmetric Eigensystems (>= GSL-1.10)

--- GSL::Eigen::Gen.alloc(n)
--- GSL::Eigen::Gen::Workspace.alloc(n)
    Allocates a workspace for computing eigenvalues of n-by-n real generalized
    nonsymmetric eigensystems. The size of the workspace is O(n).

--- GSL::Eigen::Gen::params(compute_s, compute_t, balance, w)
--- GSL::Eigen::gen_params(compute_s, compute_t, balance, w)
   Set some parameters which determine how the eigenvalue problem is solved 
   in subsequent calls to (({GSL::Eigen::gen})).

   If ((|compute_s|)) is set to 1, the full Schur form ((|S|)) will be 
   computed by (({GSL::Eigen::gen}). If it is set to 0, ((|S|)) will 
   not be computed (this is the default setting). ((|S|)) is a quasi upper 
   triangular matrix with 1-by-1 and 2-by-2 blocks on its diagonal. 
   1-by-1 blocks correspond to real eigenvalues, and 2-by-2 blocks 
   correspond to complex eigenvalues. 

   If ((|compute_t|)) is set to 1, the full Schur form ((|T|)) will 
   be computed by (({GSL::Eigen::gen}). If it is set to 0, ((|T|)) 
   will not be computed (this is the default setting). ((|T|)) 
   is an upper triangular matrix with non-negative elements on its diagonal. 
   Any 2-by-2 blocks in ((|S|)) will correspond to a 2-by-2 diagonal block 
   in ((|T|)). 

   The ((|balance|)) parameter is currently ignored, since generalized 
   balancing is not yet implemented. 

--- GSL::Eigen::gen(A, B, w)
    Computes the eigenvalues of the real generalized nonsymmetric matrix pair
    ((|A, B|)), and returns them as pairs in (alpha, beta), 
    where alpha is (({GSL::Vector::Complex})) and beta is (({GSL::Vector})). 
    If beta_i is non-zero, then lambda = alpha_i / beta_i is an eigenvalue.
    Likewise, if alpha_i is non-zero, then mu = beta_i / alpha_i is an 
    eigenvalue of the alternate problem mu A y = B y. 
    The elements of ((|beta|)) are normalized to be non-negative. 

    If ((|S|)) is desired, it is stored in ((|A|)) on output. 
    If ((|T|)) is desired, it is stored in ((|B|)) on output. 
    The ordering of eigenvalues in ((|alpha, beta|)) 
    follows the ordering of the diagonal blocks in the Schur forms ((|S|))
    and ((|T|)). 

--- GSL::Eigen::gen_QZ(A, B, w)
    This method is identical to (({GSL::Eigen::gen})) except it also computes 
    the left and right Schur vectors and returns them.

--- GSL::Eigen::Genv.alloc(n)
--- GSL::Eigen::Genv::Workspace.alloc(n)
    Allocatesa workspace for computing eigenvalues and eigenvectors of 
    n-by-n real generalized nonsymmetric eigensystems. 
    The size of the workspace is O(7n). 

--- GSL::Eigen::genv(A, B, w)
    Computes eigenvalues and right eigenvectors of the n-by-n real generalized 
    nonsymmetric matrix pair ((|A, B|)). The eigenvalues and eigenvectors 
    are returned in ((|alpha, beta, evec|)). 
    On output, ((|A, B|)) contains the generalized Schur form ((|S, T|)). 

--- GSL::Eigen::genv_QZ(A, B, w)
    This method is identical to (({GSL::Eigen::genv})) except it also computes 
    the left and right Schur vectors and returns them.

== Sorting Eigenvalues and Eigenvectors
--- GSL::Eigen::symmv_sort(eval, evec, type = GSL::Eigen::SORT_VAL_ASC)
--- GSL::Eigen::Symmv::sort(eval, evec, type = GSL::Eigen::SORT_VAL_ASC)
    These methods simultaneously sort the eigenvalues stored in the vector 
    ((|eval|)) and the corresponding real eigenvectors stored in the 
    columns of the matrix ((|evec|)) into ascending or descending order 
    according to the value of the parameter ((|type|)),
      * (({GSL::Eigen::SORT_VAL_ASC}))
        ascending order in numerical value
      * (({GSL::Eigen::SORT_VAL_DESC}))
        escending order in numerical value
      * (({GSL::Eigen::SORT_ABS_ASC}))
        scending order in magnitude
      * (({GSL::Eigen::SORT_ABS_DESC}))
        descending order in magnitude
    The sorting is carried out ((|in-place|)).

--- GSL::Eigen::hermv_sort(eval, evec, type = GSL::Eigen::SORT_VAL_ASC)
--- GSL::Eigen::Hermv::sort(eval, evec, type = GSL::Eigen::SORT_VAL_ASC)
    These methods simultaneously sort the eigenvalues stored in the vector 
    ((|eval|)) and the corresponding complex eigenvectors stored in the columns 
    of the matrix ((|evec|)) into ascending or descending order according 
    to the value of the parameter ((|type|)) as shown above.

--- GSL::Eigen::nonsymmv_sort(eval, evec, type = GSL::Eigen::SORT_VAL_ASC)
--- GSL::Eigen::Nonsymmv::sort(eval, evec, type = GSL::Eigen::SORT_VAL_ASC)
    Sorts the eigenvalues stored in the vector ((|eval|)) and the corresponding 
    complex eigenvectors stored in the columns of the matrix ((|evec|)) 
    into ascending or descending order according to the value of the 
    parameter ((|type|)) as shown above. 
    Only (({GSL::EIGEN_SORT_ABS_ASC})) and (({GSL::EIGEN_SORT_ABS_DESC})) 
    are supported due to the eigenvalues being complex. 

--- GSL::Eigen::gensymmv_sort(eval, evec, type = GSL::Eigen::SORT_VAL_ASC)
--- GSL::Eigen::Gensymmv::sort(eval, evec, type = GSL::Eigen::SORT_VAL_ASC)
    Sorts the eigenvalues stored in the vector ((|eval|)) and the 
    corresponding real eigenvectors stored in the columns of the matrix 
    ((|evec|)) into ascending or descending order according to the value of 
    the parameter ((|type|)) as shown above. 

--- GSL::Eigen::gensymmv_sort(eval, evec, type = GSL::Eigen::SORT_VAL_ASC)
--- GSL::Eigen::Gensymmv::sort(eval, evec, type = GSL::Eigen::SORT_VAL_ASC)
    Sorts the eigenvalues stored in the vector ((|eval|)) and the 
    corresponding complex eigenvectors stored in the columns of the matrix 
    ((|evec|)) into ascending or descending order according to the value of 
    the parameter ((|type|)) as shown above. 

--- GSL::Eigen::genv_sort(alpha, beta, evec, type = GSL::Eigen::SORT_VAL_ASC)
--- GSL::Eigen::Genv::sort(alpha, beta, evec, type = GSL::Eigen::SORT_VAL_ASC)
    Sorts the eigenvalues stored in the vectors ((|alpha, beta|)) and the 
    corresponding complex eigenvectors stored in the columns of the matrix 
    ((|evec|)) into ascending or descending order according to the value of 
    the parameter ((|type|)) as shown above. Only (({GSL::EIGEN_SORT_ABS_ASC}))
    and (({GSL::EIGEN_SORT_ABS_DESC})) are supported due to the eigenvalues 
    being complex. 

((<prev|URL:linalg.html>))
((<next|URL:fft.html>))

((<Reference index|URL:ref.html>))
((<top|URL:index.html>))

=end
