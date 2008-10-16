=begin
= Discrete Hankel Transforms
This chapter describes functions for performing Discrete Hankel Transforms 
(DHTs). 

(1) ((<Definitions|URL:dht.html#1>))
(2) ((<Initialization|URL:dht.html#2>))
(3) ((<Methods|URL:dht.html#3>))

== Definitions
The discrete Hankel transform acts on a vector of sampled data, where the 
samples are assumed to have been taken at points related to the zeroes of a 
Bessel function of fixed order; compare this to the case of the discrete 
Fourier transform, where samples are taken at points related to the zeroes 
of the sine or cosine function. 

Specifically, let f(t) be a function on the unit interval. Then the finite 
\nu-Hankel transform of f(t) is defined to be the set of numbers g_m given by, 
so that, Suppose that f is band-limited in the sense that g_m=0 for m > M.
Then we have the following fundamental sampling theorem. It is this discrete 
expression which defines the discrete Hankel transform. The kernel in the 
summation above defines the matrix of the \nu-Hankel transform of size M-1. 
The coefficients of this matrix, being dependent on \nu and M, must be 
precomputed and stored; the (({GSL::Dht})) object encapsulates this data. 
The constructor (({GSL::Dht.alloc})) returns a (({GSL::Dht})) object 
which must be properly initialized with (({GSL::Dht#init})) before 
it can be used to perform transforms on data sample vectors, 
for fixed \nu and M, using the (({GSL::Dht#apply})) method. 
The implementation allows a scaling of the fundamental 
interval, for convenience, so that one can assume the function is defined on 
the interval [0,X], rather than the unit interval. 

Notice that by assumption f(t) vanishes at the endpoints of the interval, 
consistent with the inversion formula and the sampling formula given above. 
Therefore, this transform corresponds to an orthogonal expansion in 
eigenfunctions of the Dirichlet problem for the Bessel differential equation. 


== Initialization

--- GSL::Dht.alloc(size)
--- GSL::Dht.alloc(size, nu, xmax)
    These methods allocate a Discrete Hankel transform object (({GSL::Dht})) 
    of size ((|size|)).
    If three arguments are given, the object is initialized with the values of
    ((|nu, xmax|)).

--- GSL::Dht#init(nu, xmax)
    This initializes the transform ((|self|)) for the given values of ((|nu|)) and ((|xmax|)).

== Methods
--- GSL::Dht#apply(vin, vout)
--- GSL::Dht#apply(vin)
    This applies the transform ((|self|)) to the vector ((|vin|)) whose size is 
    equal to the size of the transform.

--- GSL::Dht#x_sample(n)
    This method returns the value of the n'th sample point in the unit interval, 
    (j_{nu,n+1}/j_{nu,M}) X. These are the points where the function f(t) is 
    assumed to be sampled.

--- GSL::Dht#k_sample(n)
    This method returns the value of the n'th sample point in "k-space", 
    j_{nu,n+1}/X.

--- GSL::Dht#size
    Returns the size of the sample arrays to be transformed
--- GSL::Dht#nu
    Returns the Bessel function order
--- GSL::Dht#xmax
    Returns the upper limit to the x-sampling domain 
--- GSL::Dht#kmax
    Returns the upper limit to the k-sampling domain 

--- GSL::Dht#j
    Returns an array of computed J_nu zeros, j_{nu,s} = j[s] 
    as a (({GSL::Vector::View})).

--- GSL::Dht#Jjj
    Returns an array of transform numerator, J_nu(j_i j_m / j_N)
    as a (({GSL::Vector::View})).

--- GSL::Dht#J2
    Returns an array of transform numerator, J_nu(j_i j_m / j_N).

--- GSL::Dht#coef
--- GSL::Dht#coef(n, m)
    Return the (n,m)-th transform coefficient.

((<prev|URL:sum.html>))
((<next|URL:roots.html>))

((<Reference index|URL:ref.html>))
((<top|URL:index.html>))

=end
