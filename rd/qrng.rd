=begin
= Quasi-Random Sequences

This chapter describes the quasi-random sequence generator (({GSL::QRng}))
of arbitrary dimensions. A quasi-random sequence progressively covers 
a d-dimensional space with a set of points that are uniformly distributed. 
Quasi-random sequences are also known as low-discrepancy sequences.
The quasi-random sequence generators use an interface that is similar 
to the interface for random number generators.

Contents:
(1) ((<Quasi-random number generator initialization|URL:qrng.html#1>))
(2) ((<Sampling from a quasi-random number generator|URL:qrng.html#2>))
(3) ((<Auxiliary quasi-random number generator functions|URL:qrng.html#3>))
(4) ((<Saving and resorting quasi-random number generator state|URL:qrng.html#4>))
(5) ((<Quasi-random number generator algorithms|URL:qrng.html#5>))

== Quasi-random number generator initialization

--- GSL::QRng.alloc(T, d)
    This returns a GSL::QRng object, a quasi-random sequence generator of type ((|T|)) and dimension ((|d|)). 

--- GSL::QRng::init
    This reinitializes the generator to its starting point. 

== Sampling from a quasi-random number generator 

--- GSL::QRng::get(x)
    This calculate the next point ((|x|)) from the sequence generator. Here ((|x|)) is an instance of the ((<GSL::Vector|URL:vector.html>)) class. The space available for ((|x|)) must match the dimension of the generator. The point ((|x|)) will lie in the range 0 < x_i < 1 for each x_i. 

    This is used as
      q = QRng.alloc(QRng::SOBOL, dim)
      v = Vector.alloc(dim)
      for i in 0..1024 do
        q.get(v)
        printf("%.5f %.5f\n", v[0], v[1])
      end

== Auxiliary quasi-random number generator functions 

--- GSL::QRng::name
    Returns the name of the generator ((|self|)).

--- GSL::QRng::size

== Saving and resorting quasi-random number generator state
--- GSL::QRng::clone
--- GSL::QRng::duplicate
    Return a newly created generator which is an exact copy of the generator ((|self|)).
== Quasi-random number generator algorithms
In creating a generator by the method (({GSL::QRng.alloc(T, d)})), 
the algorithm type ((|T|)) is given by a String or a Fixnum constant.
The following quasi-random sequence algorithms are available,

* "(({niederreiter_2}))" (String)
* (({GSL::QRng::NIEDERREITER_2})) (Fixnum)

  The generator of this type uses the algorithm described in Bratley, Fox, Niederreiter, ACM Trans. Model. Comp. Sim. 2, 195 (1992). It is valid up to 12 dimensions.

* "(({sobol}))" (String)
* (({GSL::QRng::SOBOL})) (Fixnum)

  This generator uses the Sobol sequence described in Antonov, Saleev, USSR Comput. Maths. Math. Phys. 19, 252 (1980). It is valid up to 40 dimensions.

((<prev|URL:rng.html>))
((<next|URL:randist.html>))

((<Reference index|URL:ref.html>))
((<top|URL:index.html>))

=end
