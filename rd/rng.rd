=begin
= Random Number Generation
The library provides a large collection of random number generators which 
can be accessed through a uniform interface. Environment variables allow you 
to select different generators and seeds at runtime, so that you can easily 
switch between generators without needing to recompile your program. 
Each instance of a generator keeps track of its own state, allowing the 
generators to be used in multi-threaded programs. Additional functions are 
available for transforming uniform random numbers into samples from 
continuous or discrete probability distributions such as the Gaussian, 
log-normal or Poisson distributions. 



Contents:
(1) ((<General comments on random numbers|URL:rng.html#1>))
(2) ((<The Random Number Generator Interface: GSL::Rng class|URL:rng.html#2>))
(3) ((<Random number generator initialization|URL:rng.html#3>))
(4) ((<Sampling from a random number generator|URL:rng.html#4>))
(5) ((<Auxiliary random number generator functions|URL:rng.html#5>))
(6) ((<Random number environment variables|URL:rng.html#6>))

== General comments on random numbers
In 1988, Park and Miller wrote a paper entitled "Random number generators: 
good ones are hard to find." [Commun. ACM, 31, 1192-1201]. Fortunately, some
excellent random number generators are available, though poor ones are still 
in common use. You may be happy with the system-supplied random number 
generator on your computer, but you should be aware that as computers get 
faster, requirements on random number generators increase. Nowadays, a 
simulation that calls a random number generator millions of times can often 
finish before you can make it down the hall to the coffee machine and back. 

A very nice review of random number generators was written by Pierre L'Ecuyer, 
as Chapter 4 of the book: Handbook on Simulation, Jerry Banks, ed.
(Wiley, 1997). The chapter is available in postscript from L'Ecuyer's 
ftp site (see references). Knuth's volume on Seminumerical Algorithms 
(originally published in 1968) devotes 170 pages to random number generators,
and has recently been updated in its 3rd edition (1997). It is brilliant, 
a classic. If you don't own it, you should stop reading right now, run to the 
nearest bookstore, and buy it. 

A good random number generator will satisfy both theoretical and statistical 
properties. Theoretical properties are often hard to obtain (they require real 
math!), but one prefers a random number generator with a long period, 
low serial correlation, and a tendency not to "fall mainly on the planes." 
Statistical tests are performed with numerical simulations. Generally,
a random number generator is used to estimate some quantity for which the 
theory of probability provides an exact answer. Comparison to this exact 
answer provides a measure of "randomness". 

== The Random Number Generator Interface 
It is important to remember that a random number generator is not a "real" 
function like sine or cosine. Unlike real functions, successive calls to a 
random number generator yield different return values. Of course that is just 
what you want for a random number generator, but to achieve this effect, 
the generator must keep track of some kind of "state" variable. Sometimes this 
state is just an integer (sometimes just the value of the previously generated 
random number), but often it is more complicated than that and may involve a 
whole array of numbers, possibly with some indices thrown in. To use the 
random number generators, you do not need to know the details of what 
comprises the state, and besides that varies from algorithm to algorithm. 

The random number generator library uses (({GSL::Rng})) class for the interface.
== Random number generator initialization 

--- GSL::Rng.alloc(rng_type[, seed])
    This method returns a GSL::Rng object of a random number generator of type 
    ((|rng_type|)) with a seed ((|seed|)). These two arguments can be omitted, 
    and the generator 'gsl_rng_mt19937' and a seed 0 are used as defaults.
    The GSL library provides a number of random number generator types, 
    and one can choose with a constant (({GSL::RNG_xxx})) or a String, as

      * (({GSL::Rng::MT19937})) or (({"gsl_rng_mt19937"})) or (({"mt19937"}))
      * (({GSL::Rng::RANLXS0}))  or (({"gsl_rng_ranlsx0"}))  or (({"ranlxs0"})) 
      * (({GSL::Rng::ZUF}))  or (({"gsl_rng_zuf"}))  or (({"zuf"})) 
      * ...

    See the ((<GSL reference manual|URL:http://www.gnu.org/software/gsl/manual/html_node/Random-number-generator-algorithms.html#Random-number-generator-algorithms>)) for the complete list. 
    The following demonstrates how to use this class,

      require 'gsl'

      r = Rng.alloc(Rng::TAUS, 1)
      r2 = Rng.alloc("ran0", 2)
      p r.name                <- "taus"
      p r.get                 <- get an integer
      p r2.uniform            <- get a float of [0, 1)

    A generator of the type ((|gsl_rng_taus|)) is created with seed 1,
    and ((|gsl_rng_ran0|)) with seed 2. The
    method ((|get|)) returns a random integer. 
    The method ((|uniform|)) returns
    a floating number uniformly distributed in the range [0, 1).

    If the package ((<rngextra|URL:http://www.network-theory.co.uk/download/rngextra/>)) is installed, additional
    two generators are available, 
      * (({GSL::Rng::RNGEXTRA_RNG1})), (({"rngextra_rng1"}))
      * (({GSL::Rng::RNGEXTRA_RNG2})), (({"rngextra_rng2"}))

    They are created as

      r1 = Rng.alloc(Rng::RNGEXTRA_RNG1)
      p r1.name                <- "rng1"

      r2 = Rng.alloc("rngextra_rng2")
      p r2.name                <- "rng2"

--- GSL::Rng.default_seed
    Returns the default seed
--- GSL::Rng.set_default_seed(seed)
--- GSL::Rng.default_seed=(seed)
    Override the default seed by ((|seed|)).

--- GSL::Rng.types_setup()
    Returns an array of all the available generators' names. 

--- GSK::Rng.memcpy(dest, src)
    Copies the random number generator ((|src|))) into the pre-existing generator 
    ((|dest|)), making dest into an exact copy of ((|src|)). 
    The two generators must be of the same type.

--- GSL::Rng#set(s)
    This method initializes the random number generator with a given seed ((|s|)).

== Sampling from a random number generator

--- GSL::Rng#get
    This returns a random integer from the reciever generator.

--- GSL::Rng#uniform
    This method returns a double precision floating point number uniformly 
    distributed in the range [0,1). 

--- GSL::Rng#uniform_pos
    This returns a positive double precision floating point number uniformly 
    distributed in the range (0,1), excluding both 0.0 and 1.0.

--- GSL::Rng#uniform_int(n)
    This method returns a random integer from 0 to n-1 inclusive. 

== Auxiliary random number generator functions

--- GSL::Rng#name
    This method returns a Sting object of the name of the generator.

--- GSL::Rng#max
--- GSL::Rng#min
    These method return the largest/smallest value that the method 
    ((|get|)) can return. 

--- GSL::Rng#clone
--- GSL::Rng#duplicate
    Return a newly created generator which is an exact copy of the generator ((|self|)).

== Random number environment variables 
The library allows you to choose a default generator and seed from the 
environment variables (({GSL_RNG_TYPE})) and (({GSL_RNG_SEED})) 
and the method (({GSL::Rng::env_setup})). 

--- GSL::Rng.env_setup()
    Reads the environment variables (({GSL_RNG_TYPE})) and 
    (({GSL_RNG_SEED})) and uses their values to set the corresponding 
    library variables.

    If you don't specify a generator for (({GSL_RNG_TYPE})) 
    then "mt19937" is used as the default. 
    The initial value of the default seed is zero. 


((<prev|URL:integration.html>))
((<next|URL:qrng.html>))

((<Reference index|URL:ref.html>))
((<top|URL:index.html>))
=end
