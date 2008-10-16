=begin
= Series Acceleration
In Ruby/GSL, series acceleration functions are provided as singleton methods
for the (({GSL::Sum::Levin_u, Levin_utrunc})) classes, and methods of
an object of the ((<GSL::Vector|URL:vector.html>)) class.

== Modules and classes
* GSL
  * Sum (Module)
    * Levin_u (Class)
    * Levin_utrunc (Class)

== Methods

--- GSL::Sum::Levin_u.alloc(size)
--- GSL::Sum::Levin_utrunc.alloc(size)

--- GSL::Sum::Levin_u.accel(v)
    This method takes the terms of a series in vector ((|v|)) and computes 
    the extrapolated limit of the series using a Levin u-transform. This returns
    an array of (({[sum, abserr, sum_plain, terms_used]})),
    where ((|sum|)) is the extrapolated sum, ((|abserr|)) is an estimate of the 
    absolute error, and ((|sum_plain|)) is the actual term-by-term sum.

--- GSL::Sum::Levin_utrunc.accel(v)
    This method takes the terms of a series in vector ((|v|)) and computes 
    the extrapolated limit of the series using a Levin u-transform. This returns
    an array of (({[sum, abserr_trunc, sum_plain, terms_used]})).

--- GSL::Sum::Levin_u#accel(v)
--- GSL::Sum::Levin_u#sum_plain
--- GSL::Sum::Levin_u#terms_used
--- GSL::Sum::Levin_utrunc#accel(v)
--- GSL::Sum::Levin_utrunc#sum_plain
--- GSL::Sum::Levin_utrunc#terms_used

--- GSL::Vector#accel
--- GSL::Vector#accel_sum
--- GSL::Vector#sum_accel
--- GSL::Vector#sum
    These calculate the "extrapolated" sum of the terms contained in a 
    GSL::Vector object, using a Levin u-transform. The returned values is a 
    Ruby array with 4 elements, as [((|sum_accel, err, sum_plain, terms_used|))],
    where ((|sum_accel|)) is the extraplated sum, ((|err|)) is the absolute error,
    ((|sum_plain|)) is the term-by-term sum, and ((|terms_used|)) is the number of
    terms actually used in the calculation.

((<prev|URL:cheb.html>))
((<next|URL:dht.html>))

((<Reference index|URL:ref.html>))
((<top|URL:index.html>))

=end

