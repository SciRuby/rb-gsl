=begin
= Sorting
Contents:
(1) ((<Heapsort of vectors|URL:sort.html#1>))
(2) ((<Sorting vectors|URL:sort.html#2>))
(3) ((<Selecting the k smallest or largest elements|URL:sort.html#3>))

== Heapsort

--- GSL::Vector#heapsort
--- GSL::Vector::Complex#heapsort
--- GSL.heapsort(v)
    These method sort the elements of the vector ((|self|)) 
    using the comparison function given by a block, and return the result
    as a new vector object. The vector ((|self|)) is not changed.

    Example: Sorting a complex vector in ascending order by magnitudes.

        p v.heapsort { |a, b|
          a.abs <=> b.abs
        }

--- GSL::Vector#heapsort!
--- GSL::Vector::Complex#heapsort!
--- GSL.heapsort!(v)
    These method sort the elements of the vector ((|self|)) in-place.

--- GSL::Vector#heapsort_index
--- GSL::Vector::Complex#heapsort_index
--- GSL.heapsort_index(v)
    These method indirectly sort the elements of the vector ((|self|)) 
    using the comparison 
    function given by a block, and return the result
    as a permutation object. The vector itself is not changed.
    
== Sorting vectors
--- GSL::Vector#sort!
    This method sorts the elements of the vector ((|self|)) into 
    ascending numerical order. The vector itself is changed.

--- GSL::Vector#sort
    This returns a new vector whose elements are sorted into ascending 
    numerical order. The vector ((|self|)) is not changed.

--- GSL::Vector#sort_index
    This method indirectly sorts the elements of the vector ((|self|)) 
    into ascending order, 
    and returns the result as a (({GSL::Permutation})) object. 
    The elements of the returned permutation give the index of the vector 
    element which would
    have been stored in that position if the vector had been sorted in place. 
    The first element of the permutation gives the index of the least element 
    in ((|self|)),  and the last element of the permutation gives the index 
    of the greatest element in 
    ((|self|)). The vector ((|self|)) is not changed.

== Selecting the k smallest or largest elements
--- GSL::Vector#sort_smallest(k)
--- GSL::Vector#sort_largest(k)
    These functions return a new vector of the ((|k|)) smallest or 
    largest elements  of the vector ((|self|)). 
    The argument ((|k|)) must be less than or equal to the length 
    of the vector ((|self|)). 

--- GSL::Vector#sort_smallest_index(k)
--- GSL::Vector#sort_largest_index(k)
    These functions return a new (({GSL::Permutation})) object of the indices of the 
    ((|k|)) smallest or largest elements of the vector ((|self|)). 
    ((|k|)) must be less than or equal to the length of the vector.


((<prev|URL:combi.html>))
((<next|URL:blas.html>))

((<Reference index|URL:ref.html>))
((<top|URL:index.html>))

=end
