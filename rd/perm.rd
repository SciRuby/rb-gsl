= Permutations
Contents:
(1) ((<Permuation allocations|URL:perm.html#1>))
(2) ((<Methods|URL:perm.html#2>))
    (1) ((<Accessing permutation elements|URL:perm.html#2.1>))
    (2) ((<Permuation properties|URL:perm.html#2.2>))
    (3) ((<Permuation functions|URL:perm.html#2.3>))
    (4) ((<Reading and writing permutations|URL:perm.html#2.4>))
    (5) ((<Permutations in cyclic form|URL:perm.html#2.5>))
(3) ((<Applying Permutations|URL:perm.html#3>))

== Permuation allocations
--- GSL::Permutation.alloc(n)
    These functions create a new permutation of size ((|n|)). 
    The permutation is not initialized and its elements are undefined. 
    Use (({GSL::Permutation.calloc})) if you want to create a permutation 
    which is initialized to the identity. 

--- GSL::Permutation.calloc(n)
    This creates a new permutation of size ((|n|)) and initializes it to the identity. 

== Methods
--- GSL::Permutation#init()
    This initializes the permutation to the identity, i.e. (0,1,2,...,n-1). 

--- GSL::Permutation.memcpy(dest, src)
    This method copies the elements of the permutation ((|src|)) 
    into the permutation ((|dest|)). The two permutations must have the same size.

--- GSL::Permutation#clone
    This creates a new permutation with the same elements of ((|self|)).

=== Accessing permutation elements

--- GSL::Permutation#get(i)
    Returns the value of the ((|i|))-th element of the permutation. 

--- GSL::Permutation#swap(i, j)
    This exchanges the ((|i|))-th and ((|j|))-th elements of the permutation.

=== Permutation properties
--- GSL::Permutation#size
    Returns the size of the permutation.
--- GSL::Permutation#valid
    This checks that the permutation ((|self|)) is valid. 
    The n elements should contain each of the numbers 0 .. n-1 once and only once.

--- GSL::Permutation#valid?
    This returns true if the permutation ((|self|)) is valid, and false otherwise.

=== Permutation functions

--- GSL::Permutation#reverse
    This reverses the elements of the permutation ((|self|)).
--- GSL::Permutation#inverse
    This computes the inverse of the permutation ((|self|)), and returns
    as a new permutation.

--- GSL::Permutation#next
    This method advances the permutation ((|self|)) to the next permutation in 
    lexicographic order and returns (({GSL::SUCCESS})). If no further permutations 
    are available it returns (({GSL::FAILURE})) and leaves ((|self|)) unmodified. 
    Starting with the identity permutation and repeatedly applying this function 
    will iterate through all possible permutations of a given order.
--- GSL::Permutation#prev
    This method steps backwards from the permutation ((|self|)) to the previous 
    permutation in lexicographic order, returning (({GSL_SUCCESS})). 
    If no previous permutation is available it returns (({GSL_FAILURE})) 
    and leaves ((|self|)) unmodified.

=== Reading and writing permutations
--- GSL::Permutation#fwrite(io)
--- GSL::Permutation#fwrite(filename)
--- GSL::Permutation#fread(io)
--- GSL::Permutation#fread(filename)
--- GSL::Permutation#fprintf(io, format = "%u\n")
--- GSL::Permutation#fprintf(filename, format = "%u\n")
--- GSL::Permutation#fscanf(io)
--- GSL::Permutation#fscanf(filename)

=== Permutations in cyclic Form
A permutation can be represented in both ((|linear|)) and 
((|cyclic|)) notations. The functions described in this section convert 
between the two forms. The linear notation is an index mapping, and has 
already been described above. The cyclic notation expresses a 
permutation as a series of circular rearrangements of groups 
of elements, or ((|cycles|)). 

For example, under the cycle (1 2 3), 1 is replaced by 2, 2 is replaced 
by 3 and 3 is replaced by 1 in a circular fashion. Cycles of different 
sets of elements can be combined independently, for example (1 2 3) (4 5) 
combines the cycle (1 2 3) with the cycle (4 5), which is an exchange of 
elements 4 and 5. A cycle of length one represents an element which is 
unchanged by the permutation and is referred to as a ((|singleton|)). 

It can be shown that every permutation can be decomposed into combinations 
of cycles. The decomposition is not unique, but can always be rearranged 
into a standard ((|canonical form|)) by a reordering of elements. 
The library uses the canonical form defined in Knuth's 
((|Art of Computer Programming|)) (Vol 1, 3rd Ed, 1997) Section 1.3.3, p.178. 

The procedure for obtaining the canonical form given by Knuth is, 


(1) Write all singleton cycles explicitly 
(2) Within each cycle, put the smallest number first 
(3) Order the cycles in decreasing order of the first number in the cycle. 

For example, the linear representation (2 4 3 0 1) is represented as 
(1 4) (0 2 3) in canonical form. The permutation corresponds to an 
exchange of elements 1 and 4, and rotation of elements 0, 2 and 3. 

The important property of the canonical form is that it can be reconstructed 
from the contents of each cycle without the brackets. In addition, by removing 
the brackets it can be considered as a linear representation of a different 
permutation. In the example given above the permutation (2 4 3 0 1) would 
become (1 4 0 2 3). This mapping has many applications in the theory of 
permutations. 

--- GSL::Permutation#linear_to_canonical
--- GSL::Permutation#to_canonical
    Computes the canonical form of the permutation ((|self|)) and 
    returns it as a new (({GSL::Permutation})).

--- GSL::Permutation#canonical_to_linear
--- GSL::Permutation#to_linear
    Converts a permutation ((|self|)) in canonical form back into linear 
    form and returns it as a new (({GSL::Permutation})).


--- GSL::Permutation#inversions
    Counts the number of inversions in the permutation ((|self|)). 
    An inversion is any pair of elements that are not in order. 
    For example, the permutation 2031 has three inversions, corresponding
    to the pairs (2,0) (2,1) and (3,1). 
    The identity permutation has no inversions. 

--- GSL::Permutation#linear_cycles
    Counts the number of cycles in the permutation ((|self|)), 
    given in linear form. 

--- GSL::Permutation#canonical_cycles
    Counts the number of cycles in the permutation ((|self|)), 
    given in canonical form. 

== Applying Permutations
--- GSL::Permutation::permute(v)
    Applies the permutation ((|self|)) to the elements of the vector ((|v|)), 
    considered as a row-vector acted on by a permutation matrix from the 
    right, v' = v P. The j-th column of the permutation matrix P is 
    given by the p_j-th column of the identity matrix. 
    The permutation ((|self|)) and the vector ((|v|)) must have the same length.
--- GSL::Permutation::permute_inverse(v)
    Applies the inverse of the permutation ((|self|)) to the elements of 
    the vector ((|v|)), considered as a row-vector acted on by an inverse 
    permutation matrix from the right, v' = v P^T. 
    Note that for permutation matrices the inverse is the same as the 
    transpose. The j-th column of the permutation matrix P is given by 
    the p_j-th column of the identity matrix. 
    The permutation ((|self|)) and the vector ((|v|)) must have the same length.
--- GSL::Permutation.mul(pa, pb)
    Combines the two permutations ((|pa|)) and ((|pb|)) into a single 
    permutation ((|p|)) and returns it. 
    The permutation ((|p|)) is equivalent to applying ((|pb|)) first 
    and then ((|pa|)). 

right
((<prev|URL:matrix.html>))
((<next|URL:combi.html>))

((<Reference index|URL:ref.html>))
((<top|URL:index.html>))

