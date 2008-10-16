=begin
= GSL::Vector class

Contents:
(1) ((<Class methods|URL:vector.html#1>))
(2) ((<Notes|URL:vector.html#2>))
(3) ((<Methods|URL:vector.html#3>))
    (1) ((<Accessing vector elements|URL:vector.html#3.1>))
    (2) ((<Initializing vector elements|URL:vector.html#3.2>))
    (3) ((<Iterators|URL:vector.html#3.3>))
    (4) ((<IO|URL:vector.html#3.4>))
    (5) ((<Copying vectors|URL:vector.html#3.5>))
    (6) ((<Vector views|URL:vector.html#3.6>))
    (7) ((<Vector operations|URL:vector.html#3.7>))
    (8) ((<Vector operations with size changes|URL:vector.html#3.8>))
    (9) ((<Finding maximum and minimum elements of vectors|URL:vector.html#3.9>))
    (10) ((<Vector properties|URL:vector.html#3.10>))
    (11) ((<Element-wise vector comparison|URL:vector.html#3.11>))
    (12) ((<Histogram|URL:vector.html#3.12>))
    (13) ((<Sorting|URL:vector.html#3.13>))
    (14) ((<BLAS methods|URL:vector.html#3.14>))
    (15) ((<Data type conversions|URL:vector.html#3.15>))
    (16) ((<NArray|URL:vector.html#3.16>))
    (17) ((<GNU graph interface|URL:vector.html#3.17>))

See also ((<GSL::Vector::Complex|URL:vector_complex.html>)).

== Class methods

--- GSL::Vector.alloc(ary) 
--- GSL::Vector.alloc(ary) 
--- GSL::Vector.alloc(range) 
--- GSL::Vector.alloc(size) 
--- GSL::Vector.alloc(elm0, elm1, ....) 
--- GSL::Vector[elm0, elm1, ....]
    Constructors.

    Ex:
       irb(main):002:0> v1 = GSL::Vector.alloc(5)
       => GSL::Vector: [ 0.000e+00 0.000e+00 0.000e+00 0.000e+00 0.000e+00 ]
       irb(main):003:0> v2 = GSL::Vector.alloc(1, 3, 5, 2)
       => GSL::Vector: [ 1.000e+00 3.000e+00 5.000e+00 2.000e+00 ]
       irb(main):004:0> v3 = GSL::Vector[1, 3, 5, 2]
       => GSL::Vector: [ 1.000e+00 3.000e+00 5.000e+00 2.000e+00 ]
       irb(main):005:0> v4 = GSL::Vector.alloc([1, 3, 5, 2])
       => GSL::Vector: [ 1.000e+00 3.000e+00 5.000e+00 2.000e+00 ]
       irb(main):006:0> v5 = GSL::Vector[1..6]
       => GSL::Vector: [ 1.000e+00 2.000e+00 3.000e+00 4.000e+00 5.000e+00 6.000e+00 ]

--- GSL::Vector.calloc(size)
    This method creates a vector object, and initializes all the elements to zero.

--- GSL::Vector.linspace(min, max, n = 10)
    Creates an (({GSL::Vector})) with ((|n|)) linearly spaced elements 
    between ((|min|)) and ((|max|)). If ((|min|)) is greater than ((|max|)), 
    the elements are stored in decreasing order. This mimics the (({linspace}))
    function of ((<GNU Octave|URL:http://www.octave.org/>)).

    Ex: 
        irb(main):002:0> x = GSL::Vector.linspace(0, 10, 5)
        [ 0.000e+00 2.500e+00 5.000e+00 7.500e+00 1.000e+01 ]
        irb(main):003:0> y = GSL::Vector.linspace(10, 0, 5)
        [ 1.000e+01 7.500e+00 5.000e+00 2.500e+00 0.000e+00 ]

--- GSL::Vector.logspace(min, max, n)
    Similar to (({GSL::Vector#linspace})) except that the values are 
    logarithmically spaced from 10^((|min|)) to 10^((|max|)).

    Ex:
        irb(main):007:0* x = GSL::Vector.logspace(1, 3, 5)
        [ 1.000e+01 3.162e+01 1.000e+02 3.162e+02 1.000e+03 ]
        irb(main):008:0> x = GSL::Vector.logspace(3, 1, 5)
        [ 1.000e+03 3.162e+02 1.000e+02 3.162e+01 1.000e+01 ]

--- GSL::Vector.logspace2(min, max, n)
    Similar to (({GSL::Vector#linspace})) except that the values are 
    logarithmically spaced from ((|min|)) to ((|max|)).

    Ex:
        irb(main):010:0* x = GSL::Vector.logspace2(10, 1000, 5)
        [ 1.000e+01 3.162e+01 1.000e+02 3.162e+02 1.000e+03 ]
        irb(main):011:0> x = GSL::Vector.logspace2(1000, 10, 5)
        [ 1.000e+03 3.162e+02 1.000e+02 3.162e+01 1.000e+01 ]

--- GSL::Vector.indgen(n, start=0, step=1)
    This creates a vector of length ((|n|)) with elements from ((|start|))
    with interval ((|step|)) (mimics NArray#indgen).
    
    Ex:
        irb(main):019:0> v = GSL::Vector::Int.indgen(5)
        => GSL::Vector::Int: 
        [ 0 1 2 3 4 ]
        irb(main):020:0> v = GSL::Vector::Int.indgen(5, 3)
        => GSL::Vector::Int: 
        [ 3 4 5 6 7 ]
        irb(main):021:0> v = GSL::Vector::Int.indgen(5, 3, 2)
        => GSL::Vector::Int: 
        [ 3 5 7 9 11 ]

--- GSL::Vector.filescan(filename)
    Reads a formatted ascii file and returns an array of vectors.
    For a data file (({a.dat})) as
      1 5 6 5
      3 5 6 7
      5 6 7 9
    then (({a, b, c, d = Vetor.filescan("a.dat")})) yields
      a = [1, 3, 5]
      b = [5, 5, 6]
      c = [6, 6, 7]
      d = [5, 7, 9]

=== NArray Extension
If an (({NArray})) object is given, a newly allocated vector is created.

Ex:
        na = NArray[1.0, 2, 3, 4, 5]
        p na                <----- NArray.float(5): 
                                   [ 1.0, 2.0, 3.0, 4.0, 5.0]
        v = GSL::Vector.alloc(na)  
        p v                 <----- [ 1 2 3 4 5 ]


See also ((<here|URL:vector.html#3.13>)).

== NOTE:
In Ruby/GSL, vector lendth is limited within the range of Fixnum.
For 32-bit CPU, the maximum of vector length is 2^30 ~ 1e9.

== Methods

=== Accessing vector elements
--- GSL::Vector#get(indices)
--- GSL::Vector#[indices]
    Return elements(s) of the vector ((|self|)). 

--- GSL::Vector#set(i, val)
--- GSL::Vector#[]=
    Set the ((|i|))-th element of the vector ((|self|)) to ((|val|)).

    Ex:
          irb(main):001:0> require("rbgsl")
          => true
          irb(main):002:0> v = GSL::Vector[0..5]
          => GSL::Vector: [ 0.000e+00 1.000e+00 2.000e+00 3.000e+00 4.000e+00 5.000e+00 ]
          irb(main):003:0> v[2]
          => 2.0
          irb(main):004:0> v[1, 3, 4]
          => GSL::Vector: [ 1.000e+00 3.000e+00 4.000e+00 ]
          irb(main):005:0> v[1..3]
          => GSL::Vector::View: [ 1.000e+00 2.000e+00 3.000e+00 ]
          irb(main):006:0> v[3] = 9
          => 9
          irb(main):007:0> v[-1] = 123
          => 123
          irb(main):008:0> v
          => GSL::Vector: [ 0.000e+00 1.000e+00 2.000e+00 9.000e+00 4.000e+00 1.230e+02 ]

=== Initializing vector elements
--- GSL::Vector#set_all(x)
    This method sets all the elements of the vector to the value ((|x|)).

--- GSL::Vector#set_zero
    This method sets all the elements of the vector to zero.

--- GSL::Vector#set_basis!(i)
    This method makes a basis vector by setting all the elements of the vector
    to zero except for the ((|i|))-th element, which is set to one. 
    For a vector (({v})) of size 10, the method
      v.set_basis!(4)
    sets the vector ((|v|)) to a basis vector (({[0, 0, 0, 0, 1, 0, 0, 0, 0, 0]})). 

--- GSL::Vector#set_basis(i)
    This method returns a new basis vector by setting all the elements of the 
    vector to zero except for the i-th element which is set to one. 
    For a vector (({v})) of size 10, the method
      vb = v.set_basis(4)
    creates a new vector ((|vb|)) with elements (({[0, 0, 0, 0, 1, 0, 0, 0, 0, 0]})). 
    The vector ((|v|)) is not changed.

--- GSL::Vector#indgen!(start=0, step=1)
--- GSL::Vector#indgen(start=0, step=1)
    Mimics NArray#indgen!.

=== Iterators
--- GSL::Vector#each
--- GSL::Vector#reverse_each
    An iterator for each of the vector elements, used as

      v.each do |x|    # Show all the elements
        p x
      end

--- GSL::Vector#each_index
--- GSL::Vector#reverse_each_index
    Iterators 

--- GSL::Vector#collect { |item| .. }
    Creates a new vector by collecting the vector elements modified with some
    operations.

    Ex:
      irb(main):003:0> a = GSL::Vector::Int[0..5]
      => GSL::Vector::Int
      [ 0 1 2 3 4 5 ]
      irb(main):004:0> b = a.collect {|v| v*v}
      => GSL::Vector::Int
      [ 0 1 4 9 16 25 ]
      irb(main):005:0> a
      => GSL::Vector::Int
      [ 0 1 2 3 4 5 ]

--- GSL::Vector#collect! { |item| .. }
    Ex:
      irb(main):006:0> a = GSL::Vector::Int[0..5]
      => GSL::Vector::Int
      [ 0 1 2 3 4 5 ]
      irb(main):007:0> a.collect! {|v| v*v}
      => GSL::Vector::Int
      [ 0 1 4 9 16 25 ]
      irb(main):008:0> a
      => GSL::Vector::Int
      [ 0 1 4 9 16 25 ]

=== IO
--- GSL::Vector#print
--- GSL::Vector#fprintf(io, format = "%e")
--- GSL::Vector#fprintf(filename, format = "%e")
--- GSL::Vector#fscanf(io)
--- GSL::Vector#fscanf(filename)
--- GSL::Vector#fwrite(io)
--- GSL::Vector#fwrite(filename)
--- GSL::Vector#fread(io)
--- GSL::Vector#fread(filename)
    Methods for writing or reading the vector. 
    The first argument is an (({IO})) or a (({String})) object.

=== Copying vectors
--- GSL::Vector#clone
--- GSL::Vector#duplicate
    Create a new vector of the same elements.

=== Vector views
The (({GSL::Vector::View})) class is defined to be used as "references" to
vectors. Since the (({Vector::View})) class is a subclass of (({Vector})),
an instance of the (({View})) class created by slicing a (({Vector})) object 
can be used same as the original vector. A
(({View})) object shares the data with the original vector, i.e. any changes
in the elements of the (({View})) object affect to the original vector.

--- GSL::Vector#subvector
--- GSL::Vector#subvector(n)
--- GSL::Vector#subvector(offset, n)
--- GSL::Vector#subvector(offset, stride, n)
    Create a (({Vector::View})) object slicing ((|n|)) elements
    of the vector ((|self|)) from the offset ((|offset|)). If called with one
    argument ((|n|)), ((|offset|)) is set to 0. With no arguments, a view is
    created with the same length of the original vector.

    * Example:
       #!/usr/bin/env ruby
       require("rbgsl")

       v = GSL::Vector[1, 2, 3, 4, 5, 6]
       view = v.subvector(1, 4)
       p view.class         <----- GSL::Vector::View
       view.print           <----- [ 2 3 4 5 ]

       view[2] = 99
       view.print           <----- [ 2 3 99 5 ]
       v.print              <----- [ 1 2 3 99 5 6 ]

--- GSL::Vector#subvector_with_stride(offset, n, stride)
    Return a (({Vector::View})) object of a subvector of another vector ((|self|)) 
    with an additional stride argument. The subvector is formed in the same way 
    as for (({Vector#subvector})) but the new vector view has ((|n|)) elements 
    with a step-size of ((|stride|)) from one element to the next in the original vector. 

--- GSL::Vectir#matrix_view(n1, n2)
    This creates a (({Matrix::View})) object from the vector ((|self|)).
    It enables to use the vector as a ((<Matrix|URL:matrix.html>)) object.

    * Ex:

        irb(main):019:0> v = GSL::Vector::Int.alloc(1..9)
        => GSL::Vector::Int: 
        [ 1 2 3 4 5 6 7 8 9 ]
        irb(main):020:0> m = v.matrix_view(3, 3)
        => GSL::Matrix::Int::View: 
        [ 1 2 3 
          4 5 6 
          7 8 9 ]
        irb(main):021:0> m[1][2] = 99
        => 99
        irb(main):022:0> v
        => GSL::Vector::Int: 
        [ 1 2 3 4 5 99 7 8 9 ]

=== Vector operations

--- GSL::Vector#swap_elements(i, j)
    This method exchanges the i-th and j-th elements of the vector ((|in-place|)). 

--- GSL::Vector#reverse
    Reverses the order of the elements of the vector.

      irb(main):025:0> v = GSL::Vector::Int[1..5]
      => GSL::Vector::Int: 
      [ 1 2 3 4 5 ]
      irb(main):026:0> v.reverse
      => GSL::Vector::Int: 
      [ 5 4 3 2 1 ]

--- GSL::Vector#trans
--- GSL::Vector#transpose
--- GSL::Vector#col
--- GSL::Vector#row
    Transpose the vector from a row vector into a column vector and vice versa.

      irb(main):029:0> v = GSL::Vector::Int[1..5]
      => GSL::Vector::Int: 
      [ 1 2 3 4 5 ]
      irb(main):030:0> v.col
      => GSL::Vector::Int::Col: 
      [ 1 
        2 
        3 
        4 
        5 ]

--- GSL::Vector#add(b)
    Adds the elements of vector ((|b|)) to the elements 
    of the vector ((|self|)). A new vector is created, and the vector
    ((|self|)) is not changed.

--- GSL::Vector#sub(b)
    Subtracts the element of vector ((|b|)) from the elements of ((|self|)).
    A new vector is created, and the vector ((|self|)) is not changed.

--- GSL::Vector#mul(b)
    Multiplies the elements of vector ((|self|)) by the elements of vector ((|b|)).
--- GSL::Vector#div(b)
    Divides the elements of vector ((|self|)) by the elements of vector ((|b|)).

--- GSL::Vector#scale(x)
--- GSL::Vector#scale!(x)
    This method multiplies the elements of vector ((|self|)) by 
    the constant factor ((|x|)).

--- GSL::Vector#add_constant(x)
--- GSL::Vector#add_constant!(x)
    Adds the constant value ((|x|)) to the elements of the vector ((|self|)).

--- GSL::Vector#+(b)
    For ((|b|)),
      * a Number: ---> (({self.add_constanb(b)}))
      * a Vector: ---> (({self.add(b)}))
--- GSL::Vector#-(b)
    For ((|b|)),
      * a Number: ---> (({self.add_constanb(-b)}))
      * a Vector: ---> (({self.sub(b)}))
--- GSL::Vector#/(b)
    For ((|b|)),
      * a Number: ---> (({self.scale(1/b)}))
      * a Vector: ---> (({self.div(b)}))

--- GSL::Vector#*(b)
    Vector multiplication.

    (1) Scale
          irb(main):027:0> v = GSL::Vector[1, 2]
          [ 1 2 ]
          irb(main):028:0> v*2
          [ 2 4 ]                           
    (2) Element-by-element multiplication
          irb(main):018:0> a = GSL::Vector[1, 2]; b = GSL::Vector[3, 4]
          [ 3 4 ]
          irb(main):020:0> a*b
          [ 3 8 ]                             
    (3) Inner product
          irb(main):023:0> a = GSL::Vector[1, 2]; b = GSL::Vector[3, 4]
          [ 3 
            4 ]
          irb(main):024:0> a*b.col
          => 11.0                        
    (4) GSL::Vector::Col*Vector -> GSL::Matrix
          irb(main):025:0> a = GSL::Vector::Col[1, 2]; b = GSL::Vector[3, 4]
          [ 3 4 ]
          irb(main):026:0> a*b
          [ 3 4 
            6 8 ]
    (5) GSL::Matrix*Vector::Col -> GSL::Vector::Col
          irb(main):029:0> a = GSL::Vector[1, 2]; m = GSL::Matrix[[2, 3], [4, 5]]
          [ 2 3 
            4 5 ]
          irb(main):030:0> m*a          <--- Error
          TypeError: Operation with GSL::Vector is not defined (GSL::Vector::Col expected)
                  from (irb):30:in `*'
                  from (irb):30
          irb(main):031:0> m*a.col
          [ 8 
            14 ]

--- GSL::Vector#add!(b)
--- GSL::Vector#sub!(b)
--- GSL::Vector#mul!(b)
--- GSL::Vector#div!(b)
    In-place operations with a vector ((|b|)).

--- GSL::Vector#pow(p)
--- GSL::Vector#**(p)
--- GSL::Vector#pow!(p)
    Element-wise calculation of power p.

    Ex)
       irb(main):001:0> require("rbgsl")
       irb(main):002:0> v = GSL::Vector[1, 2, 3]
       => GSL::Vector
       [ 1.000e+00 2.000e+00 3.000e+00 ]
       irb(main):003:0> v.pow(2)
       => GSL::Vector
       [ 1.000e+00 4.000e+00 9.000e+00 ]
       irb(main):004:0> v**2
       => GSL::Vector
       [ 1.000e+00 4.000e+00 9.000e+00 ]
       irb(main):005:0> v
       => GSL::Vector
       [ 1.000e+00 2.000e+00 3.000e+00 ]
       irb(main):006:0> v.pow!(2)
       => GSL::Vector
       [ 1.000e+00 4.000e+00 9.000e+00 ]
       irb(main):007:0> v
       => GSL::Vector
       [ 1.000e+00 4.000e+00 9.000e+00 ]

--- GSL::Vector#swap_elements(i, j)
    This exchanges the ((|i|))-th and ((|j|))-th elements of the vector ((|self|)) in-place.
--- GSL::Vector#clone
--- GSL::Vector#duplicate
    These create a copy of the vector ((|self|)).
    
--- GSL::Vector.connect(v1, v2, v3, ...)
--- GSL::Vector#connect(v2, v3, ...)
    Creates a new vector by connecting all the elements of the given vectors.

      irb(main):031:0> v1 = GSL::Vector::Int[1, 3]
      => GSL::Vector::Int: 
      [ 1 3 ]
      irb(main):032:0> v2 = GSL::Vector::Int[4, 3, 5]
      => GSL::Vector::Int: 
      [ 4 3 5 ]
      irb(main):033:0> v1.connect(v2)
      => GSL::Vector::Int: 
      [ 1 3 4 3 5 ]

--- GSL::Vector#abs
    Creates a new vector, with elements ((|fabs(x_i)|)).

      irb(main):034:0> v = GSL::Vector::Int[-3, 2, -5, 4]
      => GSL::Vector::Int: 
      [ -3 2 -5 4 ]
      irb(main):035:0> v.abs
      => GSL::Vector::Int: 
      [ 3 2 5 4 ]

--- GSL::Vector#square
--- GSL::Vector#abs2
    Create a new vector, with elements ((|x_i*x_i|)).

      irb(main):036:0> v = GSL::Vector::Int[1..4]
      => GSL::Vector::Int: 
      [ 1 2 3 4 ]
      irb(main):037:0> v.square
      => GSL::Vector::Int: 
      [ 1 4 9 16 ]

--- GSL::Vector#sqrt
    Creates a new vector, with elements (({sqrt}))(((|x_i|))).

--- GSL::Vector#floor
--- GSL::Vector#ceil
--- GSL::Vector#round

    Ex:
      irb(main):002:0> v = GSL::Vector[1.1, 2.7, 3.5, 4.3]
      => GSL::Vector
      [ 1.100e+00 2.700e+00 3.500e+00 4.300e+00 ]
      irb(main):003:0> v.floor
      => GSL::Vector::Int
      [ 1 2 3 4 ]
      irb(main):004:0> v.ceil
      => GSL::Vector::Int
      [ 2 3 4 5 ]
      irb(main):005:0> v.round
      => GSL::Vector::Int
      [ 1 3 4 4 ]

--- GSL::Vector#normalize(nrm = 1.0)
    Creates a new vector of norm ((|nrm|)), by scaling the vector ((|self|)).
--- GSL::Vector#normalize!(nrm = 1.0)
    This normalizes the vector ((|self|)) in-place.

    Ex:
      tcsh> irb
      irb(main):001:0> require("rbgsl")
      => true
      irb(main):002:0> a = GSL::Vector[-1, -2, -3, -4]
      => GSL::Vector: 
      [ -1.000e+00 -2.000e+00 -3.000e+00 -4.000e+00 ]
      irb(main):003:0> b = a.abs
      => GSL::Vector: 
      [ 1.000e+00 2.000e+00 3.000e+00 4.000e+00 ]
      irb(main):004:0> b.sqrt
      => GSL::Vector: 
      [ 1.000e+00 1.414e+00 1.732e+00 2.000e+00 ]
      irb(main):005:0> b.square
      => GSL::Vector: 
      [ 1.000e+00 4.000e+00 9.000e+00 1.600e+01 ]
      irb(main):006:0> c = b.normalize(2)
      => GSL::Vector: 
      [ 2.582e-01 5.164e-01 7.746e-01 1.033e+00 ]
      irb(main):007:0> c.square.sum
      => 2.0

--- GSL::Vector#decimate(n)
    Creates a new vector by averaring every ((|n|)) 
    points of the vector ((|self|)) down to one point.

--- GSL::Vector#diff(k = 1)
    Calculate ((|k|))-th differences of a vector ((|self|)).


--- GSL::Vector#join(sep = " ")
    Converts the vector to a (({String})) by joining all the elements with a
    separator ((|sep|)).

--- GSL::Vector#zip(vec, ...)
--- GSL::Vector.zip(vec, ...)
    Create an (({Array})) of vectors by merging the elements of ((|self|))
    with corresponding elements from each arguments.

    Ex:
      irb(main):001:0> require("rbgsl")
      irb(main):002:0> a = GSL::Vector[4, 5, 6]
      irb(main):003:0> b = GSL::Vector[7, 8, 9]
      irb(main):004:0> GSL::Vector[1, 2, 3].zip(a, b)
      [[ 1.000e+00 4.000e+00 7.000e+00 ], 
       [ 2.000e+00 5.000e+00 8.000e+00 ], 
       [ 3.000e+00 6.000e+00 9.000e+00 ]]
      irb(main):005:0> GSL::Vector[1, 2].zip(a, b)
      [[ 1.000e+00 4.000e+00 7.000e+00 ], 
       [ 2.000e+00 5.000e+00 8.000e+00 ]]
      irb(main):006:0> a.zip(GSL::Vector[1, 2], GSL::Vector[8.0])
      [[ 4.000e+00 1.000e+00 8.000e+00 ],
       [ 5.000e+00 2.000e+00 0.000e+00 ],
       [ 6.000e+00 0.000e+00 0.000e+00 ]]

=== Vector operations with size changes
The methods below change vector length of ((|self|)).
--- GSL::Vector#pop
    Removes the last element from ((|self|)) and returns it, or 
    (({nil})) if empty.
--- GSL::Vector#shift
    Returns the first element from ((|self|)) and removes it. Returns 
    (({nil})) if empty.
--- GSL::Vector#push(x)
--- GSL::Vector#concat(x)
--- GSL::Vector#<<(x)
    Append ((|x|)) ((({Numeric})) or (({GSL::Vector}))) to the end of ((|self|)).
--- GSL::Vector#unshift(x)
    Prepends ((|x|)) to the front of ((|self|)).
--- GSL::Vector#delete_at(i)
    Deletes the element at the specified index ((|i|)), 
    returning that  element, or (({nil})) if the index is out of range.
--- GSL::Vector#delete_if { |x| ... }
    Deletes every element of ((|self|)) for which block evaluates to (({true}))
    and returns a new vector of deleted elements.

=== Finding maximum and minimum elements of vectors

--- GSL::Vector#max
    This method returns the maximum value in the vector.

--- GSL::Vector#min
    This method returns the minimum value in the vector.

--- GSL::Vector#minmax
    This method returns an array of two elements, the minimum and the maximum values
    in the vector ((|self|)).

--- GSL::Vector#max_index
    This method returns the index of the maximum value in the vector. When there are 
    several equal maximum elements then the lowest index is returned. 

--- GSL::Vector#min_index
    This method returns the index of the minimum value in the vector. When there are 
    several equal minimum elements then the lowest index is returned. 

--- GSL::Vector#minmax_index
    This method returns an array of two elements which has the indices 
    of the minimum and the maximum values in the vector ((|self|)).

=== Vector Properties
--- GSL::Vector#size
--- GSL::Vector#len
    Return the vector length.

--- GSL::Vector#sum
    Returns the sum of the vector elements.

--- GSL::Vector#prod
    Returns the product of the vector elements.

--- GSL::Vector#cumsum
    Calculate the cumulative sum of elements of ((|self|)) and returns as a new vector.

--- GSL::Vector#cumprod
    Calculate the cumulative product of elements of ((|self|)) and returns as a new vector.

--- GSL::Vector#isnull
    Returns 1 if all the elements of the vector ((|self|)) 
    are zero, and 0 otherwise.
--- GSL::Vector#isnull?
    Return (({true})) if all the elements of the vector ((|self|)) 
    are zero, and (({false})) otherwise.
--- GSL::Vector#ispos
--- GSL::Vector#ispos?
--- GSL::Vector#isneg
--- GSL::Vector#isneg?
    (GSL-1.9 or later) Return 1 (true) if all the elements of the vector ((|self|)) are zero, strictly positive, strictly negative respectively, and 0 (false) otherwise.

--- GSL::Vector#isnonneg
--- GSL::Vector#isnonneg?
    (GSL-1.10 or later) Return 1 (true) if all the elements of the vector ((|self|)) are non-negative , and 0 (false) otherwise. 

--- GSL::Vector#all?
    Returns (({true})) if all the vector elements are non-zero, and (({false})) 
    otherwise. If a block is given, the method returns (({true})) if the
    tests are true for all the elements.
--- GSL::Vector#any?
    Returns (({true})) if any the vector elements are non-zero, and (({false})) 
    otherwise. If a block is given, the method returns (({true})) if the
    tests are true for any of the elements.
--- GSL::Vector#none?
    Returns (({true})) if all the elements of the vector ((|self|)) 
    are zero, and (({false})) otherwise (just as (({GSL::Vector#isnull?}))).
    If a block is given, the method returns (({true})) if the
    tests are false for all the elements.

    Ex:
      irb(main):009:0> a = GSL::Vector[1, 2, 3]
      irb(main):010:0> b = GSL::Vector[1, 2, 0]
      irb(main):011:0> c = GSL::Vector[0, 0, 0]
      irb(main):012:0> a.all?
      => true
      irb(main):013:0> b.all?
      => false
      irb(main):014:0> b.any?
      => true
      irb(main):015:0> c.any?
      => false
      irb(main):016:0> a.none?
      => false
      irb(main):017:0> c.none?
      => true

--- GSL::Vector#all
--- GSL::Vector#any
--- GSL::Vector#none
    Returns 1 or 0.

--- GSL::Vector#equal?(other, eps = 1e-10)
--- GSL::Vector#==(other, eps = 1e-10)
    Returns (({true})) if the vectors have same size and elements
    equal to absolute accurary ((|eps|)) for all the indices, 
    and (({false})) otherwise.

=== Element-wise vector comparison
--- GSL::Vector#eq(other)
--- GSL::Vector#ne(other)
--- GSL::Vector#gt(other)
--- GSL::Vector#ge(other)
--- GSL::Vector#lt(other)
--- GSL::Vector#le(other)
    Return a (({Block::Byte})) object with elements 0/1 by comparing the two vectors
    ((|self|)) and ((|other|)). Note that the values returned are 0/1, 
    not (({true/false})), thus all of the elements are "true" in Ruby.

    Ex:
       irb(main):003:0> a = GSL::Vector[1, 2, 3]
       irb(main):004:0> b = GSL::Vector[1, 2, 5]
       irb(main):005:0> a.eq(b)
       [ 1 1 0 ]
       irb(main):006:0> a.ne(b)
       [ 0 0 1 ]
       irb(main):007:0> a.gt(b)
       [ 0 0 0 ]
       irb(main):008:0> a.ge(b)
       [ 1 1 0 ]
       irb(main):009:0> a.eq(3)
       [ 0 0 1 ]
       irb(main):010:0> a.ne(2)
       [ 1 0 1 ]
       irb(main):011:0> a.ge(2)
       [ 0 1 1 ]

--- GSL::Vector#and(other)
--- GSL::Vector#or(other)
--- GSL::Vector#xor(other)    
--- GSL::Vector#not

    Ex:
      irb(main):033:0> a = GSL::Vector[1, 0, 3, 0]
      irb(main):034:0> b = GSL::Vector[3, 4, 0, 0]
      irb(main):035:0> a.and(b)
      [ 1 0 0 0 ]
      irb(main):036:0> a.or(b)
      [ 1 1 1 0 ]
      irb(main):037:0> a.xor(b)
      [ 0 1 1 0 ]
      irb(main):038:0> a.not
      [ 0 1 0 1 ]
      irb(main):039:0> b.not
      [ 0 0 1 1 ]

--- GSL::Vector#where
--- GSL::Vector#where { |elm| ... }
    Returns the vector indices where the tests are true. If all the test failed
    ((|nil|)) is returned.

    Ex:
      irb(main):003:0> v = GSL::Vector::Int[0, 3, 0, -2, 3, 5, 0, 3]
      irb(main):004:0> v.where
      [ 1 3 4 5 7 ]                   # where elements are non-zero
      irb(main):007:0> v.where { |elm| elm == -2 }
      [ 3 ]
      irb(main):008:0> a = GSL::Vector[0, 0, 0]
      irb(main):009:0> a.where
      => nil

=== Histogram
--- GSL::Vector#histogram(n)
--- GSL::Vector#histogram(ranges)
--- GSL::Vector#histogram(n, min, max)
--- GSL::Vector#histogram(n, [min, max])
    Creates a histogram filling the vector ((|self|)).

    Example:
       irb(main):003:0> r = GSL::Rng.alloc           # Random number generator
       => #<GSL::Rng:0x6d8594>
       irb(main):004:0> v = r.gaussian(1, 1000)    # Generate 1000 Gaussian random numbers
       => GSL::Vector
       [ 1.339e-01 -8.810e-02 1.674e+00 7.336e-01 9.975e-01 -1.278e+00 -2.397e+00 ... ]
       irb(main):005:0> h = v.histogram(50, [-4, 4])  # Creates a histogram of size 50, range [-4, 4)
       => #<GSL::Histogram:0x6d28b0>
       irb(main):006:0> h.graph("-T X -C -g 3")    # Show the histogram
       => true

    This is equivalent to
       h = Histogram.alloc(50, [-4, 4])
       h.increment(v)

=== Sorting

--- GSL::Vector#sort
--- GSL::Vector#sort!
    These methods sort the vector ((|self|)) in ascending numerical order.

--- GSL::Vector#sort_index
    This method indirectly sorts the elements of the vector ((|self|)) into 
    ascending order, and returns the resulting permutation. 
    The elements of permutation give the index of the vector element which 
    would have been stored in that position if the vector had been sorted in place. 
    The first element of permutation gives the index of the least element in the
    vector, and the last element of permutation gives the index of the greatest 
    vector element. The vector ((|self|)) is not changed.

--- GSL::Vector#sort_smallest(n)
--- GSL::Vector#sort_largest(n)
--- GSL::Vector#sort_smallest_index(n)
--- GSL::Vector#sort_largest_index(n)

    Ex:
        irb(main):005:0> v = GSL::Vector::Int[8, 2, 3, 7, 9, 1, 4]
        => GSL::Vector::Int: 
        [ 8 2 3 7 9 1 4 ]
        irb(main):006:0> v.sort
        => GSL::Vector::Int: 
        [ 1 2 3 4 7 8 9 ]
        irb(main):007:0> v.sort_index
        => GSL::Permutation: 
        [ 5 1 2 6 3 0 4 ]
        irb(main):008:0> v.sort_largest(3)
        => GSL::Vector::Int: 
        [ 9 8 7 ]
        irb(main):009:0> v.sort_smallest(3)
        => GSL::Vector::Int: 
        [ 1 2 3 ]

=== BLAS Methods
--- GSL::Vector#nrm2
--- GSL::Vector#dnrm2
    Compute the Euclidean norm ||x||_2 = sqrt {sum x_i^2} of the vector.

--- GSL::Vector#asum
--- GSL::Vector#dasum
    Compute the absolute sum \sum |x_i| of the elements of the vector.

=== Data type conversions
--- GSL::Vector#to_a 
    This method converts the vector into a Ruby array. A Ruby array also can be
    converted into a GSL::Vector object with the (({to_gv})) method. For example,

      v = GSL::Vector.alloc([1, 2, 3, 4, 5])
      a = v.to_a   -> GSL::Vector to an array
      p a          -> [1.0, 2.0, 3.0, 4.0, 5.0]
      a[2] = 12.0
      v2 = a.to_gv  -> a new GSL::Vector object
      v2.print     -> 1.0000e+00 2.0000e+00 1.2000e+01 4.0000e+00 5.0000e+00

--- GSL::Vector#to_m(nrow, ncol)
    Creates a (({GSL::Matrix})) object of ((|nrow|)) rows and ((|ncol|)) columns.

      irb(main):011:0> v = GSL::Vector::Int[1..5]
      => GSL::Vector::Int: 
      [ 1 2 3 4 5 ]
      irb(main):012:0> v.to_m(2, 3)
      => GSL::Matrix::Int: 
      [ 1 2 3 
        4 5 0 ]
      irb(main):013:0> v.to_m(2, 2)
      => GSL::Matrix::Int: 
      [ 1 2 
        3 4 ]
      irb(main):014:0> v.to_m(3, 2)
      => GSL::Matrix::Int: 
      [ 1 2 
        3 4 
        5 0 ]

--- GSL::Vector#to_m_diagonal
    Converts the vector into a diagonal matrix. 
    See also ((<GSL::Matrix.diagonal(v)|URL:matrix.html>)).

       irb(main):012:0> v = GSL::Vector[1..4].to_i
       => GSL::Vector::Int: 
       [ 1 2 3 4 ]
       irb(main):013:0> v.to_m_diagonal
       => GSL::Matrix::Int: 
       [ 1 0 0 0 
         0 2 0 0 
         0 0 3 0 
         0 0 0 4 ]

--- GSL::Vector#to_m_circulant
    Creates a circulant matrix.

      irb(main):002:0> v = GSL::Vector::Int[1..5]
      => GSL::Vector::Int: 
      [ 1 2 3 4 5 ]
      irb(main):003:0> v.to_m_circulant
      => GSL::Matrix::Int: 
      [ 5 1 2 3 4   
        4 5 1 2 3   
        3 4 5 1 2   
        2 3 4 5 1   
        1 2 3 4 5 ]

--- GSL::Vector#to_complex
--- GSL::Vector#to_complex2
    Example:

      irb(main):002:0> v = GSL::Vector[1..4]
      => GSL::Vector
      [ 1.000e+00 2.000e+00 3.000e+00 4.000e+00 ]
      irb(main):003:0> v.to_complex
      [ [1.000e+00 0.000e+00] [2.000e+00 0.000e+00] [3.000e+00 0.000e+00] [4.000e+00 0.000e+00] ]
      => #<GSL::Vector::Complex:0x6d7d24>
      irb(main):004:0> v.to_complex2
      [ [1.000e+00 2.000e+00] [3.000e+00 4.000e+00] ]
      => #<GSL::Vector::Complex:0x6d6424>

--- GSL::Vector#to_tensor(rank, dimension)

=== (({GSL::Vector <---> NArray}))

--- GSL::Vector#to_na
    Converts a vector ((|self|)) into an (({NArray})) object. 
    The data are copied to newly allocated memory.

--- GSL::Vector#to_na2
--- GSL::Vector#to_na_ref
    Create an (({NArray})) reference of the vector ((|self|)).

    Example:
      irb(main):020:0> v = GSL::Vector::Int[1, 2, 3, 4]
      => GSL::Vector::Int
      [ 1 2 3 4 ]
      irb(main):021:0> na = v.to_na
      => NArray.int(4): 
      [ 1, 2, 3, 4 ]
      irb(main):022:0> na2 = v.to_na2
      => NArray(ref).int(4): 
      [ 1, 2, 3, 4 ]
      irb(main):023:0> na[1] = 99
      => 99
      irb(main):024:0> v              # na and v are independent
      => GSL::Vector::Int
      [ 1 2 3 4 ]
      irb(main):025:0> na2[1] = 99    # na2 points to the data of v
      => 99
      irb(main):026:0> v
      => GSL::Vector::Int
      [ 1 99 3 4 ]

--- NArray#to_gv
--- NArray#to_gslv
    Create (({GSL::Vector})) object from the (({NArray})) object ((|self|)).

--- NArray#to_gv_view
--- NArray#to_gv2
--- NArray#to_gslv_view
    A (({GSL::Vector::View})) object is created from the NArray object ((|self|)). 
    This method does not allocate memory for the data: the data of ((|self|)) 
    are not copied, but shared with the (({View})) object created, thus 
    any modifications to the (({View})) object affect on the original NArray 
    object. In other words, the (({View})) object can be used as a ((|reference|)) 
    to the NArray object. 

    Ex:
       tcsh> irb
       irb(main):001:0> require("rbgsl")
       => true
       irb(main):002:0> na = NArray[1.0, 2, 3, 4, 5]    
       => NArray.float(5): 
       [ 1.0, 2.0, 3.0, 4.0, 5.0 ]
       irb(main):003:0> vv = na.to_gv_view   # Create a view sharing the memory
       => GSL::Vector::View 
       [ 1.000e+00 2.000e+00 3.000e+00 4.000e+00 5.000e+00 ]
       irb(main):004:0> vv[3] = 9
       => 9
       irb(main):005:0> na
       => NArray.float(5): 
       [ 1.0, 2.0, 3.0, 9.0, 5.0 ]           # The data are changed
       irb(main):006:0> v = na.to_gv         # A vector with newly allocated memory
       => GSL::Vector 
       [ 1.000e+00 2.000e+00 3.000e+00 9.000e+00 5.000e+00 ]
       irb(main):007:0> v[1] = 123
       => 123
       irb(main):008:0> v
       => GSL::Vector 
       [ 1.000e+00 1.230e+02 3.000e+00 9.000e+00 5.000e+00 ]
       irb(main):009:0> na                   
       => NArray.float(5): 
       [ 1.0, 2.0, 3.0, 9.0, 5.0 ]           # v and na are independent 
       irb(main):010:0> na = NArray[1.0, 2, 3, 4, 5, 6]
       => NArray.float(6): 
       [ 1.0, 2.0, 3.0, 4.0, 5.0, 6.0 ]
       irb(main):011:0> m = na.to_gv_view.matrix_view(2, 3)
       => GSL::Matrix::View
       [  1.000e+00  2.000e+00  3.000e+00 
          4.000e+00  5.000e+00  6.000e+00 ]
       irb(main):012:0> m[1][2] = 9
       => 9
       irb(main):013:0> na
       => NArray.float(6): 
       [ 1.0, 2.0, 3.0, 4.0, 5.0, 9.0 ]

=== Graphics
--- GSL::Vector.graph(y)
--- GSL::Vector.graph(y, options)
--- GSL::Vector.graph(x, y)
--- GSL::Vector.graph(x, y, options)
--- GSL::Vector#graph(options)
--- GSL::Vector#graph(x, options)
    These methods use the GNU plotutils (({graph})) application to plot
    a vector ((|self|)). The options of (({graph})) as "-T X -C" can be given by a String.

    Example:
      irb(main):008:0> x = GSL::Vector.linspace(0, 2.0*M_PI, 20)
      irb(main):009:0> c = GSL::Sf::cos(x)
      irb(main):010:0> s = GSL::Sf::sin(x)
      irb(main):011:0> GSL::Vector.graph(x, c, s, "-T X -C -L 'cos(x), sin(x)'")

((<prev|URL:sf.html>))
((<next|URL:matrix.html>))

((<Reference index|URL:ref.html>))
((<top|URL:index.html>))
    
=end
