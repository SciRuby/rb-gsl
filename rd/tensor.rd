=begin
= Tensor manipulations
The tensor library is developed by J. Burguet and distributed 
as an add-on package of GSL. See ((<here|URL:http://sources.redhat.com/ml/gsl-discuss/2004-q4/msg00053.html>)) and ((<here|URL:http://sources.redhat.com/ml/gsl-discuss/2004-q4/msg00055.html>)).

== Class methods
--- GSL::Tensor.new(rank, dimention)
--- GSL::Tensor.alloc(rank, dimention)
--- GSL::Tensor[rank, dimention]
    Create a tensor of rank ((|rank|)) and dimension ((|dimention|)).

--- GSL::Tensor.calloc(rank, dimention)
    Creates a tensor of rank ((|rank|)) and dimension ((|dimention|)),
    and initializes all the elements to zero.
 
--- GSL::Tensor.copy(tensor)
    Create a tensor copying the existing tensor ((|tensor|)).

--- GSL::Tensor.memcpy(dest, src)
    Copies the tensor ((|src|)) to another ((|dest|)). The two
    tensors must have the same shape.

--- GSL::Tensor.swap(a, b)
    Exchanges the elements of the tensor ((|a|)) and ((|b|)).

== Instance methods
=== Accessing tensor elements
--- GSL::Tensor#set_zero
    Sets all the element of the tensor ((|self|)) to zero.
--- GSL::Tensor#set_all(x)
    Sets all the element of the tensor ((|self|)) to ((|x|)).
--- GSL::Tensor#set(indices, x)
--- GSL::Tensor#[indices]=x
    Sets the element of the given indices to ((|x|)).

--- GSL::Tensor#get(indices)
--- GSL::Tensor#[indices]
    Returns the tensor element. If the number of indices given is smaller than the
    rank of the tensor, the method GSL::Tensor#subtensor is called.

    Ex:
        irb(main):002:0> t = Tensor.new(2, 3)
        => #<GSL::Tensor:0x762ae8>
        irb(main):003:0> t.set(1, 2, 2, 123)
        => #<GSL::Tensor:0x762ae8>
        irb(main):004:0> t.get(1, 2, 2)
        => 123.0
        irb(main):005:0> t[0, 0, 2] = 456
        => 456
        irb(main):006:0> t[0, 0, 2]
        => 456.0

--- GSL::Tensor#subtensor(indices)
--- GSL::Tensor#[indices]
    Return a subtensor.
 
    Ex:
        irb(main):001:0> require("gsl")
        => true
        irb(main):002:0> t = Vector[1..125].to_tensor(3, 5)
        => GSL::Tensor: 
        [ 1.000e+00 2.000e+00 3.000e+00 4.000e+00 5.000e+00 6.000e+00 7.000e+00 ... ]
        irb(main):003:0> t[0]
        => GSL::Tensor::View: 
        [ 1.000e+00 2.000e+00 3.000e+00 4.000e+00 5.000e+00 
          6.000e+00 7.000e+00 8.000e+00 9.000e+00 1.000e+01 
          1.100e+01 1.200e+01 1.300e+01 1.400e+01 1.500e+01 
          1.600e+01 1.700e+01 1.800e+01 1.900e+01 2.000e+01 
          2.100e+01 2.200e+01 2.300e+01 2.400e+01 2.500e+01 ]
        irb(main):004:0> t[0,2]
        => GSL::Tensor::View: 
        [ 1.100e+01 1.200e+01 1.300e+01 1.400e+01 1.500e+01 ]
        irb(main):005:0> t[3,1]
        => GSL::Tensor::View: 
        [ 8.100e+01 8.200e+01 8.300e+01 8.400e+01 8.500e+01 ]
        irb(main):006:0> t[1][2]
        => GSL::Tensor::View: 
        [ 3.600e+01 3.700e+01 3.800e+01 3.900e+01 4.000e+01 ]

--- GSL::Tensor#swap_indices(i, j)
--- GSL::Tensor#data
    Returns the data as (({GSL::Vector::View})).
--- GSL::Tensor#to_v
    Creates a new vector from the tensor.

--- GSL::Tensor#to_vector
    Converts the tensor of rank 1 into a (({GSL::Vector::View})) object.
--- GSL::Tensor#to_matrix
    Converts the tensor of rank 2 into a (({GSL::Matrix::View})) object.
    
=== IO
--- GSL::Tensor#fwrite(io)
--- GSL::Tensor#fwrite(filename)
--- GSL::Tensor#fread(io)
--- GSL::Tensor#fread(filename)
--- GSL::Tensor#fprintf(io, format="%g")
--- GSL::Tensor#fprintf(filename, format="%g")
--- GSL::Tensor#fscanf(io)
--- GSL::Tensor#fscanf(filename)

=== Max, min
--- GSL::Tensor#max
--- GSL::Tensor#min
--- GSL::Tensor#minmax
--- GSL::Tensor#max_index
--- GSL::Tensor#min_index
--- GSL::Tensor#minmax_index

=== Tensor operations
--- GSL::Tensor#add(b)
--- GSL::Tensor#+(b)
    Creates a new tensor adding two tensors ((|self|)) and ((|b|)).
--- GSL::Tensor#add!(b)
    Adds the element of tensor ((|b|)) to the elements of ((|self|)) , ((|in-place|)).
--- GSL::Tensor#sub(b)
--- GSL::Tensor#+(b)
    Creates a new tensor subtracting the tensors ((|b|)) from ((|self|)).
--- GSL::Tensor#sub!(b)
    Subtracts the element of tensor ((|b|)) from the elements of ((|self|)) , ((|in-place|)).
--- GSL::Tensor#mul_elements(b)
    This calculate element-by-element multiplication of ((|self|)) and ((|b|)), 
    and returns a new tensor.
--- GSL::Tensor#mul_elements!(b)
    Multiplies the elements of tensor ((|self|)) to the elements of ((|b|)) , ((|in-place|)).
--- GSL::Tensor#div_elements(b)
--- GSL::Tensor#/(b)
    This calculate element-by-element division of ((|self|)) and ((|b|)), 
    and returns a new tensor.
    Multiplies the elements of tensor ((|b|)) to the elements of ((|self|)) , ((|in-place|)).
--- GSL::Tensor#div_elements!(b)
    Divides the elements of tensor ((|self|)) to the elements of ((|b|)) , ((|in-place|)).
--- GSL::Tensor#add_constant(x)
    Creates a new tensor adding the constant ((|x|)) to the tensor ((|self|)).
--- GSL::Tensor#add_constant!(x)
    Adds the constant ((|x|)) to the elements of tensor ((|self|)) , ((|in-place|)).
--- GSL::Tensor#scale(x)
    Creates a new tensor scaling the tensor ((|self|)) by the constant ((|x|)).
--- GSL::Tensor#scale!(x)
    Multiplies the constant ((|x|)) to the elements of tensor ((|self|)) , ((|in-place|)).
--- GSL::Tensor#add_diagonal(x)
    Creates a new tensor adding the constant ((|x|)) to the diagonal elements
    of the tensor ((|self|)).
--- GSL::Tensor#add_diagonal!(x)
    Adds the constant ((|x|)) to the diagonal elements of tensor ((|self|)) , ((|in-place|)).
--- GSL::Tensor#product(b)
--- GSL::Tensor#*(b)
    Calculate tensorian product of ((|self|)) and ((|b|)).
--- GSL::Tensor#contract(i, j)

--- GSL::Tensor#equal?(b, eps = 1e-10)
--- GSL::Tensor#==(b)
    Returns (({true})) if the tensors have same size and elements
    equal to absolute accurary ((|eps|)) for all the indices, 
    and (({false})) otherwise.

=== Tensor properties
--- GSL::Tensor#isnull
    Returns 1 if all the elements of the tensor are zero, and 0 otherwise.
--- GSL::Tensor#isnull?
    Returns (({true})) if all the elements of the tensor are zero, and (({false})) otherwise.

--- GSL::Tensor#rank
    Returns the rank
--- GSL::Tensor#dimension
    Returns the dimension
--- GSL::Tensor#size
    Returns the size

((<prev|URL:rngextra.html>))
((<next|URL:narray.html>))

((<Reference index|URL:ref.html>))
((<top|URL:index.html>))

=end
