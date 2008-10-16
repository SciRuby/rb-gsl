=begin
= Matrices
Contents:
(1) ((<Class methods|URL:matrix.html#1>))
(2) ((<Instance methods|URL:matrix.html#2>))
    (1) ((<Accessing matrix elements|URL:matrix.html#2.1>))
    (2) ((<Initializing matrix elements|URL:matrix.html#2.2>))
    (3) ((<IO|URL:matrix.html#2.3>))
    (4) ((<Matrix views|URL:matrix.html#2.4>))
    (5) ((<Creating row and column views|URL:matrix.html#2.5>))
    (6) ((<Iterators|URL:matrix.html#2.6>))
    (7) ((<Copying matrices|URL:matrix.html#2.7>))
    (8) ((<Copying rows and columns|URL:matrix.html#2.8>))
    (9) ((<Exchanging rows and columns|URL:matrix.html#2.9>))
    (10) ((<Matrix operations|URL:matrix.html#2.10>))
    (11) ((<Finding maximum and minimum elements of matrices|URL:matrix.html#2.11>))
    (12) ((<Matrix properties|URL:matrix.html#2.12>))
(3) ((<NArray|URL:matrix.html#3>))
(4) ((<Special matrices|URL:matrix.html#4>))

== Class methods

--- GSL::Matrix.alloc(n)
--- GSL::Matrix.alloc(size1, size2)
--- GSL::Matrix.alloc(array) 
--- GSL::Matrix.alloc(arrays)
--- GSL::Matrix.alloc(...)
--- GSL::Matrix[...]
    These methods create a (({GSL::Matrix})) object.

    (1) From arrays
         irb(main):001:0> require("rbgsl")
         => true
         irb(main):002:0> m = GSL::Matrix[[1, 2, 3], [4, 5, 6], [7, 8, 9]]
         => GSL::Matrix 
         [ 1.000e+00 2.000e+00 3.000e+00 
           4.000e+00 5.000e+00 6.000e+00 
           7.000e+00 8.000e+00 9.000e+00 ]
      
    (2) With an array and rows&cols,
            m = GSL::Matrix.alloc([1, 2, 3, 4, 5, 6, 7, 8, 9], 3, 3)

    (3) With Range objects,
            irb(main):002:0> m = GSL::Matrix.alloc(1..3, 4..6, 7..9)
            [ 1.000e+00 2.000e+00 3.000e+00 
              4.000e+00 5.000e+00 6.000e+00 
              7.000e+00 8.000e+00 9.000e+00 ]
            irb(main):004:0> m2 = GSL::Matrix[1..6, 2, 3]
            [ 1.000e+00 2.000e+00 3.000e+00 
              4.000e+00 5.000e+00 6.000e+00 ]

--- GSL::Matrix.eye(n)
--- GSL::Matrix.eye(n1, n2)
    Examples:
      irb(main):008:0> m = GSL::Matrix::Int.eye(3)
      => GSL::Matrix::Int
      [ 1 0 0 
        0 1 0 
        0 0 1 ]
      irb(main):009:0> m = GSL::Matrix::Int.eye(2, 4)
      => GSL::Matrix::Int
      [ 1 0 0 0 
        0 1 0 0 ]
 
--- GSL::Matrix.identity(n)
--- GSL::Matrix.scalar(n)
--- GSL::Matrix.unit(n)
--- GSL::Matrix.I(n)
    Create diagonal matrices of dimensions n*n, of values 1.0.
    
--- GSL::Matrix.diagonal(a, b, c, ...)
--- GSL::Matrix.diagonal(Ary)
--- GSL::Matrix.diagonal(Range)
--- GSL::Matrix.diagonal(Vector)
    Creates a diagonal matrix of given elements.

    Example:
      irb(main):011:0> GSL::Matrix::Int.diagonal(1..4)
      => GSL::Matrix::Int
      [ 1 0 0 0 
        0 2 0 0 
        0 0 3 0 
        0 0 0 4 ]
      irb(main):012:0> GSL::Matrix::Int.diagonal(2, 5, 3)
      => GSL::Matrix::Int
      [ 2 0 0 
        0 5 0 
        0 0 3 ]

--- GSL::Matrix.ones(n)
--- GSL::Matrix.ones(n1, n2)
    Create a matrix of all the elements 1.

--- GSL::Matrix.zeros(n)
--- GSL::Matrix.zeros(n1, n2)
    Create a matrix of all the elements 1.

--- GSL::Matrix.indgen(n1, n2, start=0, step=1)
    Example:

      irb(main):016:0> m = GSL::Matrix::Int.indgen(3, 5)
      => GSL::Matrix::Int 
      [  0  1  2  3  4 
         5  6  7  8  9 
        10 11 12 13 14 ]
      irb(main):017:0> m = GSL::Matrix::Int.indgen(3, 5, 2)
      => GSL::Matrix::Int 
      [  2  3  4  5  6 
         7  8  9 10 11 
        12 13 14 15 16 ]
      irb(main):018:0> m = GSL::Matrix::Int.indgen(3, 5, 2, 3)
      => GSL::Matrix::Int 
      [  2  5  8 11 14 
        17 20 23 26 29 
        32 35 38 41 44 ]

=== NOTE:
Matrix dimensions are limited within the range of Fixnum.
For 32-bit CPU, the maximum of matrix dimension is 2^30 ~ 1e9.

== Instance Methods
=== Accessing matrix elements 

--- GSL::Matrix#size1
    Returns the number of rows of matrix ((|self|)).
--- GSL::Matrix#size2
    Returns the number of columns of matrix ((|self|)).
--- GSL::Matrix#shape
    Returns the number of rows and columns as an array.

    Ex:

       irb(main):005:0> m.size1
       => 3
       irb(main):006:0> m.size2
       => 5
       irb(main):007:0> m.shape
       => [3, 5]

--- GSL::Matrix#set(argv)
--- GSL::Matrix#[]=
    This method sets elements of the matrix.

--- GSL::Matrix#get(i, j)
--- GSL::Matrix#[]
    This method returns the ((|(i,j)|))-th element of the matrix ((|self|)). 

    Examples:
      irb(main):002:0> m = GSL::Matrix[1..9, 3, 3]
      => GSL::Matrix
      [ 1.000e+00 2.000e+00 3.000e+00 
        4.000e+00 5.000e+00 6.000e+00 
        7.000e+00 8.000e+00 9.000e+00 ]
      irb(main):003:0> m[1, 2]
      => 6.0
      irb(main):004:0> m[1, 2] = 123     # m.set(1, 2, 123)
      irb(main):005:0> m
      => GSL::Matrix
      [ 1.000e+00 2.000e+00 3.000e+00 
        4.000e+00 5.000e+00 1.230e+02 
        7.000e+00 8.000e+00 9.000e+00 ]
      irb(main):006:0> m[1]
      => GSL::Vector::View
      [ 4.000e+00 5.000e+00 6.000e+00 ]
      irb(main):007:0> m.set([3, 5, 2], [4, 5, 3], [7, 1, 5])
      => GSL::Matrix
      [ 3.000e+00 5.000e+00 2.000e+00 
        4.000e+00 5.000e+00 3.000e+00 
        7.000e+00 1.000e+00 5.000e+00 ]

=== Initializing matrix elements 
--- GSL::Matrix#set_all(x)
    This method sets all the elements of the matrix ((|self|)) to the value x.

--- GSL::Matrix#set_zero
    This method sets all the elements of the matrix to zero. 

--- GSL::Matrix#set_identity
    This method sets the elements of the matrix to the corresponding 
    elements of the identity matrix, i.e. a unit diagonal with all off-diagonal 
    elements zero. This applies to both square and rectangular matrices. 

=== IO
--- GSL::Matrix#fwrite(io)
--- GSL::Matrix#fwrite(filename)
--- GSL::Matrix#fread(io)
--- GSL::Matrix#fread(filename)
--- GSL::Matrix#fprintf(io, format = "%e")
--- GSL::Matrix#fprintf(filename, format = "%e")
--- GSL::Matrix#fscanf(io)
--- GSL::Matrix#fscanf(filename)

=== Matrix views
The (({GSL::Matrix::View})) class is defined to be used as "references" to
matrices. The (({Matrix::View})) class is a subclass of (({Matrix})),
and an instance of the (({View})) class created by slicing a (({Matrix})) object 
can be used same as the original matrix. The
(({View})) object shares the data with the original matrix, i.e. any changes
in the elements of the (({View})) object affect to the original.

--- GSL::Matrix#submatrix(k1, k2, n1, n2)
--- GSL::Matrix#view(k1, k2, n1, n2)
    This returns a (({GSL::Matirx::View})) object, a submatrix of the matrix 
    ((|self|)). The upper-left element of the submatrix is the element ((|(k1,k2)|)) 
    of the original matrix. The submatrix has ((|n1|)) rows and ((|n2|)) columns. 

--- GSL::Vectir#matrix_view(n1, n2)
    This creates a (({Matrix::View})) object from the vector ((|self|)).

    Ex:
       irb(main):002:0> v = Vector[1..9]
       => GSL::Vector 
       [ 1.000e+00 2.000e+00 3.000e+00 4.000e+00 5.000e+00 6.000e+00 7.000e+00 8.000e+00 9.000e+00 ]
       irb(main):003:0> m = v.matrix_view(3, 3)
       => GSL::Matrix::View 
       [ 1.000e+00 2.000e+00 3.000e+00 
         4.000e+00 5.000e+00 6.000e+00 
         7.000e+00 8.000e+00 9.000e+00 ]
       irb(main):004:0> m[1][1] = 99.99
       => 99.99
       irb(main):005:0> v
       => GSL::Vector 
       [ 1.000e+00 2.000e+00 3.000e+00 4.000e+00 9.999e+01 6.000e+00 7.000e+00 8.000e+00 9.000e+00 ]
       irb(main):006:0> 


=== Creating row and column views

--- GSL::Matrix#row(i)
    These methods return ((|i|))-th row of the matrix as a (({Vector::View}))
    object. Any modifications to the (({Vectror::View})) object returned by this method
    propagate to the original matrix. 
    
    By using this property, it is possible to access to matrix elements 
    and also modify them, with the expressions as (({a = m[i][j]})) 
    and (({m[i][j] = val})). 
    The results are just same as using the methods (({Matrix#get(i, j)}))
    and (({Matrix#set(i, j, val)})), but different in manner. 
    The methods (({get})) and (({set})) use
    the GSL C functions (({gsl_matrix_get(), gsl_matrix_set()})), i.e. access
    to the (i, j)-element directly. In the expression (({m[i][j]})), first a
    (({Vector::View})) object which points to the data of the i-th row of the
    matrix (({m})) is created, (({m[i]})), 
    and then the j-th element of the (({Vector::View})) object (({m[i]})),
    expressed as (({m[i][j]})), is extracted/modified.

--- GSL::Matrix#column(i)
--- GSL::Matrix#col(i) 
    These methods return a vector view of the ((|j|))-th column of the matrix.

--- GSL::Matrix#subrow(i, offset, n)
    Returns a vector view of the ((|i|))-th row of the matrix ((|self|)) 
    beginning at ((|offset|)) elements past the first column 
    and containing ((|n|)) elements. (>= GSL-1.10)

--- GSL::Matrix#subcolumn(j, offset, n)
    Returns a vector view of the ((|j|))-th column of the matrix ((|self|)) 
    beginning at ((|offset|)) elements past the first row 
    and containing ((|n|)) elements. (>= GSL-1.10)

--- GSL::Matrix#diagonal 
    This method returns a (({Vector::View})) of the diagonal of the matrix.
    The matrix is not required to be square. For a rectangular matrix the 
    length of the diagonal is the same as the smaller dimension of the matrix.


    Ex:
      irb(main):017:0> m = GSL::Matrix[1..9, 3, 3]
      => GSL::Matrix 
      [ 1.000e+00 2.000e+00 3.000e+00 
        4.000e+00 5.000e+00 6.000e+00 
        7.000e+00 8.000e+00 9.000e+00 ]
      irb(main):018:0> m.row(1)
      => GSL::Vector::View 
      [ 4.000e+00 5.000e+00 6.000e+00 ]
      irb(main):019:0> m.col(2)
      => GSL::Vector::Col::View 
      [ 3.000e+00 
        6.000e+00 
        9.000e+00 ]
      irb(main):020:0> m.col(2)[2] = 123
      => 123
      irb(main):021:0> m
      => GSL::Matrix 
      [ 1.000e+00 2.000e+00 3.000e+00 
        4.000e+00 5.000e+00 6.000e+00 
        7.000e+00 8.000e+00 1.230e+02 ]
      irb(main):022:0> m.diagonal
      => GSL::Vector::View: 
      [ 1.000e+00 5.000e+00 1.230e+02 ]

--- GSL::Matrix#subdiagonal(k)
    Returns a vector view view of the ((|k|))-th subdiagonal 
    of the matrix ((|self|)). 
    The matrix is not required to be square. The diagonal of the matrix 
    corresponds to k = 0.

--- GSL::Matrix#superdiagonal(k)
    Returns a vector view of the ((|k|))-th superdiagonal of the matrix ((|self|)). 
    The matrix is not required to be square. The diagonal of the matrix 
    corresponds to k = 0.

--- GSL::Matrix#to_v
    Creates a (({GSL::Vector})) object "flattening" the rows of the matrix ((|self|)).

        irb(main):002:0> m = GSL::Matrix[1..6, 2, 3]
        => GSL::Matrix 
        [ 1.000e+00 2.000e+00 3.000e+00 
          4.000e+00 5.000e+00 6.000e+00 ]
        irb(main):003:0> m.to_v
        => GSL::Vector 
        [ 1.000e+00 2.000e+00 3.000e+00 4.000e+00 5.000e+00 6.000e+00 ]

=== Iterators
--- GSL::Matrix#each_row
    Iterator for each of rows in the matrix ((|self|)). 
--- GSL::Matrix#each_col
    Iterator for each of columns in the matrix ((|self|)). 

--- GSL::Matrix#collect { |item| .. }

=== Copying matrices
--- GSL::Matrix#clone
--- GSL::Matrix#duplicate
    Create a new matrix of the same elements.

--- GSL::Matrix.memcpy(dest, src)
--- GSL::Matrix.swap(dest, src)

=== Copying rows and columns

--- GSL::Matrix#get_row(i)
    This method returns a new vector (not a view) which contains the elements 
    of the ((|i|))-th row of the matrix ((|self|)).

--- GSL::Matrix#get_col(j)
    This method returns a new vector (not a view) which contains the elements of the ((|j|))-th 
    column of the matrix ((|self|)).

--- GSL::Matrix#set_row(i, v)
    This method copies the elements of the vector ((|v|)) into the ((|i|))-th 
    row of the matrix. 
    The length of the vector must be the same as the length of the row. 

--- GSL::Matrix#set_col(j, v)
    This method copies the elements of the vector ((|v|)) into the ((|j|))-th 
    column of the matrix. The length of the vector must be the same as the length 
    of the column. 

=== Exchanging rows and columns
--- GSL::Matrix#swap_rows!(i, j)
    This method exchanges the ((|i|))-th and ((|j|))-th rows of the matrix ((|in-place|)).
--- GSL::Matrix#swap_rows(i, j)
    This method creates a new matrix exchanging the ((|i|))-th and ((|j|))-th rows of the matrix ((|self|)).

--- GSL::Matrix#swap_columns!(i, j)
    This method exchanges the ((|i|))-th and ((|j|))-th columns of the matrix ((|in-place|)). 
--- GSL::Matrix#swap_columns(i, j)
    This method creates a new matrix exchanging the ((|i|))-th and ((|j|))-th columns of the matrix ((|self|)).

--- GSL::Matrix#swap_rowcol(i, j)
    This method exchanges the ((|i|))-th row and ((|j|))-th column of the matrix.
    The matrix must be square for this operation to be possible.

--- GSL::Matrix#transpose_memcpy
--- GSL::Matrix#transpose
    This method returns a matrix of a transpose of the matrix. The matrix
    ((|self|)) is not modified.

--- GSL::Matrix#transpose!
    This method replaces the matrix by its transpose by copying the 
    elements of the matrix ((|in-place|)). The matrix must be square for this
    operation to be possible. 

--- GSL::Matrix#reverse_rows
--- GSL::Matrix#flipud
    Example:
        irb(main):018:0> m = GSL::Matrix::Int[1..9, 3, 3]
        => GSL::Matrix::Int 
        [ 1 2 3 
          4 5 6 
          7 8 9 ]
        irb(main):019:0> m.reverse_rows
        => GSL::Matrix::Int 
        [ 7 8 9 
          4 5 6 
          1 2 3 ]

--- GSL::Matrix#reverse_columns
--- GSL::Matrix#fliplr
    Example:
        irb(main):018:0> m = GSL::Matrix::Int[1..9, 3, 3]
        => GSL::Matrix::Int
        [ 1 2 3 
          4 5 6 
          7 8 9 ]
        irb(main):020:0> m.reverse_rows.reverse_columns
        => GSL::Matrix::Int 
        [ 9 8 7 
          6 5 4 
          3 2 1 ]

--- GSL::Matrix#rot90(n = 1)
    Return a copy of ((|self|)) with the elements rotated 
    counterclockwise in 90-degree increments. The argument ((|n|)) is 
    optional, and specifies how many 90-degree rotations are to be applied 
    (the default value is 1). 
    Negative values of ((|n|)) rotate the matrix in a clockwise direction.

    Examples:
      irb(main):014:0> m = GSL::Matrix::Int[1..6, 2, 3]
      => GSL::Matrix::Int
      [ 1 2 3 
        4 5 6 ]
      irb(main):015:0> m.rot90
      => GSL::Matrix::Int
      [ 3 6 
        2 5 
        1 4 ]
      irb(main):016:0> m.rot90(2)
      => GSL::Matrix::Int
      [ 6 5 4 
        3 2 1 ]
      irb(main):017:0> m.rot90(3)
      => GSL::Matrix::Int
      [ 4 1 
        5 2 
        6 3 ]
      irb(main):018:0> m.rot90(-1)
      => GSL::Matrix::Int
      [ 4 1 
        5 2 
        6 3 ]

--- GSL::Matrix#upper
    This creates a matrix copying the upper half part of the matrix 
    ((|self|)), including the diagonal elements.
--- GSL::Matrix#lower
    This creates a matrix copying the lower half part of the matrix 
    ((|self|)), including the diagonal elements.

         irb(main):014:0> m = GSL::Matrix[1..9, 3, 3]
         => GSL::Matrix
         [ 1.000e+00 2.000e+00 3.000e+00 
           4.000e+00 5.000e+00 6.000e+00 
           7.000e+00 8.000e+00 9.000e+00 ]
         irb(main):015:0> m.upper
         => GSL::Matrix 
         [ 1.000e+00 2.000e+00 3.000e+00 
           0.000e+00 5.000e+00 6.000e+00 
           0.000e+00 0.000e+00 9.000e+00 ]
         irb(main):016:0> m.lower
         => GSL::Matrix 
         [ 1.000e+00 0.000e+00 0.000e+00 
           4.000e+00 5.000e+00 0.000e+00 
           7.000e+00 8.000e+00 9.000e+00 ]

--- GSL::Matrix#horzcat(other)
    Returns the horizontal concatenation of ((|self|)) and ((|other|)).

    Ex:
      irb(main):001:0> require("rbgsl")
      => true
      irb(main):002:0> a = GSL::Matrix::Int[1..4, 2, 2]
      => GSL::Matrix::Int
      [ 1 2 
        3 4 ]
      irb(main):003:0> b = GSL::Matrix::Int[5..10, 2, 3]
      => GSL::Matrix::Int
      [  5  6  7 
         8  9 10 ]
      irb(main):004:0> a.horzcat(b)
      => GSL::Matrix::Int
      [  1  2  5  6  7 
         3  4  8  9 10 ]

--- GSL::Matrix#vertcat(other)
    Returns the vertical concatenation of ((|self|)) and ((|other|)).
 
    Ex:
      irb(main):002:0> a = GSL::Matrix::Int[1..4, 2, 2]
      => GSL::Matrix::Int
      [ 1 2 
        3 4 ]
      irb(main):003:0> b = GSL::Matrix::Int[5..10, 3, 2]
      => GSL::Matrix::Int
      [  5  6 
         7  8 
         9 10 ]
      irb(main):004:0> a.vertcat(b)
      => GSL::Matrix::Int
      [  1  2 
         3  4 
         5  6 
         7  8 
         9 10 ]

=== Matrix operations

--- GSL::Matrix#add(b)
--- GSL::Matrix#+(b)
    This method adds the elements of matrix ((|b|)) 
    to the elements of the  matrix. 
    The two matrices must have the same dimensions. 

    If ((|b|)) is a scalar, these methods add it to all the elements
    of the matrix ((|self|)) (equivalent to the method (({add_constant}))).

--- GSL::Matrix#sub(b)
--- GSL::Matrix#-(b)
    This method subtracts the elements of matrix ((|b|)) 
    from the elements of the  
    matrix. The two matrices must have the same dimensions. 

--- GSL::Matrix#mul_elements(b)
    This method multiplies the elements of the matrix by the elements of 
    matrix ((|b|)). The two matrices must have the same dimensions. 
    If ((|b|)) is a scalar, the method (({scale})) (see below) 
    is called.

--- GSL::Matrix#div_elements(b)

    This method divides the elements of the  matrix by the elements of 
    matrix ((|b|)). The two matrices must have the same dimensions. 

--- GSL::Matrix#scale(x)
    This method multiplies the elements of the  matrix by the constant 
    factor ((|x|)).

--- GSL::Matrix#add_constant(x)
    This method adds the constant value ((|x|)) to the elements of the matrix.

--- GSL::Matrix#*(b)
    Matrix multiplication. 

    Ex:

       irb(main):002:0> a = GSL::Matrix[1..4, 2, 2]
       => GSL::Matrix 
       [ 1.000e+00 2.000e+00 
         3.000e+00 4.000e+00 ]
       irb(main):003:0> b = GSL::Matrix[5..8, 2, 2]
       => GSL::Matrix 
       [ 5.000e+00 6.000e+00 
         7.000e+00 8.000e+00 ]
       irb(main):004:0> a*b
       => GSL::Matrix 
       [ 1.900e+01 2.200e+01 
         4.300e+01 5.000e+01 ]
       irb(main):005:0> a*2
       => GSL::Matrix 
       [ 2.000e+00 4.000e+00 
         6.000e+00 8.000e+00 ]
       irb(main):006:0> c = Vector[1, 2]
       => GSL::Vector 
       [ 1.000e+00 2.000e+00 ]
       irb(main):007:0> a*c.col
       => GSL::Vector::Col 
       [ 5.000e+00 
         1.100e+01 ]

--- GSL::Matrix#/(b)
    If ((|b|)) is a scalar or a (({Matrix})), this method calculates the
    element-by-element divisions.
    If a (({Vector::Col})) is given, this method solves the linear system
    by using LU decomposition.

    Ex: 
        irb(main):002:0> m = GSL::Matrix[1..4, 2, 2]
        => GSL::Matrix 
        [ 1.000e+00 2.000e+00   
          3.000e+00 4.000e+00 ]
        irb(main):003:0> m/3
        => GSL::Matrix 
        [ 3.333e-01 6.667e-01         <--- 1/3, 2/3
          1.000e+00 1.333e+00 ]       <--- 3/3, 4/3
        irb(main):004:0> b = Vector[5, 6].col
        => GSL::Vector::Col 
        [ 5.000e+00   
          6.000e+00 ]
        irb(main):005:0> x = m/b          <--- Solve m (x,y) = b
        => GSL::Vector::Col 
        [ -4.000e+00                  <--- x = -4
          4.500e+00 ]                 <--- y = 4.5 
        irb(main):006:0> m*x
        => GSL::Vector::Col
        [  5.000e+00 
           6.000e+00 ]

--- GSL::Matrix#^(b)
    Computes matrix power of ((|b|)).

=== Finding maximum and minimum elements of matrices

--- GSL::Matrix#max
--- GSL::Matrix#min
    These methods return the max/min value in the  matrix.

--- GSL::Matrix#minmax 
    This method returns a two elements array [min, max], 
    which contains the minimum
    and the maximum values in the  matrix.

--- GSL::Matrix#max_index 
--- GSL::Matrix#min_index 
    These methods return the index of the max/min value in the matrix.

--- GSL::Matrix#minmax_index 
    This method returns a two elements array [imin, imax], 
    which contains the indices
    of the minimum and the maximum value in the matrix.

=== Matrix properties
--- GSL::Matrix#isnull
    This returns 1 if all the elements of the matrix ((|self|)) are zero, 
    and 0 otherwise.

--- GSL::Matrix#isnull?
    This returns (({true})) if all the elements of the matrix ((|self|)) 
    are zero, and (({false})) otherwise.

--- GSL::Matrix#ispos
--- GSL::Matrix#ispos?
    (GSL-1.9 or later) Return 1 (true) if all the elements of the matrix ((|self|)) are strictly positive, and 0 (false) otherwise. 

--- GSL::Matrix#isneg
--- GSL::Matrix#isneg?
    (GSL-1.9 or later) Return 1 (true) if all the elements of the matrix ((|self|)) are strictly negative, and 0 (false) otherwise. 

--- GSL::Matrix#isnonneg
--- GSL::Matrix#isnonneg?
    (GSL-1.10 or later) Return 1 (true) if all the elements of the matrix ((|self|)) are non-negative , and 0 (false) otherwise. 

--- GSL::Matrix#any
    Returns a Vector of ones and zeros with each element indicating 
    whether any of the elements of the corresponding column of the 
    matrix are nonzero.

--- GSL::Matrix#all
    Behaves like the method (({any})), except that it returns 1 only if 
    all the elements of the matrix.

--- GSL:Matrix#trace
    This returns trace of the matrix ((|self|)), the sum of the diagonal 
    elements.

--- GSL:Matrix#norm
    Returns matrix norm, sqrt(sum_{ij} m_{ij}^2).

--- GSL:Matrix#abs
    Example:
      irb(main):004:0> m = GSL::Matrix::Int[-5..4, 3, 3]
      => GSL::Matrix::Int
      [ -5 -4 -3 
        -2 -1  0 
         1  2  3 ]
      irb(main):005:0> m.abs
      => GSL::Matrix::Int
      [ 5 4 3 
        2 1 0 
        1 2 3 ]

--- GSL::Matrix#equal?(other, eps = 1e-10)
--- GSL::Matrix#==(other, eps = 1e-10)
    Returns (({true})) if the matrices have same size and elements
    equal to absolute accurary ((|eps|)) for all the indices, 
    and (({false})) otherwise.

== NArray

--- GSL::Matrix#to_na
    The Matrix object ((|self|)) is converted into an (({NMatrix})) object. 
    The matrix data are copied to newly allocated memory.

--- NArray#to_gm
--- NArray#to_gslm
    Convert (({NArray})) object into (({GSL::Matrix})).

--- NArray#to_gm_view
--- NArray#to_gslm_view
    A (({GSL::Matrix::View})) object is created from the NArray object ((|na|)). 
    The data of ((|na|)) are 
    not copied, thus any modifications to the View object affect on the original 
    NArray object ((|na|)). 
    The View object can be used as a reference to the NMatrix object.

== Special matrices
--- GSL::Matrix.hirbert(n)
    Returns the Hilbert matrix of order ((|n|)). The ((|ij|)) element is
    defined as 1/(i+j+1).

--- GSL::Matrix.invhirbert(n)
    Returns the inverse of a Hilbert matrix of order ((|n|)).

    Ex:
      irb(main):009:0> m = GSL::Matrix.hilbert(4)
      => GSL::Matrix 
      [ 1.000e+00 5.000e-01 3.333e-01 2.500e-01 
        5.000e-01 3.333e-01 2.500e-01 2.000e-01 
        3.333e-01 2.500e-01 2.000e-01 1.667e-01 
        2.500e-01 2.000e-01 1.667e-01 1.429e-01 ]
      irb(main):010:0> invm = GSL::Matrix.invhilbert(4)
      => GSL::Matrix 
      [ 1.600e+01 -1.200e+02 2.400e+02 -1.400e+02 
        -1.200e+02 1.200e+03 -2.700e+03 1.680e+03 
        2.400e+02 -2.700e+03 6.480e+03 -4.200e+03 
        -1.400e+02 1.680e+03 -4.200e+03 2.800e+03 ]
      irb(main):011:0> invm2 = m.inv
      => GSL::Matrix 
      [ 1.600e+01 -1.200e+02 2.400e+02 -1.400e+02 
        -1.200e+02 1.200e+03 -2.700e+03 1.680e+03 
        2.400e+02 -2.700e+03 6.480e+03 -4.200e+03 
        -1.400e+02 1.680e+03 -4.200e+03 2.800e+03 ]
      irb(main):012:0> m*invm
      => GSL::Matrix 
      [ 1.000e+00 5.684e-14 -2.274e-13 1.137e-13 
        1.998e-15 1.000e+00 -4.663e-14 3.109e-14 
        3.664e-15 -7.239e-14 1.000e+00 -1.017e-13 
        -2.442e-15 1.510e-14 -8.038e-14 1.000e+00 ]
      irb(main):013:0> m*invm2
      => GSL::Matrix 
      [ 1.000e+00 0.000e+00 0.000e+00 0.000e+00 
        -1.554e-15 1.000e+00 -2.389e-14 8.349e-15 
        1.295e-15 3.405e-15 1.000e+00 -6.957e-15 
        1.110e-15 1.916e-14 1.707e-14 1.000e+00 ]

--- GSL::Matrix.pascal(n)
    Returns the Pascal matrix of order ((|n|)), created from Pascal's triangle.

      irb(main):002:0> GSL::Matrix::Int.pascal(10)
      => GSL::Matrix::Int 
      [     1     1     1     1     1     1     1     1     1     1 
            1     2     3     4     5     6     7     8     9    10 
            1     3     6    10    15    21    28    36    45    55 
            1     4    10    20    35    56    84   120   165   220 
            1     5    15    35    70   126   210   330   495   715 
            1     6    21    56   126   252   462   792  1287  2002 
            1     7    28    84   210   462   924  1716  3003  5005 
            1     8    36   120   330   792  1716  3432  6435 11440 
            1     9    45   165   495  1287  3003  6435 12870 24310 
            1    10    55   220   715  2002  5005 11440 24310 48620 ]

--- GSL::Matrix.vandermonde(v)
    Creates a Vendermonde matrix from a vector or an array ((|v|)).
 
       irb(main):002:0> GSL::Matrix.vander([1, 2, 3, 4])
       => GSL::Matrix
       [ 1.000e+00 1.000e+00 1.000e+00 1.000e+00
         8.000e+00 4.000e+00 2.000e+00 1.000e+00
         2.700e+01 9.000e+00 3.000e+00 1.000e+00
         6.400e+01 1.600e+01 4.000e+00 1.000e+00 ]

--- GSL::Matrix.toeplitz(v)
    Creates a Toeplitz matrix from a vector or an array ((|v|)).
 
        irb(main):004:0> GSL::Matrix::Int.toeplitz([1, 2, 3, 4, 5])
        => GSL::Matrix::Int 
        [ 1 2 3 4 5 
          2 1 2 3 4 
          3 2 1 2 3 
          4 3 2 1 2 
          5 4 3 2 1 ]

--- GSL::Matrix.circulant(v)
    Creates a circulant matrix from a vector or an array ((|v|)).
 
        irb(main):005:0> GSL::Matrix::Int.circulant([1, 2, 3, 4])
        => GSL::Matrix::Int 
        [ 4 1 2 3 
          3 4 1 2 
          2 3 4 1 
          1 2 3 4 ]

((<prev|URL:vector.html>))
((<next|URL:perm.html>))

((<Reference index|URL:ref.html>))
((<top|URL:index.html>))
    
=end
