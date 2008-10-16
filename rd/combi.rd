=begin
= Combinations
Contents:
(1) ((<Combination allocation|URL:combi.html#1>))
(2) ((<Methods|URL:combi.html#2>))
    (1) ((<Accessing combination elements|URL:combi.html#2.1>))
    (2) ((<Combination properties|URL:combi.html#2.2>))
    (3) ((<Combination functions|URL:combi.html#2.3>))
    (4) ((<Reading and writing combinations|URL:combi.html#2.4>))

== Combination allocation
--- GSL::Combination.alloc(n, k)
    These create a new combination with parameters ((|n, k|)). 
    The combination is not initialized and its elements are undefined. 
    Use the method (({GSL::Combination.calloc})) if you want to create a 
    combination which is initialized to the lexicographically first combination. 

--- GSL::Combination.calloc(n, k)
    This creates a new combination with parameters ((|n, k|)) and initializes 
    it to the lexicographically first combination. 

== Methods

--- GSL::Combination#init_first
    This method initializes the combination ((|self|)) to the lexicographically 
    first combination, i.e. (0,1,2,...,k-1).

--- GSL::Combination#init_last
    This method initializes the combination ((|self|)) to the lexicographically last 
    combination, i.e. (n-k,n-k+1,...,n-1).

=== Accessing combination elements
--- GSL::Combination#get(i)
--- GSL::Combination#[i]
    This returns the value of the ((|i|))-th element of the combination ((|self|)). 

=== Combination properties
--- GSL::Combination#n
    Returns the (({n})) parameter of the combination ((|self|)).

--- GSL::Combination#k
    Returns the (({k})) parameter of the combination ((|self|)).

--- GSL::Combination#data
    Returns the vector of elements in the combination ((|self|)).

--- GSL::Combination#valid
    This method checks that the combination ((|self|)) is valid. 
    The (({k})) elements should contain numbers from range 0 .. n-1, 
    each number at most once. The numbers have to be in increasing order.

--- GSL::Combination#valid?
    Thie returns true if the combination is valid, and false otherwise.

=== Combination functions
--- GSL::Combination#next
    This method advances the combination ((|self|)) to the next combination in 
    lexicographic order and returns (({GSL::SUCCESS})). If no further combinations are 
    available it returns (({GSL::FAILURE})) and leaves ((|self|)) unmodified. 
    Starting with the first combination and repeatedly applying this function will 
    iterate through all possible combinations of a given order.

--- GSL::Combination#prev
    This method steps backwards from the combination ((|self|)) to the previous 
    combination in lexicographic order, returning (({GSL::SUCCESS})). 
    If no previous combination is available it returns (({GSL::FAILURE})) 
    and leaves ((|self|)) unmodified.

=== Reading and writing combinations
--- GSL::Combination#fwrite(filename)
--- GSL::Combination#fwrite(io)
--- GSL::Combination#fread(filename)
--- GSL::Combination#fread(io)
--- GSL::Combination#fprintf(filename, format = "%u")
--- GSL::Combination#fprintf(io, format = "%u")
--- GSL::Combination#fscanf(filename)
--- GSL::Combination#fscanf(io)

== Example
     #!/usr/bin/env ruby
     require("rbgsl")

     printf("All subsets of {0,1,2,3} by size:\n") ;
     for i in 0...4 do
       c = GSL::Combination.calloc(4, i);
       begin
         printf("{");
         c.fprintf(STDOUT, " %u");
         printf(" }\n");
       end while c.next == GSL::SUCCESS
     end

((<prev|URL:perm.html>))
((<next|URL:sort.html>))

((<Reference index|URL:ref.html>))
((<top|URL:index.html>))

=end
