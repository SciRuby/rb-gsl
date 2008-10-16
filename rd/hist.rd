=begin
= Histograms
(1) ((<Histogram allocation|URL:hist.html#1>))
(2) ((<Copying histograms|URL:hist.html#2>))
(3) ((<Updating and accessing histogram elements|URL:hist.html#3>))
(4) ((<Searching histogram ranges|URL:hist.html#4>))
(5) ((<Histogram Statistics|URL:hist.html#5>))
(6) ((<Histogram Operations|URL:hist.html#6>))
(7) ((<Reading and writing histograms|URL:hist.html#7>))
(8) ((<Extensions|URL:hist.html#8>))
    (1) ((<Histogram Operations|URL:hist.html#8.1>))
    (2) ((<Graph interface|URL:hist.html#8.2>))
    (3) ((<Histogram Fittings|URL:hist.html#8.3>))
(9) ((<The histogram probability distribution|URL:hist.html#9>))

== Histogram allocation
--- GSL::Histogram.alloc(n)
--- GSL::Histogram.alloc(n, [xmin, xmax])
--- GSL::Histogram.alloc(n, xmin, xmax)
--- GSL::Histogram.alloc(n)
--- GSL::Histogram.alloc(array)
--- GSL::Histogram.alloc(vector)
    Constructor for a histogram object with ((|n|)) bins. 

    Examples:

    (1) With an integer:
         h = Histogram.alloc(4)  <--- Histogram of 4 bins.
                                      The range is not defined yet.

                   [ bin[0] )[ bin[1] )[ bin[2] )[ bin[3] )
                   |---------|---------|---------|---------|
                range[0]  range[1]  range[2]  range[3]  range[4]
  
    (2) With an array or a vector:
         h = Histogram.alloc([1, 3, 7, 9, 20])  <--- Histogram of 4 bins.
                                                   The range is initialized as
                                                   range[0] = 1, range[1] = 3, ..., range[4] = 20.
   
    (3) With size and the range [min, max]:

          irb(main):004:0> h = Histogram.alloc(5, [0, 5])
          irb(main):005:0> h.range
          => GSL::Histogram::Range: 
          [ 0.000e+00 1.000e+00 2.000e+00 3.000e+00 4.000e+00 5.000e+00 ]
          irb(main):006:0> h.bin
          => GSL::Histogram::Bin: 
          [ 0.000e+00 0.000e+00 0.000e+00 0.000e+00 0.000e+00 ]
          irb(main):007:0> h.increment(2.5)
          irb(main):008:0> h.bin
          => GSL::Histogram::Bin: 
          [ 0.000e+00 0.000e+00 1.000e+00 0.000e+00 0.000e+00 ]

--- GSL::Histogram.alloc_uniform(n, min, max)
--- GSL::Histogram.alloc_uniform(n, [min, max])
--- GSL::Histogram.equal_bins_p(h1, h2)
--- GSL::Histogram.equal_bins(h1, h2)
    Return 1 if the all of the individual bin ranges of the two histograms 
    are identical, and 0 otherwise.
--- GSL::Histogram.equal_bins_p?(h1, h2)
--- GSL::Histogram.equal_bins?(h1, h2)
    Return ((|true|)) if the all of the individual bin ranges of the two histograms 
    are identical, and ((|false|)) otherwise.

--- GSL::Histogram#set_ranges(v)
    This sets the ranges of the existing histogram using a ((<GSL::Vector|URL:vector.html>)) object.
--- GSL::Histogram#set_ranges_uniform(xmin, xmax)
--- GSL::Histogram#set_ranges_uniform([xmin, xmax])
    This method sets the ranges of the existing histogram ((|self|)) 
    to cover the range ((|xmin|)) to ((|xmax|)) uniformly. 
    The values of the histogram bins are reset to zero. 
    The bin ranges are shown as below,
       bin[0] corresponds to xmin <= x < xmin + d
       bin[1] corresponds to xmin + d <= x < xmin + 2 d
         ......
       bin[n-1] corresponds to xmin + (n-1)d <= x < xmax
    where d is the bin spacing, d = (xmax-xmin)/n.

== Copying Histograms
--- GSL::Histogram.memcpy(dest, src)
    Copies the histogram ((|src|)) into the pre-existing histogram ((|dest|)), 
    making dest into an exact copy of ((|src|)). 
    The two histograms must be of the same size. 
--- GSL::Histogram#clone
    Returns a newly created histogram which is an exact copy of the histogram 
    ((|self|)).

== Updating and accessing histogram elements
--- GSL::Histogram#increment(x, weight = 1)
--- GSL::Histogram#fill(x, weight = 1)
--- GSL::Histogram#accumulate(x, weight = 1)
    These methods updates the histogram ((|self|)) by adding ((|weight|))
    (default = 1) to the bin whose range contains the coordinate ((|x|)).
    If ((|x|)) is an instance of (({GSL::Vector})) or (({Array})), 
    all the elements are filled into the histogram.
    If ((|x|)) is less than (greater than) the lower limit (upper limit) 
    of the histogram then none of bins are modified. 

--- GSL::Histogram#increment2(x, weight = 1)
--- GSL::Histogram#fill2(x, weight = 1)
--- GSL::Histogram#accumulate2(x, weight = 1)
    These methods updates the histogram ((|self|)) by adding ((|weight|))
    to the bin whose range contains the coordinate ((|x|)).    
    If ((|x|)) is less than the lower limit, the lowest bin is incremented.
    If ((|x|)) is greater than the upper limit, the highest bin is incremented.

--- GSL::Histogram#get(i)
--- GSL::Histogram#[i]
    These methods return the contents of the ((|i|))-th bin of the histogram 
    ((|self|)).

--- GSL::Hiatogram#get_range(i)
    This method finds the upper and lower range limits of the ((|i|))-th bin 
    of the histogram ((|self|)), and returns an array [((|lower, upper|))].

--- GSL::Histogram#range
    This returns a (({Vector::View})) object as a reference to the pointer 
    (({double *range})) in the (({gsl_histogram})) struct.

--- GSL::Histogram#bin
    This returns a (({Vector::View})) object to access the pointer (({double *bin})) in the (({gsl_histogram})) struct.

--- GSL::Histogram#max
--- GSL::Histogram#min
--- GSL::Histogram#bins
    These methods return the maximum upper and minimum lower range 
    limits and the number of bins of the histogram ((|self|)).

--- GSL::Histogram#reset
    This method resets all the bins in the histogram ((|self|)) to zero.

== Searching histogram ranges 
--- GSL::Histogram#find(x)
    This method finds and sets the index i to the bin number which 
    covers the coordinate ((|x|)) in the histogram ((|self|)).

== Histogram Statistics 
--- GSL::Histogram#max_val
    This returns the maximum value contained in the histogram bins.

--- GSL::Histogram#max_bin
    This returns the index of the bin containing the maximum value. 
    In the case where several bins contain the same maximum value the 
    smallest index is returned.

--- GSL::Histogram#min_val
    This returns the minimum value contained in the histogram bins.

--- GSL::Histogram#min_bin
    This returns the index of the bin containing the minimum value. 
    In the case where several bins contain the same maximum value 
    the smallest index is returned.

--- GSL::Histogram#mean
    This returns the mean of the histogrammed variable, 
    where the histogram is regarded as a probability distribution. 
    Negative bin values are ignored for the purposes of this calculation. 
    The accuracy of the result is limited by the bin width.

--- GSL::Histogram#sigma
    This function returns the standard deviation of the histogrammed variable, 
    where the histogram is regarded as a probability distribution. 
    Negative bin values are ignored for the purposes of this calculation. 
    The accuracy of the result is limited by the bin width.

--- GSL::Histogram#sum(istart = 0, iend = n-1)
    The sum of values of the histogram ((|self|)) from the ((|istart|))-th bin
    to the ((|iend|))-th bin.


== Histogram Operations

--- GSL::Histogram#add(h2)
--- GSL::Histogram#sub(h2)
--- GSL::Histogram#mul(h2)
--- GSL::Histogram#div(h2)
--- GSL::Histogram#scale(val)
--- GSL::Histogram#shift(val)


== Reading and writing histograms
--- GSL::Histogram#fwrite(io)
--- GSL::Histogram#fwrite(filename)
--- GSL::Histogram#fread(io)
--- GSL::Histogram#fread(filename)
--- GSL::Histogram#fprintf(io, range_format = "%e", bin_format = "%e")
--- GSL::Histogram#fprintf(filename, range_format = "%e", bin_format = "%e")
--- GSL::Histogram#fscanf(io)
--- GSL::Histogram#fscanf(filename)

== Extentions
=== Histogram operations
--- GSL::Histogram#normalize
    This methods scales the contents of the bins of histogram ((|self|)) by its
    maximum value.

--- GSL::Histogram#rebin(m = 2)
    This method creates a new histogram merging ((|m|)) bins in one in the
    histogram ((|self|)). This method cannot be used for histograms of 
    non-uniform bin size. If ((|m|)) is not an exact divider of the number 
    of bins of ((|self|)), the range of the rebinned histogram is extended
    not to lose the entries in the last ((|m-1|)) (at most) bins.

    Example: a histogram ((|h|)) of size 5 with the range [0, 5), binned as

        GSL::Histogram::Range: 
        [ 0.000e+00 1.000e+00 2.000e+00 3.000e+00 4.000e+00 5.000e+00 ]
        GSL::Histogram::Bin: 
        [ 0.000e+00 3.000e+00 1.000e+00 1.000e+00 3.000e+00 ]
 
    When a new histogram is created merging two bins into one as 
    (({h2 = h.rebin})), then (({h2})) looks like

        GSL::Histogram::Range: 
        [ 0.000e+00 2.000e+00 4.000e+00 6.000e+00 ]
        GSL::Histogram::Bin: 
        [ 3.000e+00 2.000e+00 3.000e+00 ]

--- GSL::Histogram#reverse
    This method create a new histogram reversing the order of the range and the bin of 
    histogram ((|self|)).
    
--- GSL::Histogram#integrate(istart = 0, iend = n-1)
--- GSL::Histogram#integrate([istart, iend])
--- GSL::Histogram#integrate(direction = 1 or -1)
    This method calculates cumulative counts of the histogram ((|self|))
    from the ((|istart|))-th bin to the ((|iend|))-th bin (((|iend|)) inclusive), 
    and returns a (({GSL::Histogram::Integral})) 
    object. If ((|istart <= iend|)) (or ((|direction == 1|))), 
    the ((|i|))-th bin value of a 
    (({GSL::Histogram::Integral})) object ((|hi|)) created from a 
    (({GSL::Histogram})) ((|h|)) is given by (({hi[i] = hi[i-1] + h[i]})).
    If ((|istart > iend|)) (or ((|direction == -1|))), (({hi[i] = hi[i+1] = h[i]})).

--- GSL::Histogram::Integral#differentiate
--- GSL::Histogram::Integral#diff
=== Graphics
--- GSL::Histogram#graph(options)
    This method uses the GNU plotutils (({graph})) to draw the histogram ((|self|)).
    The options as "-T X -C -l x" etc are given by a String.

=== Fitting
--- GSL::Histogram#fit_exponential(binstart = 0, binend = n-1)
    This method fits the histogram ((|self|)) to an exponential model 
    (({h[n] = a exp(b x[n])})) using  the bins of indices from ((|binstart|)) 
    to ((|binend|)). The result is returned as an Array of 6 elements,
    ((|[a, b, erra, errb, sumsq, dof]|)), where
    * ((|a|)): scale factor
    * ((|b|)): exponent
    * ((|erra, errb|)): fitting errors
    * ((|sumsq|)): fitting chi-squared (not reduced chi-squared)
    * ((|dof|)): degree-of-freedom, the number of bins used minus the number of parameters (2)

--- GSL::Histogram#fit_power(binstart = 0, binend = n-1)
    This method fits the histogram ((|self|)) to a power-law model 
    (({h[n] = a x[n]^b})) using  the bins of indices from ((|binstart|)) 
    to ((|binend|)). The result is returned as an Array of 6 elements,
    ((|[a, b, erra, errb, sumsq, dof]|)).

--- GSL::Histogram#fit_gaussian(binstart = 0, binend = n-1)
    This method fits the histogram ((|self|)) to Gaussian distribution
    using the bins of indices from ((|binstart|)) to ((|binend|)),
    and returns an Array of 8 elements, 
    ((|[sigma, mean, height, errsig, errmean, errhei, sumsq, dof]|)). 

    Example:
      #!/usr/bin/env ruby
      require("gsl")

      N = 10000
      MAX = 8
      rng = Rng.alloc

      data = Ran.gaussian(rng, 1.5, N) + 2
      h = Histogram.alloc(100, [-MAX, MAX])
      h.increment(data)

      sigma, mean, height, = h.fit_gaussian
      x = Vector.linspace(-MAX, MAX, 100)
      y = height*Ran::gaussian_pdf(x-mean, sigma)
      GSL::graph(h, [x, y], "-T X -C -g 3")

== The histogram probability distribution
The probability distribution function for a histogram consists of a set of bins 
which measure the probability of an event falling into a given range of 
a continuous variable x. A probability distribution function is defined
by the following class, which actually stores the cumulative probability 
distribution function. This is the natural quantity for generating samples 
via the inverse transform method, because there is a one-to-one mapping 
between the cumulative probability distribution and the range [0,1]. 
It can be shown that by taking a uniform random number in this 
range and finding its corresponding coordinate in the cumulative probability 
distribution we obtain samples with the desired probability distribution.

=== GSL::Histogram::Pdf class
--- GSL::Histogram::Pdf.alloc(n)
--- GSL::Histogram::Pdf.alloc(h)
    Constructors. If a histogram ((|h|)) is given, 
    the probability distribution is initialized with the contents of ((|h|)).

--- GSL::Histogram::Pdf#init(h)
    This initializes the probability distribution ((|self|)) with the contents 
    of the histogram ((|h|)). 

--- GSL::Histogram::Pdf#sample(r)
    This method uses ((|r|)), a uniform random number between zero and one, 
    to compute a single random sample from the probability distribution ((|self|)). 
    The algorithm used to compute the sample s is given by the following formula,
       s = range[i] + delta * (range[i+1] - range[i])
    where i is the index which satisfies 
    (({sum[i] <= r < sum[i+1]})) and (({delta})) is (({(r - sum[i])/(sum[i+1] - sum[i])})).

--- GSL::Histogram::Pdf#n
    This returns the number of bins of the probability distribution function.

--- GSL::Histogram:Pdf#range
    This returns a (({Vector::View})) object as a reference to the pointer 
    (({double *range})) in the (({gsl_histogram_pdf})) struct.
--- GSL::Histogram:Pdf#sum
    This returns a (({Vector::View})) object as a reference to the pointer 
    (({double *sum})) in the (({gsl_histogram_pdf})) struct.

((<prev|URL:stats.html>))
((<next|URL:hist2d.html>))

((<Reference index|URL:ref.html>))
((<top|URL:index.html>))

=end
