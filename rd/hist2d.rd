=begin
= GSL::Histogram2d class

== Class methods
--- GSL::Histogram2d.alloc(nx, ny)
--- GSL::Histogram2d.alloc(xrange, yrange)
--- GSL::Histogram2d.alloc(nx, xmin, xmax, ny, ymin, ymax)
--- GSL::Histogram2d.alloc(nx, [xmin, xmax], ny, [ymin, ymax])
    Constructors.

    (1) With sizes
          irb(main):002:0> h2 = GSL::Histogram2d.alloc(2, 3)  # size 6
          irb(main):003:0> h2.xrange
          => GSL::Histogram::Range: 
          [ 0.000e+00 1.000e+00 2.000e+00 ]
          irb(main):004:0> h2.yrange
          => GSL::Histogram::Range: 
          [ 0.000e+00 1.000e+00 2.000e+00 3.000e+00 ]
          irb(main):005:0> h2.bin
          => GSL::Histogram::Bin: 
          [ 0.000e+00 0.000e+00 0.000e+00 0.000e+00 0.000e+00 0.000e+00 ]
          irb(main):006:0> h2.xrange.class.superclass
          => GSL::Vector::View::ReadOnly
          irb(main):007:0> h2.bin.class.superclass
          => GSL::Vector::View::ReadOnly

    (2) With ranges
           irb(main):006:0> h2 = Histogram2d.alloc([3, 6, 7], [2, 3, 8])
           irb(main):007:0> h2.xrange
           => GSL::Histogram::Range: 
           [ 3.000e+00 6.000e+00 7.000e+00 ]
           irb(main):008:0> h2.yrange
           => GSL::Histogram::Range: 
           [ 2.000e+00 3.000e+00 8.000e+00 ]
           irb(main):009:0> h2.bin
           => GSL::Histogram::Bin: 
           [ 0.000e+00 0.000e+00 0.000e+00 0.000e+00 ]

    (3) With sizes and ranges
           irb(main):010:0> h2 = Histogram2d.alloc(4, [0, 4], 3, [1, 5])
           irb(main):011:0> h2.xrange
           => GSL::Histogram::Range: 
           [ 0.000e+00 1.000e+00 2.000e+00 3.000e+00 4.000e+00 ]
           irb(main):012:0> h2.yrange
           => GSL::Histogram::Range: 
           [ 1.000e+00 2.333e+00 3.667e+00 5.000e+00 ]

--- GSL::Histogram2d.alloc_uniform(nx, xmin, xmax, ny, ymin, ymax)
--- GSL::Histogram2d.equal_bins_p(h1, h2)

== Methods
--- GSL::Histogram2d#set_ranges(vx, vy)
    This sets the ranges of the existing histogram using 
    two ((<GSL::Vector|URL:vector.html>)) objects or arrays.
--- GSL::Histogram2d#set_ranges_uniform(xmin, xmax, ymin, ymax)
--- GSL::Histogram2d#set_ranges_uniform([xmin, xmax], [ymin, ymax])
    This sets the ranges of the existing histogram ((|self|)) to cover 
    the ranges ((|xmin|)) to ((|xmax|)) and ((|ymin|)) to ((|ymax|)) uniformly. 
    The values of the histogram bins are reset to zero.

--- GSL::Histogram2d#clone
--- GSL::Histogram2d#duplicate
    Returns a newly created histogram which is an exact copy of the 
    histogram ((|self|)).

--- GSL::Histogram2d#increment(x, y, weight = 1)
--- GSL::Histogram2d#accumulate(x, y, weight = 1)
--- GSL::Histogram2d#fill(x, y, weight = 1)
    These method update the histogram ((|self|)) by adding ((|weight|)) 
    to the bin whose ((|x|)) and ((|y|)) ranges contain the coordinates (x,y).
    If (x,y) lies outside the limits of the histogram then none of the 
    bins are modified. 

--- GSL::Histogram2d#increment2(x, y, weight = 1)
--- GSL::Histogram2d#accumulate2(x, y, weight = 1)
--- GSL::Histogram2d#fill2(x, y, weight = 1)
    These are similar to (({GSL::Histogram#increment/accumulate/fill})),
    but when (x,y) lies outside the limits of the histogram then the edge
    bin is incremented.

--- GSL::Histogram2d#get(i[, j])
--- GSL::Histogram2d#get([i, j])
--- GSL::Histogram2d#[i, j]
    Return the contents of the ((|(i,j)|))-th bin of the histogram ((|self|)). 

--- GSL::Histogram2d#get_xrange(i)
--- GSL::Histogram2d#get_yrange(j)
    These methods find the upper and lower range limits of the ((|i|))-th 
    and ((|j|))-th bins in the x and y directions of the histogram ((|self|)).

    Ex: 
       irb(main):030:0> h2 = Histogram2d.alloc(2, [0, 2], 3, [1, 4])
       irb(main):031:0> h2.get_xrange(1)
       => [1.0, 2.0]
       irb(main):032:0> h2.get_yrange(1)
       => [2.0, 3.0]

--- GSL::Histogram2d#xrange
--- GSL::Histogram2d#yrange
    These methods returns the range of histogram as (({Vector::View})) objects.

--- GSL::Histogram2d#xmax
--- GSL::Histogram2d#xmin
--- GSL::Histogram2d#ymax
--- GSL::Histogram2d#ymin
--- GSL::Histogram2d#nx
--- GSL::Histogram2d#ny
    These methodss return the maximum upper and minimum lower range limits 
    and the number of bins for the x and y directions of the histogram ((|self|)).

--- GSL::Histogram2d#reset
    Resets all the bins of the histogram ((|self|)) to zero.

--- GSL::Histogram2d#find(x, y)
    Finds the indices i and j to the to the bin which covers the 
    coordinates (((|x,y|))).

--- GSL::Histogram2d#max_val
    Returns the maximum value contained in the histogram bins.

--- GSL::Histogram2d#max_bin
    Returns the indices (i,j) of the bin containing the maximum value in 
    the histogram ((|self|)).

--- GSL::Histogram2d#min_val
    Returns the minimum value contained in the histogram bins.
--- GSL::Histogram2d#min_bin
    Returns the indices (i,j) of the bin containing the minimum value
    in the histogram ((|self|)).

--- GSL::Histogram2d#xmean
    Returns the mean of the histogrammed x variable, where the histogram 
    is regarded as a probability distribution. 
    Negative bin values are ignored for the purposes of this calculation.

--- GSL::Histogram2d#xsigma
    Returns the standard deviation of the histogrammed x variable, 
    where the histogram is regarded as a probability distribution. 
    Negative bin values are ignored for the purposes of this calculation.

--- GSL::Histogram2d#ymean
    Returns the mean of the histogrammed y variable, where the histogram 
    is regarded as a probability distribution. Negative bin values are 
    ignored for the purposes of this calculation.

--- GSL::Histogram2d#ysigma
    Returns the standard deviation of the histogrammed y variable, 
    where the histogram is regarded as a probability distribution. 
    Negative bin values are ignored for the purposes of this calculation.

--- GSL::Histogram2d#cov
    Returns the covariance of the histogrammed x and y variables, 
    where the histogram is regarded as a probability distribution. 
    Negative bin values are ignored for the purposes of this calculation.

--- GSL::Histogram2d#sum
--- GSL::Histogram2d#integral
    Return the sum of all bin values. Negative bin values are included in the sum.
--- GSL::Histogram2d#add
--- GSL::Histogram2d#sub
--- GSL::Histogram2d#mul
--- GSL::Histogram2d#div
--- GSL::Histogram2d#scale(val)
--- GSL::Histogram2d#shift(val)
--- GSL::Histogram2d#fwrite(io)
--- GSL::Histogram2d#fwrite(filename)
--- GSL::Histogram2d#fread(io)
--- GSL::Histogram2d#fread(filename)
--- GSL::Histogram2d#fprintf(io, range_format, bin_format)
--- GSL::Histogram2d#fprintf(filename, range_format, bin_format)
--- GSL::Histogram2d#fscanf(io)
--- GSL::Histogram2d#fscanf(filename)

--- GSL::Histogram2d#xproject(jstart = 0, jend = ny-1)
    Creates a (({GSL::Histogram})) object by projecting the 2D histogram 
    ((|self|)) onto the x-axis over the y-range from ((|jstart|)) to ((|jend|)).

--- GSL::Histogram2d#yproject(istart = 0, iend = nx-1)
    Creates a (({GSL::Histogram})) object by projecting the 2D histogram 
    ((|self|)) onto the y-axis over the x-range from ((|istart|)) to ((|iend|)).

== GSL::Histogram2d::Pdf class
--- GSL::Histogram2d::Pdf.alloc(nx, ny)
    Constructors

--- GSL::Histogram2d::Pdf#init(h)
    Initializer with a GSL::Histogram2d object ((|h|)).

--- GSL::Histogram2d::Pdf#sample(r1, r2)

== Example
  #!/usr/bin/env ruby
  require("gsl")

  N = 10000
  BINS = 100

  rng = Rng.new(2)
  h2 = Histogram2d.alloc(BINS, [-5, 5], BINS, [-8, 8])

  sig1 = 1.0
  sig2 = 2.0

  for i in 0...N do
    r1 = rng.gaussian(sig1)
    r2 = rng.gaussian(sig2)
    h2.increment(r1, r2)
  end

  hx = h2.xproject
  hy = h2.yproject
  printf("%f %f %f %f\n", h2.xsigma, h2.ysigma, hx.sigma, hy.sigma)
  GSL::graph(hx, hy, "-T X -C -g 3")


((<prev|URL:hist.html>))
((<next|URL:hist3d.html>))

((<Reference index|URL:ref.html>))
((<top|URL:index.html>))
=end
