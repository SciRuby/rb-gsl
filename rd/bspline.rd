=begin
= Basis Splines
This chapter describes functions for the computation of smoothing basis splines (B-splines). This is only for GSL-1.9 or later.

(1) ((<Overview|URL:bspline.html#1>))
(2) ((<Initializing the B-splines solver|URL:bspline.html#2>))
(3) ((<Constructing the knots vector|URL:bspline.html#3>))
(4) ((<Evaluation of B-splines|URL:bspline.html#4>))

== Overview

B-splines are commonly used as basis functions to fit smoothing curves to large data sets. To do this, the abscissa axis is broken up into some number of intervals, where the endpoints of each interval are called breakpoints. These breakpoints are then converted to knots by imposing various continuity and smoothness conditions at each interface. Given a nondecreasing knot vector t = \{t_0, t_1, \dots, t_{n+k-1\, the n basis splines of order k are defined by for i = 0, \dots, n-1. The common case of cubic B-splines is given by k = 4. The above recurrence relation can be evaluated in a numerically stable way by the de Boor algorithm. 

If we define appropriate knots on an interval [a,b] then the B-spline basis functions form a complete set on that interval. Therefore we can expand a smoothing function as given enough (x_j, f(x_j)) data pairs. The c_i can be readily obtained from a least-squares fit. 

== Initializing the B-splines solver 
--- GSL::BSpline.alloc(k, nbreak)
    This method creates a workspace for computing B-splines of order ((|k|)). The number of breakpoints is given by ((|nbreak|)). This leads to ((|n = nbreak + k - 2|)) basis functions. Cubic B-splines are specified by ((|k = 4|)). The size of the workspace is ((|O(5k + nbreak)|)). 

== Constructing the knots vector 
--- GSL::BSpline#knots(breakpts)
    This method computes the knots associated with the given breakpoints ((|breakpts|)) and returns the knots as a (({GSL::Vector::View})) object.
--- GSL::BSpline#knots_uniform(a, b)
    This method assumes uniformly spaced breakpoints on [((|a,b|))] and constructs the corresponding knot vector using the previously specified ((|nbreak|)) parameter. 
== Evaluation of B-splines 
--- GSL::BSpline#eval(x[, B])
    This method evaluates all B-spline basis functions at the position ((|x|)) and stores them in ((|B|)) (if given), so that the ith element of ((|B|)) is ((|B_i(x)|)). ((|B|)) must be of length ((|n = nbreak + k - 2|)). If ((|B|)) is not given, a newly created vector is returned.It is far more efficient to compute all of the basis functions at once than to compute them individually, due to the nature of the defining recurrence relation. 

((<prev|URL:nonlinearfit.html>))
((<next|URL:const.html>))

((<Reference index|URL:ref.html>))
((<top|URL:index.html>))
=end
